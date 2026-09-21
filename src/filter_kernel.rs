use crate::minimizers::{Buffers, KmerHasher, decode_u64, decode_u128, fill_minimizers_unchecked};
use crate::{MinimizerSet, MinimizerVec, RapidHashSet, validate_unit_interval};

#[derive(Clone, Copy, Debug)]
pub struct FilterParams {
    pub deplete: bool,
    pub abs_threshold: usize,
    pub rel_threshold: f64,
    pub prefix_length: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FilterDecision {
    pub keep: bool,
    pub hit_count: usize,
    /// Distinct minimizer count, or an upper-bound positional count if the relative
    /// threshold wasn't evaluated. debug forces distinct calculation
    pub total_minimizers: usize,
    pub hit_kmers: Vec<String>,
}

/// Reused per record: holds matching minimizers while counting hits, then every minimizer
/// if `distinct_minimizers` has to compute an exact denominator.
#[derive(Clone)]
enum SeenMinimizers {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl SeenMinimizers {
    fn new(kmer_length: u8) -> Self {
        if kmer_length <= 32 {
            Self::U64(RapidHashSet::default())
        } else {
            Self::U128(RapidHashSet::default())
        }
    }

    fn clear(&mut self) {
        match self {
            Self::U64(set) => set.clear(),
            Self::U128(set) => set.clear(),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::U64(set) => set.len(),
            Self::U128(set) => set.len(),
        }
    }
}

#[derive(Clone)]
pub struct FilterKernel {
    params: FilterParams,
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    buffers: Buffers,
    seen_minimizers: SeenMinimizers,
}

impl FilterKernel {
    pub fn new(kmer_length: u8, window_size: u8, params: FilterParams) -> anyhow::Result<Self> {
        validate_unit_interval("relative threshold", params.rel_threshold)?;
        Ok(Self {
            params,
            kmer_length,
            window_size,
            hasher: KmerHasher::new(kmer_length as usize),
            buffers: if kmer_length <= 32 {
                Buffers::new_u64()
            } else {
                Buffers::new_u128()
            },
            seen_minimizers: SeenMinimizers::new(kmer_length),
        })
    }

    #[inline]
    pub fn params(&self) -> FilterParams {
        self.params
    }

    /// Min hits needed to hit rel threshold. Never decreases as `total_minimizers` grows,
    /// so `classify_seqs` can safely substitute the positional count for the distinct count.
    #[inline]
    fn rel_required_hits(&self, total_minimizers: usize) -> usize {
        if total_minimizers == 0 {
            0
        } else {
            let lower = (self.params.rel_threshold * total_minimizers as f64) as usize;
            (lower
                + usize::from(lower as f64 / (total_minimizers as f64) < self.params.rel_threshold))
            .max(1)
        }
    }

    #[inline]
    fn required_hits(&self, total_minimizers: usize) -> usize {
        self.params
            .abs_threshold
            .max(self.rel_required_hits(total_minimizers))
    }

    /// True if hit_count is between the absolute and relative floors, where swapping in
    /// the positional count for the distinct count could change outcome.
    #[inline]
    fn needs_distinct(&self, hit_count: usize, positions: usize) -> bool {
        hit_count >= self.params.abs_threshold && hit_count < self.rel_required_hits(positions)
    }

    #[inline]
    fn keep_from_counts(&self, hit_count: usize, total_minimizers: usize) -> bool {
        let required = self.required_hits(total_minimizers);
        if self.params.deplete {
            hit_count < required
        } else {
            hit_count >= required
        }
    }

    #[inline]
    fn seq_for_filter<'a>(&self, seq: &'a [u8]) -> &'a [u8] {
        let seq = if self.params.prefix_length > 0 && seq.len() > self.params.prefix_length {
            &seq[..self.params.prefix_length]
        } else {
            seq
        };
        seq.strip_suffix(b"\n").unwrap_or(seq)
    }

    pub fn classify_read(
        &mut self,
        index: &MinimizerSet,
        seq: &[u8],
        debug: bool,
    ) -> FilterDecision {
        self.classify_seqs(index, &[seq], debug)
    }

    pub fn classify_pair(
        &mut self,
        index: &MinimizerSet,
        seq1: &[u8],
        seq2: &[u8],
        debug: bool,
    ) -> FilterDecision {
        self.classify_seqs(index, &[seq1, seq2], debug)
    }

    /// Pool distinct minimizer hits across mates
    fn classify_seqs(
        &mut self,
        index: &MinimizerSet,
        seqs: &[&[u8]],
        debug: bool,
    ) -> FilterDecision {
        self.seen_minimizers.clear();
        let mut hit_count = 0;
        let mut positions = 0;
        let mut filled = 0;
        let mut hit_kmers = Vec::new();

        for &seq in seqs {
            let seq = self.seq_for_filter(seq);
            if seq.len() < self.kmer_length as usize {
                continue;
            }
            fill_minimizers_unchecked(
                seq,
                &self.hasher,
                self.kmer_length,
                self.window_size,
                &mut self.buffers,
            );
            positions += self.buffers.minimizers.len();
            filled += 1;
            hit_count += self.count_buffer_hits(index, debug, &mut hit_kmers);
        }

        // Avoid costly distinct counting where poss
        let total_minimizers = if debug || self.needs_distinct(hit_count, positions) {
            self.distinct_minimizers(seqs, filled)
        } else {
            positions
        };

        let keep = if total_minimizers == 0 {
            // No minis (too short / all-ambiguous) means no possible index match: keep only when depleting
            self.params.deplete
        } else {
            self.keep_from_counts(hit_count, total_minimizers)
        };

        FilterDecision {
            keep,
            hit_count,
            total_minimizers,
            hit_kmers,
        }
    }

    /// Exact distinct minimizer count pooled across mates. Expensive but should be rare
    fn distinct_minimizers(&mut self, seqs: &[&[u8]], filled: usize) -> usize {
        self.seen_minimizers.clear();
        if filled == 1 {
            // Exactly one sequence produced minimizers, so the buffer still holds all of them
            self.insert_buffer_minimizers();
        } else {
            for &seq in seqs {
                let seq = self.seq_for_filter(seq);
                if seq.len() < self.kmer_length as usize {
                    continue;
                }
                fill_minimizers_unchecked(
                    seq,
                    &self.hasher,
                    self.kmer_length,
                    self.window_size,
                    &mut self.buffers,
                );
                self.insert_buffer_minimizers();
            }
        }
        self.seen_minimizers.len()
    }

    fn insert_buffer_minimizers(&mut self) {
        match (&self.buffers.minimizers, &mut self.seen_minimizers) {
            (MinimizerVec::U64(vec), SeenMinimizers::U64(seen)) => {
                for &minimizer in vec {
                    seen.insert(minimizer);
                }
            }
            (MinimizerVec::U128(vec), SeenMinimizers::U128(seen)) => {
                for &minimizer in vec {
                    seen.insert(minimizer);
                }
            }
            _ => unreachable!("minimizer width does not match seen set variant"),
        }
    }

    fn count_buffer_hits(
        &mut self,
        index: &MinimizerSet,
        debug: bool,
        hit_kmers: &mut Vec<String>,
    ) -> usize {
        let mut hit_count = 0;
        match (&self.buffers.minimizers, index, &mut self.seen_minimizers) {
            (MinimizerVec::U64(vec), MinimizerSet::U64(set), SeenMinimizers::U64(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) {
                        hit_count += 1;
                        if debug {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u64(minimizer, self.kmer_length))
                                    .to_string(),
                            );
                        }
                    }
                }
            }
            (MinimizerVec::U64(vec), MinimizerSet::Fuse(filter), SeenMinimizers::U64(seen)) => {
                for &minimizer in vec {
                    if filter.contains(minimizer) && seen.insert(minimizer) {
                        hit_count += 1;
                        if debug {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u64(minimizer, self.kmer_length))
                                    .to_string(),
                            );
                        }
                    }
                }
            }
            (MinimizerVec::U128(vec), MinimizerSet::U128(set), SeenMinimizers::U128(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) {
                        hit_count += 1;
                        if debug {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u128(minimizer, self.kmer_length))
                                    .to_string(),
                            );
                        }
                    }
                }
            }
            _ => unreachable!("minimizer width does not match index variant"),
        }
        hit_count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RapidHashSet;

    #[test]
    fn short_read_behavior_matches_filter_modes() {
        let index = MinimizerSet::U64(RapidHashSet::default());
        let params = FilterParams {
            deplete: true,
            abs_threshold: 2,
            rel_threshold: 0.01,
            prefix_length: 0,
        };
        let mut kernel = FilterKernel::new(31, 15, params).unwrap();
        assert!(kernel.classify_read(&index, b"ACGT", false).keep);

        let mut kernel = FilterKernel::new(
            31,
            15,
            FilterParams {
                deplete: false,
                ..params
            },
        )
        .unwrap();
        assert!(!kernel.classify_read(&index, b"ACGT", false).keep);
    }

    #[test]
    fn relative_threshold_is_a_minimum_proportion() {
        let kernel = FilterKernel::new(
            31,
            15,
            FilterParams {
                deplete: false,
                abs_threshold: 1,
                rel_threshold: 0.49,
                prefix_length: 0,
            },
        )
        .unwrap();

        assert_eq!(kernel.required_hits(5), 3);
        assert!(!kernel.keep_from_counts(2, 5));
        assert!(kernel.keep_from_counts(3, 5));
    }

    #[test]
    fn relative_threshold_avoids_decimal_over_ceiling() {
        let kernel = FilterKernel::new(
            31,
            15,
            FilterParams {
                deplete: false,
                abs_threshold: 1,
                rel_threshold: 0.1,
                prefix_length: 0,
            },
        )
        .unwrap();

        assert_eq!(kernel.required_hits(30), 3);
    }

    /// Builds a minimizer index from `reference` for use in tests.
    fn index_of(reference: &[u8], kmer_length: u8, window_size: u8) -> MinimizerSet {
        let hasher = KmerHasher::new(kmer_length as usize);
        match crate::minimizers::compute_minimizers(reference, &hasher, kmer_length, window_size) {
            MinimizerVec::U64(values) => MinimizerSet::U64(values.into_iter().collect()),
            MinimizerVec::U128(values) => MinimizerSet::U128(values.into_iter().collect()),
        }
    }

    /// Using positions instead of distinct count must never change outcome. `debug`
    /// forces the real count so both paths agree, and the repeat read makes sure the
    /// substitution actually happens.
    #[test]
    fn positional_denominator_never_changes_the_verdict() {
        const K: u8 = 7;
        const W: u8 = 3;
        let reference: &[u8] = b"ACGTTGCAAGGCTTAACCGGTTACGATCGATCGGATCCTAGCTAGCTTAACCGGATCGTA";
        let index = index_of(reference, K, W);

        let reads: Vec<Vec<u8>> = vec![
            reference.to_vec(),                                        // in index
            reference[..20].repeat(8), // repeat: distinct << positions
            [&reference[..25], &b"TTTTTTTTTTTTTTTTTTTT"[..]].concat(), // partial match
            b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec(), // homopolymer
            b"ACG".to_vec(),           // shorter than k
            Vec::new(),                // empty
        ];

        let mut substitution_observed = false;
        for &deplete in &[false, true] {
            for &abs_threshold in &[1usize, 2, 3, 5] {
                for &rel_threshold in &[0.0, 0.01, 0.1, 0.25, 0.5, 0.9, 1.0] {
                    let mut kernel = FilterKernel::new(
                        K,
                        W,
                        FilterParams {
                            deplete,
                            abs_threshold,
                            rel_threshold,
                            prefix_length: 0,
                        },
                    )
                    .unwrap();
                    let context =
                        format!("deplete={deplete} abs={abs_threshold} rel={rel_threshold}");

                    for (i, read) in reads.iter().enumerate() {
                        let fast = kernel.classify_read(&index, read, false);
                        let exact = kernel.classify_read(&index, read, true);
                        assert_eq!(fast.keep, exact.keep, "read {i}, {context}");
                        assert_eq!(fast.hit_count, exact.hit_count, "read {i}, {context}");
                        substitution_observed |= fast.total_minimizers > exact.total_minimizers;
                    }

                    for (i, first) in reads.iter().enumerate() {
                        for (j, second) in reads.iter().enumerate() {
                            let fast = kernel.classify_pair(&index, first, second, false);
                            let exact = kernel.classify_pair(&index, first, second, true);
                            assert_eq!(fast.keep, exact.keep, "pair {i}/{j}, {context}");
                            assert_eq!(fast.hit_count, exact.hit_count, "pair {i}/{j}, {context}");
                            substitution_observed |= fast.total_minimizers > exact.total_minimizers;
                        }
                    }
                }
            }
        }
        assert!(
            substitution_observed,
            "positional count never stood in for the distinct count, so nothing was tested"
        );
    }

    /// `needs_distinct` relies on `rel_required_hits` being monotone -- check that here.
    #[test]
    fn rel_required_hits_is_monotone() {
        for &rel_threshold in &[0.0, 0.01, 0.1, 0.5, 0.99, 1.0] {
            let kernel = FilterKernel::new(
                31,
                15,
                FilterParams {
                    deplete: false,
                    abs_threshold: 1,
                    rel_threshold,
                    prefix_length: 0,
                },
            )
            .unwrap();
            let mut previous = 0;
            for total in 0..2000 {
                let required = kernel.rel_required_hits(total);
                assert!(
                    required >= previous,
                    "rel_required_hits not monotone at total={total} (rel={rel_threshold})"
                );
                previous = required;
            }
        }
    }

    #[test]
    fn relative_threshold_uses_distinct_minimizer_denominator() {
        let index = MinimizerSet::U64(RapidHashSet::from_iter([0]));
        let mut kernel = FilterKernel::new(
            3,
            1,
            FilterParams {
                deplete: false,
                abs_threshold: 1,
                rel_threshold: 1.0,
                prefix_length: 0,
            },
        )
        .unwrap();

        let decision = kernel.classify_pair(&index, b"AAAAA", b"AAAAA", false);
        assert_eq!(decision.hit_count, 1);
        assert_eq!(decision.total_minimizers, 1);
        assert!(decision.keep);
    }
}
