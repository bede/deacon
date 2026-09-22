use crate::minimizers::{
    Buffers, KmerHasher, decode_u64, decode_u128, extend_minimizers_unchecked,
};
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
    /// Distinct index hits (a sufficient lower bound without diagnostics)
    pub hit_count_lower_bound: usize,
    /// Distinct minimizers (a positional upper bound when the exact count is unnecessary)
    pub minimizer_count_upper_bound: usize,
    /// Canonical hit k-mers, populated only with diagnostics
    pub hit_kmers: Vec<String>,
}

/// Reused per record for hit dedup, then for `count_distinct_minimizers`
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

    /// Min hits needed to hit rel threshold. Never decreases as `minimizer_count` grows
    #[inline]
    fn rel_required_hits(&self, minimizer_count: usize) -> usize {
        if minimizer_count == 0 {
            0
        } else {
            let lower = (self.params.rel_threshold * minimizer_count as f64) as usize;
            (lower
                + usize::from(lower as f64 / (minimizer_count as f64) < self.params.rel_threshold))
            .max(1)
        }
    }

    #[inline]
    fn required_hits(&self, minimizer_count: usize) -> usize {
        self.params
            .abs_threshold
            .max(self.rel_required_hits(minimizer_count))
    }

    /// True if hits are between the absolute and relative floors, so positions can't decide
    #[inline]
    fn needs_exact_minimizer_count(
        &self,
        hit_count_lower_bound: usize,
        positional_minimizer_count: usize,
    ) -> bool {
        hit_count_lower_bound >= self.params.abs_threshold
            && hit_count_lower_bound < self.rel_required_hits(positional_minimizer_count)
    }

    #[inline]
    fn is_index_match(
        &self,
        hit_count_lower_bound: usize,
        minimizer_count_upper_bound: usize,
    ) -> bool {
        minimizer_count_upper_bound > 0
            && hit_count_lower_bound >= self.required_hits(minimizer_count_upper_bound)
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

    /// Classify with early exit and possibly bounded counts
    pub fn classify_read(&mut self, index: &MinimizerSet, seq: &[u8]) -> FilterDecision {
        self.classify_seqs(index, &[seq], false)
    }

    /// Classify with exact counts and hit k-mers
    pub fn classify_read_with_diagnostics(
        &mut self,
        index: &MinimizerSet,
        seq: &[u8],
    ) -> FilterDecision {
        self.classify_seqs(index, &[seq], true)
    }

    /// Classify a pair with early exit and possibly bounded counts
    pub fn classify_pair(
        &mut self,
        index: &MinimizerSet,
        seq1: &[u8],
        seq2: &[u8],
    ) -> FilterDecision {
        self.classify_seqs(index, &[seq1, seq2], false)
    }

    /// Classify a pair with exact counts and hit k-mers
    pub fn classify_pair_with_diagnostics(
        &mut self,
        index: &MinimizerSet,
        seq1: &[u8],
        seq2: &[u8],
    ) -> FilterDecision {
        self.classify_seqs(index, &[seq1, seq2], true)
    }

    /// Pool distinct minimizer hits across mates
    fn classify_seqs(
        &mut self,
        index: &MinimizerSet,
        seqs: &[&[u8]],
        exact: bool,
    ) -> FilterDecision {
        self.seen_minimizers.clear();
        // Records with no long enough mate skip the loop below, so clear here or they
        // inherit the previous record's minimizers
        self.buffers.minimizers.clear();
        let mut hit_kmers = Vec::new();

        // Pool mates first: a per-mate `stop_at` would understate the denominator
        for &seq in seqs {
            let seq = self.seq_for_filter(seq);
            if seq.len() < self.kmer_length as usize {
                continue;
            }
            extend_minimizers_unchecked(
                seq,
                &self.hasher,
                self.kmer_length,
                self.window_size,
                &mut self.buffers,
            );
        }

        // Positions bound the distinct count and required_hits never decreases, so no
        // denominator can ask for more hits than this: stop probing once it's reached
        let positional_minimizer_count = self.buffers.minimizers.len();
        let stop_at = if exact {
            usize::MAX
        } else {
            self.required_hits(positional_minimizer_count)
        };
        let hit_count_lower_bound = self.count_hits_until(index, exact, stop_at, &mut hit_kmers);

        // Avoid costly distinct counting where possible
        let minimizer_count_upper_bound = if exact
            || self.needs_exact_minimizer_count(hit_count_lower_bound, positional_minimizer_count)
        {
            self.count_distinct_minimizers()
        } else {
            positional_minimizer_count
        };

        let is_match = self.is_index_match(hit_count_lower_bound, minimizer_count_upper_bound);
        let keep = is_match != self.params.deplete;

        FilterDecision {
            keep,
            hit_count_lower_bound,
            minimizer_count_upper_bound,
            hit_kmers,
        }
    }

    /// Exact distinct minimizer count pooled across mates. Expensive but should be rare
    fn count_distinct_minimizers(&mut self) -> usize {
        self.seen_minimizers.clear();
        self.insert_buffer_minimizers();
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

    /// Counts distinct hits, stopping at `stop_at` where the match is already decided
    fn count_hits_until(
        &mut self,
        index: &MinimizerSet,
        exact: bool,
        stop_at: usize,
        hit_kmers: &mut Vec<String>,
    ) -> usize {
        let mut distinct_hits = 0;
        match (&self.buffers.minimizers, index, &mut self.seen_minimizers) {
            (MinimizerVec::U64(vec), MinimizerSet::U64(set), SeenMinimizers::U64(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) {
                        distinct_hits += 1;
                        if exact {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u64(minimizer, self.kmer_length))
                                    .to_string(),
                            );
                        }
                        if distinct_hits >= stop_at {
                            break;
                        }
                    }
                }
            }
            (MinimizerVec::U64(vec), MinimizerSet::Fuse(filter), SeenMinimizers::U64(seen)) => {
                for &minimizer in vec {
                    if filter.contains(minimizer) && seen.insert(minimizer) {
                        distinct_hits += 1;
                        if exact {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u64(minimizer, self.kmer_length))
                                    .to_string(),
                            );
                        }
                        if distinct_hits >= stop_at {
                            break;
                        }
                    }
                }
            }
            (MinimizerVec::U128(vec), MinimizerSet::U128(set), SeenMinimizers::U128(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) {
                        distinct_hits += 1;
                        if exact {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u128(minimizer, self.kmer_length))
                                    .to_string(),
                            );
                        }
                        if distinct_hits >= stop_at {
                            break;
                        }
                    }
                }
            }
            _ => unreachable!("minimizer width does not match index variant"),
        }
        distinct_hits
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
        assert!(kernel.classify_read(&index, b"ACGT").keep);

        let mut kernel = FilterKernel::new(
            31,
            15,
            FilterParams {
                deplete: false,
                ..params
            },
        )
        .unwrap();
        assert!(!kernel.classify_read(&index, b"ACGT").keep);
    }

    #[test]
    fn zero_minimizers_never_match() {
        let index = MinimizerSet::U64(RapidHashSet::default());
        for deplete in [false, true] {
            let mut kernel = FilterKernel::new(
                31,
                15,
                FilterParams {
                    deplete,
                    abs_threshold: 0,
                    rel_threshold: 0.0,
                    prefix_length: 0,
                },
            )
            .unwrap();
            assert_eq!(kernel.classify_read(&index, b"").keep, deplete);
        }
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
        assert!(!kernel.is_index_match(2, 5));
        assert!(kernel.is_index_match(3, 5));
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

    /// Neither using positions instead of the distinct count nor exiting early once the
    /// match is certain may change outcome. `exact` forces the exact path, and the repeat
    /// read makes sure the substitution actually happens.
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
        let mut early_exit_observed = false;
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
                        let fast = kernel.classify_read(&index, read);
                        let exact = kernel.classify_read_with_diagnostics(&index, read);
                        assert_eq!(fast.keep, exact.keep, "read {i}, {context}");
                        assert_hits_sufficient(
                            &kernel,
                            &fast,
                            &exact,
                            &format!("read {i}, {context}"),
                        );
                        substitution_observed |=
                            fast.minimizer_count_upper_bound > exact.minimizer_count_upper_bound;
                        early_exit_observed |=
                            fast.hit_count_lower_bound < exact.hit_count_lower_bound;
                    }

                    for (i, first) in reads.iter().enumerate() {
                        for (j, second) in reads.iter().enumerate() {
                            let fast = kernel.classify_pair(&index, first, second);
                            let exact =
                                kernel.classify_pair_with_diagnostics(&index, first, second);
                            assert_eq!(fast.keep, exact.keep, "pair {i}/{j}, {context}");
                            assert_hits_sufficient(
                                &kernel,
                                &fast,
                                &exact,
                                &format!("pair {i}/{j}, {context}"),
                            );
                            substitution_observed |= fast.minimizer_count_upper_bound
                                > exact.minimizer_count_upper_bound;
                            early_exit_observed |=
                                fast.hit_count_lower_bound < exact.hit_count_lower_bound;
                        }
                    }
                }
            }
        }
        assert!(
            substitution_observed,
            "positional count never stood in for the distinct count, so nothing was tested"
        );
        assert!(
            early_exit_observed,
            "hit counting never exited early, so nothing was tested"
        );
    }

    /// The fast path may stop counting hits early, but only once it has enough of them to
    /// match under the exact distinct denominator, and it must never invent hits
    fn assert_hits_sufficient(
        kernel: &FilterKernel,
        fast: &FilterDecision,
        exact: &FilterDecision,
        context: &str,
    ) {
        assert!(
            fast.hit_count_lower_bound <= exact.hit_count_lower_bound,
            "{context}"
        );
        assert!(
            fast.hit_count_lower_bound == exact.hit_count_lower_bound
                || fast.hit_count_lower_bound
                    >= kernel.required_hits(exact.minimizer_count_upper_bound),
            "{context}"
        );
    }

    /// A record without minimizers must not inherit the previous record's pooled buffer
    #[test]
    fn short_record_after_long_record_has_no_hits() {
        const K: u8 = 7;
        const W: u8 = 3;
        let reference: &[u8] = b"ACGTTGCAAGGCTTAACCGGTTACGATCGATCGGATCCTAGCTAGCTTAACCGGATCGTA";
        let index = index_of(reference, K, W);
        let mut kernel = FilterKernel::new(
            K,
            W,
            FilterParams {
                deplete: true,
                abs_threshold: 1,
                rel_threshold: 0.0,
                prefix_length: 0,
            },
        )
        .unwrap();

        assert!(!kernel.classify_read(&index, reference).keep);
        let decision = kernel.classify_read(&index, b"ACG");
        assert_eq!(decision.hit_count_lower_bound, 0);
        assert_eq!(decision.minimizer_count_upper_bound, 0);
        assert!(decision.keep);

        assert!(!kernel.classify_pair(&index, reference, reference).keep);
        let decision = kernel.classify_pair(&index, b"ACG", b"");
        assert_eq!(decision.hit_count_lower_bound, 0);
        assert_eq!(decision.minimizer_count_upper_bound, 0);
        assert!(decision.keep);
    }

    /// `needs_exact_minimizer_count` relies on `rel_required_hits` never decreasing
    #[test]
    fn rel_required_hits_never_decreases() {
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
            for total in 0..20000 {
                let required = kernel.rel_required_hits(total);
                assert!(
                    required >= previous,
                    "rel_required_hits decreased at total={total} (rel={rel_threshold})"
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

        let decision = kernel.classify_pair(&index, b"AAAAA", b"AAAAA");
        assert_eq!(decision.hit_count_lower_bound, 1);
        assert_eq!(decision.minimizer_count_upper_bound, 1);
        assert!(decision.keep);
    }
}
