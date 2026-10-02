use crate::index::{IndexStorage, RapidHashSet};
use crate::minimizers::{decode_u64, decode_u128};
use crate::{Index, MinimizerVec, Minimizers, validate_unit_interval};
use std::sync::Arc;

/// Matching thresholds. Defaults to search, two hits and a 0.01 hit fraction.
#[derive(Clone, Copy, Debug)]
pub struct FilterParams {
    pub deplete: bool,
    pub abs_threshold: usize,
    pub rel_threshold: f64,
    pub prefix_length: usize,
}

impl Default for FilterParams {
    fn default() -> Self {
        Self {
            deplete: false,
            abs_threshold: 2,
            rel_threshold: 0.01,
            prefix_length: 0,
        }
    }
}

impl FilterParams {
    pub fn validate(&self) -> anyhow::Result<()> {
        anyhow::ensure!(self.abs_threshold > 0, "abs_threshold must be at least 1");
        validate_unit_interval("relative threshold", self.rel_threshold)
    }
}

/// Exact distinct counts. Fuse hits may include false positives.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FilterScore {
    pub is_match: bool,
    pub keep: bool,
    pub hit_count: usize,
    pub minimizer_count: usize,
}

/// Filtering decision and count bounds.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FilterDecision {
    /// Thresholds met, regardless of search/deplete mode.
    pub is_match: bool,
    pub keep: bool,
    /// Lower bound on distinct hits after early exit.
    pub hit_count_lower_bound: usize,
    /// Distinct minimizers, or a positional upper bound.
    pub minimizer_count_upper_bound: usize,
}

/// Exact score and decoded canonical hits. Fuse hits may include false positives.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FilterDiagnostics {
    pub score: FilterScore,
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

/// Filtering worker with a shared index and reusable buffers.
///
/// Use one worker per concurrent task. Calls create no threads. Input is borrowed
/// ASCII DNA, case-insensitive. Ambiguous windows are skipped. Empty, short and
/// entirely ambiguous reads never match. Pairs pool distinct minimizers, applying
/// the prefix limit to each mate.
#[derive(Clone)]
pub struct FilterKernel {
    index: Arc<Index>,
    params: FilterParams,
    minimizers: Minimizers,
    seen_minimizers: SeenMinimizers,
}

impl FilterKernel {
    pub fn new(index: Arc<Index>, params: FilterParams) -> anyhow::Result<Self> {
        params.validate()?;
        let kmer_length = index.kmer_length();
        Ok(Self {
            minimizers: Minimizers::new(kmer_length, index.window_size())?,
            seen_minimizers: SeenMinimizers::new(kmer_length),
            index,
            params,
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
    pub fn classify_read(&mut self, seq: &[u8]) -> FilterDecision {
        self.classify_seqs::<false, false>(&[seq]).0
    }

    /// Classify with exact counts and hit k-mers
    pub fn classify_read_with_diagnostics(&mut self, seq: &[u8]) -> FilterDiagnostics {
        self.diagnostics(&[seq])
    }

    /// Classify a pair with early exit and possibly bounded counts
    pub fn classify_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> FilterDecision {
        self.classify_seqs::<false, false>(&[seq1, seq2]).0
    }

    /// Classify a pair with exact counts and hit k-mers
    pub fn classify_pair_with_diagnostics(
        &mut self,
        seq1: &[u8],
        seq2: &[u8],
    ) -> FilterDiagnostics {
        self.diagnostics(&[seq1, seq2])
    }

    /// Count distinct hits and minimizers without decoding.
    pub fn score_read(&mut self, seq: &[u8]) -> FilterScore {
        Self::score(self.classify_seqs::<true, false>(&[seq]).0)
    }

    /// Count distinct hits and minimizers across both mates.
    pub fn score_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> FilterScore {
        Self::score(self.classify_seqs::<true, false>(&[seq1, seq2]).0)
    }

    fn diagnostics(&mut self, seqs: &[&[u8]]) -> FilterDiagnostics {
        let (decision, hit_kmers) = self.classify_seqs::<true, true>(seqs);
        FilterDiagnostics {
            score: Self::score(decision),
            hit_kmers,
        }
    }

    fn score(decision: FilterDecision) -> FilterScore {
        FilterScore {
            is_match: decision.is_match,
            keep: decision.keep,
            hit_count: decision.hit_count_lower_bound,
            minimizer_count: decision.minimizer_count_upper_bound,
        }
    }

    /// Pool distinct minimizer hits across mates
    fn classify_seqs<const EXACT: bool, const HIT_KMERS: bool>(
        &mut self,
        seqs: &[&[u8]],
    ) -> (FilterDecision, Vec<String>) {
        self.seen_minimizers.clear();
        self.minimizers.clear();
        let mut hit_kmers = Vec::new();

        // Pool mates first: a per-mate `stop_at` would understate the denominator
        for &seq in seqs {
            let seq = self.seq_for_filter(seq);
            self.minimizers.extend(seq);
        }

        // Positions bound the distinct count and required_hits never decreases, so no
        // denominator can ask for more hits than this: stop probing once it's reached
        let positional_minimizer_count = self.minimizers.values().len();
        let stop_at = if EXACT {
            usize::MAX
        } else {
            self.required_hits(positional_minimizer_count)
        };
        let hit_count_lower_bound = self.count_hits_until::<HIT_KMERS>(stop_at, &mut hit_kmers);

        // Avoid costly distinct counting where possible
        let minimizer_count_upper_bound = if EXACT
            || self.needs_exact_minimizer_count(hit_count_lower_bound, positional_minimizer_count)
        {
            self.count_distinct_minimizers()
        } else {
            positional_minimizer_count
        };

        let is_match = self.is_index_match(hit_count_lower_bound, minimizer_count_upper_bound);
        let keep = is_match != self.params.deplete;

        (
            FilterDecision {
                is_match,
                keep,
                hit_count_lower_bound,
                minimizer_count_upper_bound,
            },
            hit_kmers,
        )
    }

    /// Exact distinct minimizer count pooled across mates. Expensive but should be rare
    fn count_distinct_minimizers(&mut self) -> usize {
        self.seen_minimizers.clear();
        self.insert_buffer_minimizers();
        self.seen_minimizers.len()
    }

    fn insert_buffer_minimizers(&mut self) {
        match (self.minimizers.values(), &mut self.seen_minimizers) {
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
    fn count_hits_until<const HIT_KMERS: bool>(
        &mut self,
        stop_at: usize,
        hit_kmers: &mut Vec<String>,
    ) -> usize {
        let mut distinct_hits = 0;
        let kmer_length = self.minimizers.kmer_length();
        match (
            self.minimizers.values(),
            &self.index.minimizers,
            &mut self.seen_minimizers,
        ) {
            (MinimizerVec::U64(vec), IndexStorage::ExactU64(set), SeenMinimizers::U64(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) {
                        distinct_hits += 1;
                        if HIT_KMERS {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u64(minimizer, kmer_length))
                                    .to_string(),
                            );
                        }
                        if distinct_hits >= stop_at {
                            break;
                        }
                    }
                }
            }
            (MinimizerVec::U64(vec), IndexStorage::Fuse(filter), SeenMinimizers::U64(seen)) => {
                for &minimizer in vec {
                    if filter.contains(minimizer) && seen.insert(minimizer) {
                        distinct_hits += 1;
                        if HIT_KMERS {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u64(minimizer, kmer_length))
                                    .to_string(),
                            );
                        }
                        if distinct_hits >= stop_at {
                            break;
                        }
                    }
                }
            }
            (MinimizerVec::U128(vec), IndexStorage::ExactU128(set), SeenMinimizers::U128(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) {
                        distinct_hits += 1;
                        if HIT_KMERS {
                            hit_kmers.push(
                                String::from_utf8_lossy(&decode_u128(minimizer, kmer_length))
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

    #[test]
    fn short_read_behavior_matches_filter_modes() {
        let index = Arc::new(Index::from_minimizers(31, 15, MinimizerVec::U64(vec![])).unwrap());
        let params = FilterParams {
            deplete: true,
            abs_threshold: 2,
            rel_threshold: 0.01,
            prefix_length: 0,
        };
        let mut kernel = FilterKernel::new(Arc::clone(&index), params).unwrap();
        assert!(kernel.classify_read(b"ACGT").keep);

        let mut kernel = FilterKernel::new(
            Arc::clone(&index),
            FilterParams {
                deplete: false,
                ..params
            },
        )
        .unwrap();
        assert!(!kernel.classify_read(b"ACGT").keep);
    }

    #[test]
    fn zero_minimizers_never_match() {
        let index = Arc::new(Index::from_minimizers(31, 15, MinimizerVec::U64(vec![])).unwrap());
        for deplete in [false, true] {
            let mut kernel = FilterKernel::new(
                Arc::clone(&index),
                FilterParams {
                    deplete,
                    abs_threshold: 1,
                    rel_threshold: 0.0,
                    prefix_length: 0,
                },
            )
            .unwrap();
            assert_eq!(kernel.classify_read(b"").keep, deplete);
        }
    }

    #[test]
    fn relative_threshold_is_a_minimum_proportion() {
        let index = Arc::new(Index::from_minimizers(31, 15, MinimizerVec::U64(vec![])).unwrap());
        let kernel = FilterKernel::new(
            Arc::clone(&index),
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
        let index = Arc::new(Index::from_minimizers(31, 15, MinimizerVec::U64(vec![])).unwrap());
        let kernel = FilterKernel::new(
            Arc::clone(&index),
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
    fn index_of(reference: &[u8], kmer_length: u8, window_size: u8) -> Arc<Index> {
        let mut minimizers = Minimizers::new(kmer_length, window_size).unwrap();
        let values = minimizers.compute(reference).clone();
        Arc::new(Index::from_minimizers(kmer_length, window_size, values).unwrap())
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
                        Arc::clone(&index),
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
                        let fast = kernel.classify_read(read);
                        let exact = kernel.classify_read_with_diagnostics(read).score;
                        assert_eq!(fast.keep, exact.keep, "read {i}, {context}");
                        assert_hits_sufficient(
                            &kernel,
                            &fast,
                            &exact,
                            &format!("read {i}, {context}"),
                        );
                        substitution_observed |=
                            fast.minimizer_count_upper_bound > exact.minimizer_count;
                        early_exit_observed |= fast.hit_count_lower_bound < exact.hit_count;
                    }

                    for (i, first) in reads.iter().enumerate() {
                        for (j, second) in reads.iter().enumerate() {
                            let fast = kernel.classify_pair(first, second);
                            let exact = kernel.classify_pair_with_diagnostics(first, second).score;
                            assert_eq!(fast.keep, exact.keep, "pair {i}/{j}, {context}");
                            assert_hits_sufficient(
                                &kernel,
                                &fast,
                                &exact,
                                &format!("pair {i}/{j}, {context}"),
                            );
                            substitution_observed |=
                                fast.minimizer_count_upper_bound > exact.minimizer_count;
                            early_exit_observed |= fast.hit_count_lower_bound < exact.hit_count;
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
        exact: &FilterScore,
        context: &str,
    ) {
        assert!(fast.hit_count_lower_bound <= exact.hit_count, "{context}");
        assert!(
            fast.hit_count_lower_bound == exact.hit_count
                || fast.hit_count_lower_bound >= kernel.required_hits(exact.minimizer_count),
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
            Arc::clone(&index),
            FilterParams {
                deplete: true,
                abs_threshold: 1,
                rel_threshold: 0.0,
                prefix_length: 0,
            },
        )
        .unwrap();

        assert!(!kernel.classify_read(reference).keep);
        let decision = kernel.classify_read(b"ACG");
        assert_eq!(decision.hit_count_lower_bound, 0);
        assert_eq!(decision.minimizer_count_upper_bound, 0);
        assert!(decision.keep);

        assert!(!kernel.classify_pair(reference, reference).keep);
        let decision = kernel.classify_pair(b"ACG", b"");
        assert_eq!(decision.hit_count_lower_bound, 0);
        assert_eq!(decision.minimizer_count_upper_bound, 0);
        assert!(decision.keep);
    }

    /// `needs_exact_minimizer_count` relies on `rel_required_hits` never decreasing
    #[test]
    fn rel_required_hits_never_decreases() {
        let index = Arc::new(Index::from_minimizers(31, 15, MinimizerVec::U64(vec![])).unwrap());
        for &rel_threshold in &[0.0, 0.01, 0.1, 0.5, 0.99, 1.0] {
            let kernel = FilterKernel::new(
                Arc::clone(&index),
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
        let index = Arc::new(Index::from_minimizers(3, 1, MinimizerVec::U64(vec![0])).unwrap());
        let mut kernel = FilterKernel::new(
            Arc::clone(&index),
            FilterParams {
                deplete: false,
                abs_threshold: 1,
                rel_threshold: 1.0,
                prefix_length: 0,
            },
        )
        .unwrap();

        let decision = kernel.classify_pair(b"AAAAA", b"AAAAA");
        assert_eq!(decision.hit_count_lower_bound, 1);
        assert_eq!(decision.minimizer_count_upper_bound, 1);
        assert!(decision.keep);
    }
}
