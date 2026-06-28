use crate::minimizers::{Buffers, KmerHasher, decode_u64, decode_u128, fill_minimizers_unchecked};
use crate::{MinimizerSet, MinimizerVec, RapidHashSet};

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
    pub total_minimizers: usize,
    pub hit_kmers: Vec<String>,
}

#[derive(Clone)]
enum SeenHits {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl SeenHits {
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
    seen_hits: SeenHits,
}

impl FilterKernel {
    pub fn new(kmer_length: u8, window_size: u8, params: FilterParams) -> Self {
        Self {
            params,
            kmer_length,
            window_size,
            hasher: KmerHasher::new(kmer_length as usize),
            buffers: if kmer_length <= 32 {
                Buffers::new_u64()
            } else {
                Buffers::new_u128()
            },
            seen_hits: SeenHits::new(kmer_length),
        }
    }

    #[inline]
    pub fn params(&self) -> FilterParams {
        self.params
    }

    #[inline]
    fn required_hits(&self, total_minimizers: usize) -> usize {
        let rel_required = if total_minimizers == 0 {
            0
        } else {
            ((self.params.rel_threshold * total_minimizers as f64).round() as usize).max(1)
        };
        self.params.abs_threshold.max(rel_required)
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
        self.seen_hits.clear();
        let mut total_minimizers = 0;
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
            total_minimizers += self.buffers.minimizers.len();
            self.count_buffer_hits_into_seen(index, debug, &mut hit_kmers);
        }

        let hit_count = self.seen_hits.len();
        // No minis (too short / all-ambiguous) means no possible index match: keep only when depleting
        let keep = if total_minimizers == 0 {
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

    fn count_buffer_hits_into_seen(
        &mut self,
        index: &MinimizerSet,
        debug: bool,
        hit_kmers: &mut Vec<String>,
    ) {
        match (&self.buffers.minimizers, index, &mut self.seen_hits) {
            (MinimizerVec::U64(vec), MinimizerSet::U64(set), SeenHits::U64(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) && debug {
                        hit_kmers.push(
                            String::from_utf8_lossy(&decode_u64(minimizer, self.kmer_length))
                                .to_string(),
                        );
                    }
                }
            }
            (MinimizerVec::U64(vec), MinimizerSet::Fuse(filter), SeenHits::U64(seen)) => {
                for &minimizer in vec {
                    if filter.contains(minimizer) && seen.insert(minimizer) && debug {
                        hit_kmers.push(
                            String::from_utf8_lossy(&decode_u64(minimizer, self.kmer_length))
                                .to_string(),
                        );
                    }
                }
            }
            (MinimizerVec::U128(vec), MinimizerSet::U128(set), SeenHits::U128(seen)) => {
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen.insert(minimizer) && debug {
                        hit_kmers.push(
                            String::from_utf8_lossy(&decode_u128(minimizer, self.kmer_length))
                                .to_string(),
                        );
                    }
                }
            }
            _ => unreachable!("minimizer width does not match index variant"),
        }
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
        let mut kernel = FilterKernel::new(31, 15, params);
        assert!(kernel.classify_read(&index, b"ACGT", false).keep);

        let mut kernel = FilterKernel::new(
            31,
            15,
            FilterParams {
                deplete: false,
                ..params
            },
        );
        assert!(!kernel.classify_read(&index, b"ACGT", false).keep);
    }
}
