use crate::{decode_u64, decode_u128, validate_unit_interval};
use hashbrown::HashSet;
use std::hash::BuildHasher;

/// BuildHasher using rapidhash with fixed seed for fast init
#[derive(Clone, Default)]
pub struct FixedRapidHasher;

impl BuildHasher for FixedRapidHasher {
    type Hasher = rapidhash::fast::RapidHasher<'static>;

    fn build_hasher(&self) -> Self::Hasher {
        rapidhash::fast::SeedableState::fixed().build_hasher()
    }
}

/// RapidHashSet using rapidhash with fixed seed for fast init
///
/// We directly use the hashbrown version for stability of the number of slots in the data structure.
pub type RapidHashSet<T> = HashSet<T, FixedRapidHasher>;

/// Binary fuse filter (BFF) index: no false negatives, k<=32.
/// 16-bit fingerprints use ~18 bits/key (FP rate ~2^-16); 32-bit use ~36 bits/key (FP rate ~2^-32).
pub enum FuseFilterKind {
    Bits16(xorf::BinaryFuse16),
    Bits32(xorf::BinaryFuse32),
}

pub struct FuseFilter {
    pub filter: FuseFilterKind,
    /// Distinct minimizers inserted (not the filter's fingerprint count)
    pub key_count: usize,
}

impl FuseFilter {
    /// Fingerprint width in bits (16 or 32)
    pub fn filter_bits(&self) -> u8 {
        match &self.filter {
            FuseFilterKind::Bits16(_) => 16,
            FuseFilterKind::Bits32(_) => 32,
        }
    }

    /// Test membership (may report false positives)
    #[inline]
    pub fn contains(&self, minimizer: u64) -> bool {
        use xorf::Filter;
        match &self.filter {
            FuseFilterKind::Bits16(f) => f.contains(&minimizer),
            FuseFilterKind::Bits32(f) => f.contains(&minimizer),
        }
    }
}

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets and the BFF filter
pub enum MinimizerSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
    Fuse(FuseFilter),
}

/// Complexity measure used by `index filter`
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[serde(rename_all = "lowercase")]
pub enum ComplexityAlgorithm {
    #[default]
    Kdust,
    Shannon,
}

impl std::fmt::Display for ComplexityAlgorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            ComplexityAlgorithm::Kdust => "kdust",
            ComplexityAlgorithm::Shannon => "shannon",
        })
    }
}

impl MinimizerSet {
    pub fn len(&self) -> usize {
        match self {
            MinimizerSet::U64(set) => set.len(),
            MinimizerSet::U128(set) => set.len(),
            MinimizerSet::Fuse(f) => f.key_count,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn is_u64(&self) -> bool {
        // BFF is k<=32, hence u64
        matches!(self, MinimizerSet::U64(_) | MinimizerSet::Fuse(_))
    }

    /// Test u64 membership (may report false positives for BFF)
    #[inline]
    pub fn contains_u64(&self, minimizer: u64) -> bool {
        match self {
            MinimizerSet::U64(set) => set.contains(&minimizer),
            MinimizerSet::Fuse(f) => f.contains(minimizer),
            MinimizerSet::U128(_) => unreachable!("u64 minimizer queried against a u128 set"),
        }
    }

    /// Test u128 membership
    #[inline]
    pub fn contains_u128(&self, minimizer: u128) -> bool {
        match self {
            MinimizerSet::U128(set) => set.contains(&minimizer),
            MinimizerSet::U64(_) => unreachable!("u128 minimizer queried against a u64 set"),
            MinimizerSet::Fuse(_) => {
                unreachable!("u128 minimizer queried against a BFF (k <= 32)")
            }
        }
    }

    /// Extend with another MinimizerSet (union operation)
    pub fn extend(&mut self, other: Self) {
        match (self, other) {
            (MinimizerSet::U64(self_set), MinimizerSet::U64(other_set)) => {
                self_set.extend(other_set);
            }
            (MinimizerSet::U128(self_set), MinimizerSet::U128(other_set)) => {
                self_set.extend(other_set);
            }
            (MinimizerSet::Fuse(_), _) | (_, MinimizerSet::Fuse(_)) => {
                panic!("Set algebra is not supported on BFF indexes; use an exact index")
            }
            _ => panic!("Cannot extend U64 set with U128 set or vice versa"),
        }
    }

    /// Remove minimizers from another set (diff operation)
    pub fn remove_all(&mut self, other: &Self) {
        match (self, other) {
            (MinimizerSet::U64(self_set), MinimizerSet::U64(other_set)) => {
                for val in other_set {
                    self_set.remove(val);
                }
            }
            (MinimizerSet::U128(self_set), MinimizerSet::U128(other_set)) => {
                for val in other_set {
                    self_set.remove(val);
                }
            }
            (MinimizerSet::Fuse(_), _) | (_, MinimizerSet::Fuse(_)) => {
                panic!("Set algebra is not supported on BFF indexes; use an exact index")
            }
            _ => panic!("Cannot remove U128 minimizers from U64 set or vice versa"),
        }
    }

    /// Keep only minimizers present in another set (intersection operation)
    pub fn intersect(&mut self, other: &Self) {
        match (self, other) {
            (MinimizerSet::U64(self_set), MinimizerSet::U64(other_set)) => {
                self_set.retain(|val| other_set.contains(val));
            }
            (MinimizerSet::U128(self_set), MinimizerSet::U128(other_set)) => {
                self_set.retain(|val| other_set.contains(val));
            }
            (MinimizerSet::Fuse(_), _) | (_, MinimizerSet::Fuse(_)) => {
                panic!("Set algebra is not supported on BFF indexes; use an exact index")
            }
            _ => panic!("Cannot intersect U64 set with U128 set or vice versa"),
        }
    }

    /// Keep minimizers with complexity >= threshold (or < if inverted)
    pub fn retain_complexity(
        &mut self,
        kmer_length: u8,
        algorithm: ComplexityAlgorithm,
        threshold: f32,
        invert: bool,
    ) -> anyhow::Result<()> {
        validate_unit_interval("complexity threshold", threshold)?;
        use crate::minimizers::{calculate_kdust, calculate_scaled_entropy};
        let keep = |c: f32| {
            if invert {
                c < threshold
            } else {
                c >= threshold
            }
        };
        match self {
            MinimizerSet::U64(set) => set.retain(|&v| {
                let c = match algorithm {
                    ComplexityAlgorithm::Kdust => calculate_kdust(v as u128, kmer_length),
                    ComplexityAlgorithm::Shannon => {
                        calculate_scaled_entropy(&decode_u64(v, kmer_length), kmer_length)
                    }
                };
                keep(c)
            }),
            MinimizerSet::U128(set) => set.retain(|&v| {
                let c = match algorithm {
                    ComplexityAlgorithm::Kdust => calculate_kdust(v, kmer_length),
                    ComplexityAlgorithm::Shannon => {
                        calculate_scaled_entropy(&decode_u128(v, kmer_length), kmer_length)
                    }
                };
                keep(c)
            }),
            MinimizerSet::Fuse(_) => anyhow::bail!(
                "Complexity filtering is not supported on BFF indexes; use an exact index"
            ),
        }
        Ok(())
    }
}

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
#[derive(Clone)]
pub enum MinimizerVec {
    U64(Vec<u64>),
    U128(Vec<u128>),
}

impl MinimizerVec {
    pub fn clear(&mut self) {
        match self {
            MinimizerVec::U64(v) => v.clear(),
            MinimizerVec::U128(v) => v.clear(),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            MinimizerVec::U64(v) => v.len(),
            MinimizerVec::U128(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            MinimizerVec::U64(v) => v.is_empty(),
            MinimizerVec::U128(v) => v.is_empty(),
        }
    }
}

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
#[derive(Clone)]
pub(crate) enum MinimizerVecVec {
    U64(Vec<Vec<u64>>),
    U128(Vec<Vec<u128>>),
}

impl MinimizerVecVec {
    pub fn len(&self) -> usize {
        match self {
            MinimizerVecVec::U64(v) => v.iter().map(|inner| inner.len()).sum(),
            MinimizerVecVec::U128(v) => v.iter().map(|inner| inner.len()).sum(),
        }
    }
}
