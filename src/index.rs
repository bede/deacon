use crate::index_format::IndexHeader;
use crate::minimizers::{decode_u64, decode_u128};
use crate::{ComplexityAlgorithm, MinimizerVec, validate_unit_interval};
use hashbrown::HashSet;
use std::hash::BuildHasher;

/// Index storage and membership guarantees.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IndexKind {
    Exact,
    /// No false negatives. False-positive rate ~2^-fingerprint_bits per query.
    Fuse { fingerprint_bits: u8 },
}

/// Minimizer index with validated k/w. Share across workers with `Arc<Index>`.
pub struct Index {
    pub(crate) minimizers: IndexStorage,
    pub(crate) header: IndexHeader,
}

impl std::fmt::Debug for Index {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Index")
            .field("k", &self.kmer_length())
            .field("w", &self.window_size())
            .field("kind", &self.kind())
            .field("len", &self.len())
            .finish()
    }
}

impl Index {
    /// Build an exact index from canonical packed k-mers, deduplicating values.
    /// Use their original k/w and u64 for k <= 32, otherwise u128.
    pub fn from_minimizers(
        kmer_length: u8,
        window_size: u8,
        minimizers: MinimizerVec,
    ) -> anyhow::Result<Self> {
        let header = IndexHeader::new(kmer_length, window_size);
        header.validate()?;
        let bits = 2 * u32::from(kmer_length);
        let minimizers = match minimizers {
            MinimizerVec::U64(values) if kmer_length <= 32 => {
                anyhow::ensure!(
                    bits == 64 || values.iter().all(|&v| v >> bits == 0),
                    "Minimizer exceeds k-mer length"
                );
                IndexStorage::ExactU64(values.into_iter().collect())
            }
            MinimizerVec::U128(values) if kmer_length > 32 => {
                anyhow::ensure!(
                    values.iter().all(|&v| v >> bits == 0),
                    "Minimizer exceeds k-mer length"
                );
                IndexStorage::ExactU128(values.into_iter().collect())
            }
            _ => anyhow::bail!("Minimizer width does not match k-mer length"),
        };
        Ok(Self { minimizers, header })
    }

    pub fn kmer_length(&self) -> u8 {
        self.header.kmer_length
    }
    pub fn window_size(&self) -> u8 {
        self.header.window_size
    }
    /// Distinct minimizers (inserted keys for fuse storage).
    pub fn len(&self) -> usize {
        self.minimizers.len()
    }
    pub fn is_empty(&self) -> bool {
        self.minimizers.is_empty()
    }
    pub fn kind(&self) -> IndexKind {
        match &self.minimizers {
            IndexStorage::Fuse(filter) => IndexKind::Fuse {
                fingerprint_bits: filter.fingerprint_bits(),
            },
            _ => IndexKind::Exact,
        }
    }

    /// Drop minimizers below the complexity threshold, or retain them if inverted.
    /// Requires an exact index.
    pub fn retain_complexity(
        &mut self,
        algorithm: ComplexityAlgorithm,
        threshold: f32,
        invert: bool,
    ) -> anyhow::Result<()> {
        self.minimizers
            .retain_complexity(self.kmer_length(), algorithm, threshold, invert)
    }
}

/// BuildHasher using rapidhash with fixed seed for fast init
#[derive(Clone, Default)]
pub(crate) struct FixedRapidHasher;

impl BuildHasher for FixedRapidHasher {
    type Hasher = rapidhash::fast::RapidHasher<'static>;

    fn build_hasher(&self) -> Self::Hasher {
        rapidhash::fast::SeedableState::fixed().build_hasher()
    }
}

/// RapidHashSet using rapidhash with fixed seed for fast init
///
/// We directly use the hashbrown version for stability of the number of slots in the data structure.
pub(crate) type RapidHashSet<T> = HashSet<T, FixedRapidHasher>;

/// Binary fuse filter (BFF) index: no false negatives, k<=32.
/// 16-bit fingerprints use ~18 bits/key (FP rate ~2^-16), 32-bit ~36 bits/key (FP rate ~2^-32).
pub(crate) enum FuseFilterKind {
    Bits16(xorf::BinaryFuse16),
    Bits32(xorf::BinaryFuse32),
}

pub(crate) struct FuseFilter {
    pub(crate) filter: FuseFilterKind,
    /// Distinct minimizers inserted (not the filter's fingerprint count)
    pub(crate) key_count: usize,
}

impl FuseFilter {
    /// Fingerprint width in bits (16 or 32)
    pub fn fingerprint_bits(&self) -> u8 {
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

/// Exact sets or binary fuse filters.
pub(crate) enum IndexStorage {
    ExactU64(RapidHashSet<u64>),
    ExactU128(RapidHashSet<u128>),
    Fuse(FuseFilter),
}

impl IndexStorage {
    pub fn len(&self) -> usize {
        match self {
            IndexStorage::ExactU64(set) => set.len(),
            IndexStorage::ExactU128(set) => set.len(),
            IndexStorage::Fuse(f) => f.key_count,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[cfg(feature = "io")]
    pub fn extend(&mut self, other: Self) {
        match (self, other) {
            (IndexStorage::ExactU64(self_set), IndexStorage::ExactU64(other_set)) => {
                self_set.extend(other_set);
            }
            (IndexStorage::ExactU128(self_set), IndexStorage::ExactU128(other_set)) => {
                self_set.extend(other_set);
            }
            (IndexStorage::Fuse(_), _) | (_, IndexStorage::Fuse(_)) => {
                panic!("Set algebra is not supported on BFF indexes; use an exact index")
            }
            _ => panic!("Cannot extend U64 set with U128 set or vice versa"),
        }
    }

    #[cfg(feature = "io")]
    pub fn remove_all(&mut self, other: &Self) {
        match (self, other) {
            (IndexStorage::ExactU64(self_set), IndexStorage::ExactU64(other_set)) => {
                for val in other_set {
                    self_set.remove(val);
                }
            }
            (IndexStorage::ExactU128(self_set), IndexStorage::ExactU128(other_set)) => {
                for val in other_set {
                    self_set.remove(val);
                }
            }
            (IndexStorage::Fuse(_), _) | (_, IndexStorage::Fuse(_)) => {
                panic!("Set algebra is not supported on BFF indexes; use an exact index")
            }
            _ => panic!("Cannot remove U128 minimizers from U64 set or vice versa"),
        }
    }

    #[cfg(feature = "io")]
    pub fn intersect(&mut self, other: &Self) {
        match (self, other) {
            (IndexStorage::ExactU64(self_set), IndexStorage::ExactU64(other_set)) => {
                self_set.retain(|val| other_set.contains(val));
            }
            (IndexStorage::ExactU128(self_set), IndexStorage::ExactU128(other_set)) => {
                self_set.retain(|val| other_set.contains(val));
            }
            (IndexStorage::Fuse(_), _) | (_, IndexStorage::Fuse(_)) => {
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
            IndexStorage::ExactU64(set) => set.retain(|&v| {
                let c = match algorithm {
                    ComplexityAlgorithm::Kdust => calculate_kdust(v as u128, kmer_length),
                    ComplexityAlgorithm::Shannon => {
                        calculate_scaled_entropy(&decode_u64(v, kmer_length), kmer_length)
                    }
                };
                keep(c)
            }),
            IndexStorage::ExactU128(set) => set.retain(|&v| {
                let c = match algorithm {
                    ComplexityAlgorithm::Kdust => calculate_kdust(v, kmer_length),
                    ComplexityAlgorithm::Shannon => {
                        calculate_scaled_entropy(&decode_u128(v, kmer_length), kmer_length)
                    }
                };
                keep(c)
            }),
            IndexStorage::Fuse(_) => anyhow::bail!(
                "Complexity filtering is not supported on BFF indexes; use an exact index"
            ),
        }
        Ok(())
    }
}
