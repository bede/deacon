//! # Deacon
//!
//! A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format,
//! built for efficient host depletion (*deacon*-tamination).
//!
//! This crate provides both a library and a binary for filtering nucleotide sequences.
//!
#![doc = include_str!("../README.md")]

// Re-export public functionality
#[cfg(feature = "cli")]
mod filter;
mod filter_kernel;
mod index;
mod minimizers;

// Public API
#[cfg(feature = "cli")]
pub use filter::{FilterRunConfig, FilterSummary, run as run_filter, run_with_index};
pub use filter_kernel::{FilterDecision, FilterKernel, FilterParams};
#[cfg(feature = "fetch")]
pub use index::fetch as index_fetch;
pub use index::{
    IndexHeader, dump_minimizers, load_index_auto, load_index_from_path_auto,
    load_minimizers_from_path,
};
#[cfg(feature = "cli")]
pub use index::{
    build as index_build, diff as index_diff, dump as index_dump, filter as index_filter,
    freeze as index_freeze, info as index_info, intersect as index_intersect, union as index_union,
};
pub use minimizers::{
    Buffers, DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, KmerHasher, compute_minimizers, decode_u64,
    decode_u128, fill_minimizers,
};

#[cfg(feature = "cli")]
use anyhow::Result;
use std::collections::HashSet;
use std::hash::BuildHasher;
#[cfg(feature = "cli")]
use std::path::{Path, PathBuf};

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
    ) {
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
            MinimizerSet::Fuse(_) => {
                panic!("Complexity filtering is not supported on BFF indexes; use an exact index")
            }
        }
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

#[cfg(feature = "cli")]
pub struct FilterConfig<'a> {
    /// Minimizer index file path
    pub minimizers_path: &'a Path,

    /// Path to input fastx file (or - for stdin)
    pub input_path: &'a str,

    /// Path to optional second paired fastx file (or - for interleaved stdin)
    pub input2_path: Option<&'a str>,

    /// Path to output fastx file (None for stdout; detects .gz and .zst)
    pub output_path: Option<&'a Path>,

    /// Path to optional second output fastx file for paired reads (detects .gz and .zst)
    pub output2_path: Option<&'a str>,

    /// Absolute threshold for filtering sequences
    pub abs_threshold: usize,

    /// Relative threshold for filtering sequences (0.0-1.0)
    pub rel_threshold: f64,

    /// Consider only the first N nucleotides per sequence (0 = entire sequence)
    pub prefix_length: usize,

    /// Discard index minimizers below this kdust complexity threshold (None = disabled)
    pub complexity_threshold: Option<f32>,

    /// Path to JSON summary file
    pub summary_path: Option<&'a PathBuf>,

    /// Deplete mode (remove sequences WITH matches, original deacon behavior)
    pub deplete: bool,

    /// Replace sequence headers with sequential numbers (1, 2, 3...)
    pub rename: bool,

    /// Replace headers with sequential numbers followed by random u64 (1-12345, 2-67890, ...)
    pub rename_random: bool,

    /// Force FASTA output (discards quality scores)
    pub output_fasta: bool,

    /// Number of execution threads (0 = auto)
    pub threads: u16,

    /// Compression level for output files (1-22 for zst, 1-9 for gz)
    pub compression_level: u8,

    /// Number of threads for compression (0 = auto-calculate as ceil(total/2))
    pub compression_threads: u16,

    /// Debug mode: output sequences with minimizer hits to stderr
    pub debug: bool,

    /// Suppress progress reporting
    pub quiet: bool,
}

#[cfg(feature = "cli")]
impl FilterConfig<'_> {
    /// Filter with this configuration
    pub fn execute(&self) -> Result<()> {
        filter::run(self)?;
        Ok(())
    }
}

#[cfg(feature = "cli")]
pub struct IndexConfig {
    /// Path to input fastx file
    pub input_path: PathBuf,

    /// K-mer length used for indexing
    pub kmer_length: u8,

    /// Minimizer window size used for indexing
    pub window_size: u8,

    /// Path to output file (None for stdout)
    pub output_path: Option<PathBuf>,

    /// Number of execution threads (0 = auto)
    pub threads: u16,

    /// Suppress per-sequence progress output
    pub quiet: bool,
}

#[cfg(feature = "cli")]
impl IndexConfig {
    /// Validate k-mer and window size constraints
    pub fn validate(&self) -> Result<()> {
        let k = self.kmer_length as usize;
        let w = self.window_size as usize;

        // Check constraints: k <= 61, k+w <= 96, k+w even (ensures k odd and k+w-1 odd)
        if k > 61 || k + w > 96 || !(k + w).is_multiple_of(2) {
            return Err(anyhow::anyhow!(
                "Invalid k-w combination: k={}, w={}, k+w={} (constraints: k<=61, k+w<=96, k+w even)",
                k,
                w,
                k + w
            ));
        }

        Ok(())
    }

    /// Execute index build with this configuration
    pub fn execute(&self) -> Result<()> {
        index::build(self)
    }
}
