//! # Deacon
//!
//! A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format,
//! built for efficient host depletion (*deacon*-tamination).
//!
//! This crate provides both a library and a binary for filtering nucleotide sequences.
//!
#![doc = include_str!("../README.md")]

// Re-export public functionality
#[cfg(feature = "io")]
mod dedupping_vec;
#[cfg(feature = "fetch")]
mod fetch;
mod filter;
#[cfg(feature = "io")]
mod filter_io;
mod index;
mod index_format;
#[cfg(feature = "io")]
mod index_ops;
mod minimizers;

// Public API
#[cfg(feature = "fetch")]
pub use fetch::fetch as index_fetch;
pub use filter::{FilterDecision, FilterKernel, FilterParams};
#[cfg(feature = "io")]
pub use filter_io::{
    DEFAULT_CBQ_BLOCK_SIZE_MIB, FilterRunConfig, FilterSummary, run as run_filter, run_with_index,
};
pub use index_format::{
    IndexHeader, dump_minimizers, load_index_auto, load_index_from_path_auto,
    load_minimizers_from_path,
};
#[cfg(feature = "io")]
pub use index_ops::{
    build as index_build, current_index_path, diff as index_diff, dump as index_dump,
    filter as index_filter, freeze as index_freeze, info as index_info,
    intersect as index_intersect, union as index_union,
};
pub use minimizers::{
    Buffers, DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, KmerHasher, compute_minimizers, decode_u64,
    decode_u128, fill_minimizers,
};

/// Validate that a threshold is within 0.0..=1.0
pub fn validate_unit_interval<T: Into<f64> + Copy + std::fmt::Display>(
    name: &str,
    value: T,
) -> anyhow::Result<()> {
    if !(0.0..=1.0).contains(&value.into()) {
        anyhow::bail!("{name} must be between 0.0 and 1.0 inclusive; got {value}");
    }
    Ok(())
}

#[cfg(feature = "io")]
use anyhow::Result;
#[cfg(feature = "io")]
use std::path::{Path, PathBuf};

#[cfg(feature = "io")]
pub struct FilterConfig<'a> {
    /// Minimizer index file path
    pub minimizers_path: &'a Path,

    /// Path to input fastx file (or - for stdin)
    pub input_path: &'a str,

    /// Path to optional second paired fastx file
    pub input2_path: Option<&'a str>,

    /// Treat input_path as an interleaved paired stream
    pub interleaved: bool,

    /// Validate paired record names (Illumina CASAVA or /1 /2 suffixes)
    pub check_pairs: bool,

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

    /// Deplete mode (remove sequences WITH matches, original deacon behaviour)
    pub deplete: bool,

    /// Replace sequence headers with incrementing numbers (1, 2, 3...)
    pub rename: bool,

    /// Emit fasta or quality-free cbq regardless of input format (implied by a .cba output)
    pub discard_quality: bool,

    /// Preserve input record ordering (deterministic, slightly slower)
    pub ordered: bool,

    /// Number of execution threads (0 = auto)
    pub threads: u16,

    /// Compression level for output files (1-22 for zst, 1-9 for gz)
    pub compression_level: u8,

    /// cbq output block size in MiB (raised to the cbq input's block size if larger)
    pub cbq_block_size: u16,

    /// Number of threads for compression (0 = auto-calculate as ceil(total/2))
    pub compression_threads: u16,

    /// Per-record hit TSV path
    pub debug: Option<&'a PathBuf>,

    /// Suppress progress reporting
    pub quiet: bool,
}

#[cfg(feature = "io")]
impl FilterConfig<'_> {
    /// Filter with this configuration
    pub fn execute(&self) -> Result<()> {
        filter_io::run(self)?;
        Ok(())
    }
}

#[cfg(feature = "io")]
#[derive(Clone)]
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

#[cfg(feature = "io")]
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
        index_ops::build(self)
    }
}

pub(crate) use index::MinimizerVecVec;
pub use index::{
    ComplexityAlgorithm, FixedRapidHasher, FuseFilter, FuseFilterKind, MinimizerSet, MinimizerVec,
    RapidHashSet,
};
