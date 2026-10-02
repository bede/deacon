//! # Deacon
//!
//! A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format,
//! built for efficient host depletion (*deacon*-tamination).
//!
//! This crate provides both a library and a binary for filtering nucleotide sequences.
//!
#![doc = include_str!("../README.md")]

#[cfg(feature = "io")]
mod dedup_vec;
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
pub use filter::{FilterDecision, FilterDiagnostics, FilterKernel, FilterParams, FilterScore};
#[cfg(feature = "io")]
pub use filter_io::{
    DEFAULT_CBQ_BLOCK_SIZE_MIB, FilterConfig, FilterSummary, filter_files, load_filter_index,
};
pub use index::{Index, IndexKind};
pub use index_format::{load_index, load_index_from_path, write_index, write_index_to_path};
#[cfg(feature = "io")]
pub use index_ops::{
    BuildConfig, build as index_build, diff as index_diff, dump as index_dump,
    filter as index_filter, freeze as index_freeze, info as index_info,
    intersect as index_intersect, union as index_union,
};
pub use minimizers::{
    ComplexityAlgorithm, DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, MinimizerVec, Minimizers,
    decode_u64, decode_u128,
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
