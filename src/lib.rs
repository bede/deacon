//! # Deacon
//!
//! A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format,
//! built for efficient host depletion (*deacon*-tamination).
//!
//! This crate provides both a library and a binary for filtering nucleotide sequences.
//!
#![doc = include_str!("../README.md")]

// Re-export public functionality
pub mod filter;
pub mod index;
pub mod index_format;
pub mod minimizers;

// Keep these modules as internal implementation details
mod index_build;
mod index_diff;
mod index_info;
mod index_union;

// Re-export the important structures and functions for library users
pub use filter::{FilterLog, run as run_filter};
pub use index::{
    build as build_index, diff as diff_index, info as index_info, union as union_index,
};
pub use index_format::IndexHeader;
pub use minimizers::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, compute_minimizer_hashes, fill_minimizer_hashes,
};

use anyhow::Result;
use rustc_hash::FxHashSet;
use std::path::{Path, PathBuf};

/// Configuration for filtering operations
pub struct FilterConfig {
    /// Minimizer index file path
    pub minimizers_path: PathBuf,

    /// Path to input fastx file (or - for stdin)
    pub input_path: String,

    /// Path to optional second paired fastx file (or - for interleaved stdin)
    pub input2_path: Option<String>,

    /// Path to output fastx file (or - for stdout; detects .gz and .zst)
    pub output_path: String,

    /// Minimum number of minimizer matches per query sequence (pair)
    pub min_matches: usize,

    /// Consider only the first N nucleotides per sequence (0 = entire sequence)
    pub prefix_length: usize,

    /// Path to JSON output log file
    pub log_path: Option<PathBuf>,

    /// Invert filtering (keep sequences WITH matches rather than those WITHOUT)
    pub invert: bool,

    /// Replace sequence headers with sequential numbers (1, 2, 3...)
    pub rename: bool,
}

impl FilterConfig {
    /// Create a new filtering configuration with the specified minimizers path
    pub fn new<P: AsRef<Path>>(minimizers_path: P) -> Self {
        Self {
            minimizers_path: minimizers_path.as_ref().to_path_buf(),
            input_path: "-".to_string(),
            input2_path: None,
            output_path: "-".to_string(),
            min_matches: 2,
            prefix_length: 0,
            log_path: None,
            invert: false,
            rename: false,
        }
    }

    /// Set the input path
    pub fn with_input<S: Into<String>>(mut self, input_path: S) -> Self {
        self.input_path = input_path.into();
        self
    }

    /// Set the second input path for paired-end reads
    pub fn with_input2<S: Into<String>>(mut self, input2_path: S) -> Self {
        self.input2_path = Some(input2_path.into());
        self
    }

    /// Set the output path
    pub fn with_output<S: Into<String>>(mut self, output_path: S) -> Self {
        self.output_path = output_path.into();
        self
    }

    /// Set the minimum matches threshold
    pub fn with_min_matches(mut self, min_matches: usize) -> Self {
        self.min_matches = min_matches;
        self
    }

    /// Set the prefix length
    pub fn with_prefix_length(mut self, prefix_length: usize) -> Self {
        self.prefix_length = prefix_length;
        self
    }

    /// Set the log file path
    pub fn with_log<P: AsRef<Path>>(mut self, log_path: P) -> Self {
        self.log_path = Some(log_path.as_ref().to_path_buf());
        self
    }

    /// Set the invert flag
    pub fn with_invert(mut self, invert: bool) -> Self {
        self.invert = invert;
        self
    }

    /// Set the rename flag
    pub fn with_rename(mut self, rename: bool) -> Self {
        self.rename = rename;
        self
    }

    /// Execute the filtering operation with this configuration
    pub fn execute(&self) -> Result<()> {
        filter::run(
            &self.minimizers_path,
            &self.input_path,
            self.input2_path.as_deref(),
            &self.output_path,
            self.min_matches,
            self.prefix_length,
            self.log_path.as_ref(),
            self.invert,
            self.rename,
        )
    }
}

/// Configuration for index building operations
pub struct IndexConfig {
    /// Path to input fastx file
    pub input_path: PathBuf,

    /// K-mer length used for indexing
    pub kmer_length: usize,

    /// Minimizer window length used for indexing
    pub window_length: usize,

    /// Path to output file (None for stdout)
    pub output_path: Option<PathBuf>,
}

impl IndexConfig {
    /// Create a new index configuration with the specified input path
    pub fn new<P: AsRef<Path>>(input_path: P) -> Self {
        Self {
            input_path: input_path.as_ref().to_path_buf(),
            kmer_length: DEFAULT_KMER_LENGTH,
            window_length: DEFAULT_WINDOW_SIZE,
            output_path: None,
        }
    }

    /// Set the k-mer length
    pub fn with_kmer_length(mut self, kmer_length: usize) -> Self {
        self.kmer_length = kmer_length;
        self
    }

    /// Set the window length
    pub fn with_window_length(mut self, window_length: usize) -> Self {
        self.window_length = window_length;
        self
    }

    /// Set the output path
    pub fn with_output<P: AsRef<Path>>(mut self, output_path: P) -> Self {
        self.output_path = Some(output_path.as_ref().to_path_buf());
        self
    }

    /// Execute the index building operation with this configuration
    pub fn execute(&self) -> Result<()> {
        index_build::build(
            &self.input_path,
            self.kmer_length,
            self.window_length,
            self.output_path.clone(),
        )
    }
}

/// Load minimizer hashes from a file
pub fn load_minimizers<P: AsRef<Path>>(path: P) -> Result<(FxHashSet<u64>, IndexHeader)> {
    index::load_minimizer_hashes(&path)
}

/// Write minimizer hashes to a file
pub fn write_minimizers(
    minimizers: &FxHashSet<u64>,
    header: &IndexHeader,
    output_path: Option<&PathBuf>,
) -> Result<()> {
    index::write_minimizers(minimizers, header, output_path)
}
