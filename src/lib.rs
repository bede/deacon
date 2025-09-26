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
pub mod minimizers;
mod hashset;

// Re-export the important structures and functions for library users
pub use filter::{FilterSummary, run as run_filter};
pub use index::{
    IndexHeader, build as build_index, diff as diff_index, info as index_info, union as union_index,
};
pub use minimizers::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, compute_minimizer_hashes, fill_minimizer_hashes,
};

use anyhow::Result;
type FxHashSet = rustc_hash::FxHashSet<u64>;
use hashset::U64HashSet as HashSet;
use std::path::{Path, PathBuf};

pub struct FilterConfig<'a> {
    /// Minimizer index file path
    pub minimizers_path: &'a Path,

    /// Path to input fastx file (or - for stdin)
    pub input_path: &'a str,

    /// Path to optional second paired fastx file (or - for interleaved stdin)
    pub input2_path: Option<&'a str>,

    /// Path to output fastx file (or - for stdout; detects .gz and .zst)
    pub output_path: &'a str,

    /// Path to optional second output fastx file for paired reads (detects .gz and .zst)
    pub output2_path: Option<&'a str>,

    /// Absolute threshold for filtering sequences
    pub abs_threshold: usize,

    /// Relative threshold for filtering sequences (0.0-1.0)
    pub rel_threshold: f64,

    /// Consider only the first N nucleotides per sequence (0 = entire sequence)
    pub prefix_length: usize,

    /// Path to JSON summary file
    pub summary_path: Option<&'a PathBuf>,

    /// Deplete mode (remove sequences WITH matches, original deacon behavior)
    pub deplete: bool,

    /// Replace sequence headers with sequential numbers (1, 2, 3...)
    pub rename: bool,

    /// Number of execution threads (0 = auto)
    pub threads: usize,

    /// Compression level for output files (1-22 for zst, 1-9 for gz)
    pub compression_level: u8,

    /// Debug mode: output sequences with minimizer hits to stderr
    pub debug: bool,

    /// Suppress progress reporting
    pub quiet: bool,
}

impl<'a> FilterConfig<'a> {
    pub fn new(minimizers_path: &'a Path) -> Self {
        Self {
            minimizers_path,
            input_path: "-",
            input2_path: None,
            output_path: "-",
            output2_path: None,
            abs_threshold: 2,
            rel_threshold: 0.01,
            prefix_length: 0,
            summary_path: None,
            deplete: false,
            rename: false,
            threads: 0,           // Use all available threads by default
            compression_level: 2, // Default compression level
            debug: false,
            quiet: false,
        }
    }

    pub fn with_input(mut self, input_path: &'a str) -> Self {
        self.input_path = input_path;
        self
    }

    pub fn with_input2(mut self, input2_path: &'a str) -> Self {
        self.input2_path = Some(input2_path);
        self
    }

    pub fn with_output(mut self, output_path: &'a str) -> Self {
        self.output_path = output_path;
        self
    }

    pub fn with_output2(mut self, output2_path: &'a str) -> Self {
        self.output2_path = Some(output2_path);
        self
    }

    pub fn with_abs_threshold(mut self, abs_threshold: usize) -> Self {
        self.abs_threshold = abs_threshold;
        self
    }

    pub fn with_rel_threshold(mut self, rel_threshold: f64) -> Self {
        self.rel_threshold = rel_threshold;
        self
    }

    pub fn with_prefix_length(mut self, prefix_length: usize) -> Self {
        self.prefix_length = prefix_length;
        self
    }

    pub fn with_summary(mut self, summary_path: &'a PathBuf) -> Self {
        self.summary_path = Some(summary_path);
        self
    }

    pub fn with_deplete(mut self, deplete: bool) -> Self {
        self.deplete = deplete;
        self
    }

    pub fn with_rename(mut self, rename: bool) -> Self {
        self.rename = rename;
        self
    }

    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    pub fn with_compression_level(mut self, compression_level: u8) -> Self {
        self.compression_level = compression_level;
        self
    }

    pub fn with_debug(mut self, debug: bool) -> Self {
        self.debug = debug;
        self
    }

    pub fn with_quiet(mut self, quiet: bool) -> Self {
        self.quiet = quiet;
        self
    }

    /// Filter with this configuration
    pub fn execute(&self) -> Result<()> {
        filter::run(self)
    }
}

pub struct IndexConfig {
    /// Path to input fastx file
    pub input_path: PathBuf,

    /// K-mer length used for indexing
    pub kmer_length: u8,

    /// Minimizer window size used for indexing
    pub window_size: u8,

    /// Path to output file (None for stdout)
    pub output_path: Option<PathBuf>,

    /// Hash table pre-allocation capacity in millions
    pub capacity_millions: usize,

    /// Number of execution threads (0 = auto)
    pub threads: usize,

    /// Suppress per-sequence progress output
    pub quiet: bool,

    /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
    pub entropy_threshold: f32,
}

impl IndexConfig {
    /// Create a new index configuration with the specified input path
    pub fn new<P: AsRef<Path>>(input_path: P) -> Self {
        Self {
            input_path: input_path.as_ref().to_path_buf(),
            kmer_length: DEFAULT_KMER_LENGTH,
            window_size: DEFAULT_WINDOW_SIZE,
            output_path: None,
            capacity_millions: 400,
            threads: 8,
            quiet: false,
            entropy_threshold: 0.0,
        }
    }

    /// Set k-mer length
    pub fn with_kmer_length(mut self, kmer_length: u8) -> Self {
        self.kmer_length = kmer_length;
        self
    }

    /// Set window size
    pub fn with_window_size(mut self, window_size: u8) -> Self {
        self.window_size = window_size;
        self
    }

    /// Set output path
    pub fn with_output<P: AsRef<Path>>(mut self, output_path: P) -> Self {
        self.output_path = Some(output_path.as_ref().to_path_buf());
        self
    }

    /// Set hash table capacity in millions
    pub fn with_capacity_millions(mut self, capacity_millions: usize) -> Self {
        self.capacity_millions = capacity_millions;
        self
    }

    /// Set threads
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set quiet mode
    pub fn with_quiet(mut self, quiet: bool) -> Self {
        self.quiet = quiet;
        self
    }

    /// Set threshold for scaled entropy filtering at indexing time
    pub fn with_entropy_threshold(mut self, threshold: f32) -> Self {
        self.entropy_threshold = threshold;
        self
    }

    /// Execute index build with this configuration
    pub fn execute(&self) -> Result<()> {
        index::build(self)
    }
}

pub fn load_minimizers<P: AsRef<Path>>(path: P) -> Result<(HashSet, index::IndexHeader)> {
    index::load_minimizer_hashes(&path)
}

pub fn write_minimizers(
    minimizers: &HashSet,
    header: &index::IndexHeader,
    output_path: Option<&PathBuf>,
) -> Result<()> {
    index::write_minimizers(minimizers.len(), minimizers.iter(), header, output_path)
}
