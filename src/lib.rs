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

// Re-export the important structures and functions for library users
pub use filter::{FilterSummary, run as run_filter};
pub use index::{
    IndexHeader, build as build_index, diff as diff_index, dump as dump_index, dump_minimizers,
    info as index_info, intersect as intersect_index, load_minimizers, union as union_index,
};
pub use minimizers::{DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, decode_u64, decode_u128};

use anyhow::Result;
use std::collections::HashSet;
use std::hash::BuildHasher;
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

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
pub enum MinimizerSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl MinimizerSet {
    pub fn len(&self) -> usize {
        match self {
            MinimizerSet::U64(set) => set.len(),
            MinimizerSet::U128(set) => set.len(),
        }
    }

    pub fn is_u64(&self) -> bool {
        matches!(self, MinimizerSet::U64(_))
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
            _ => panic!("Cannot intersect U64 set with U128 set or vice versa"),
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

    /// Path to JSON summary file
    pub summary_path: Option<&'a PathBuf>,

    /// Deplete mode (remove sequences WITH matches, original deacon behavior)
    pub deplete: bool,

    /// Replace sequence headers with sequential numbers (1, 2, 3...)
    pub rename: bool,

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

impl<'a> FilterConfig<'a> {
    pub fn new(minimizers_path: &'a Path) -> Self {
        Self {
            minimizers_path,
            input_path: "-",
            input2_path: None,
            output_path: None,
            output2_path: None,
            abs_threshold: 2,
            rel_threshold: 0.01,
            prefix_length: 0,
            summary_path: None,
            deplete: false,
            rename: false,
            threads: 0,             // Use all available threads by default
            compression_level: 2,   // Default compression level
            compression_threads: 0, // Auto-calculate as ceil(total/2)
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

    pub fn with_output(mut self, output_path: &'a Path) -> Self {
        self.output_path = Some(output_path);
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

    pub fn with_threads(mut self, threads: u16) -> Self {
        self.threads = threads;
        self
    }

    pub fn with_compression_level(mut self, compression_level: u8) -> Self {
        self.compression_level = compression_level;
        self
    }

    pub fn with_compression_threads(mut self, compression_threads: u16) -> Self {
        self.compression_threads = compression_threads;
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

    /// Number of execution threads (0 = auto)
    pub threads: u16,

    /// Suppress per-sequence progress output
    pub quiet: bool,

    /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
    pub entropy_threshold: f32,
}

impl IndexConfig {
    /// Create a new index configuration with the specified input path
    pub fn new(input_path: PathBuf) -> Self {
        Self {
            input_path: input_path,
            kmer_length: DEFAULT_KMER_LENGTH,
            window_size: DEFAULT_WINDOW_SIZE,
            output_path: None,
            threads: 8,
            quiet: false,
            entropy_threshold: 0.0,
        }
    }

    /// Validate k-mer and window size constraints
    pub fn validate(&self) -> Result<()> {
        let k = self.kmer_length as usize;
        let w = self.window_size as usize;

        // Check constraints: k <= 61, k+w <= 96, k+w even (ensures k odd and k+w-1 odd)
        if k > 61 || k + w > 96 || (k + w) % 2 != 0 {
            return Err(anyhow::anyhow!(
                "Invalid k-w combination: k={}, w={}, k+w={} (constraints: k<=61, k+w<=96, k+w even)",
                k,
                w,
                k + w
            ));
        }

        Ok(())
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
    pub fn with_output(mut self, output_path: PathBuf) -> Self {
        self.output_path = Some(output_path);
        self
    }

    /// Set threads
    pub fn with_threads(mut self, threads: u16) -> Self {
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
