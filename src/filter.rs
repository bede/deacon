use crate::index::load_minimizers_cached;
use crate::minimizers::{KmerHasher, decode_u64, decode_u128};
use crate::{FilterConfig, MinimizerSet};
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use packed_seq::{PackedNSeqVec, SeqVec, u32x8};
use paraseq::Record;
use paraseq::fastx::Reader;
use paraseq::parallel::{PairedParallelProcessor, ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::sync::Arc;
use std::time::Instant;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer
const DEFAULT_BUFFER_SIZE: usize = 64 * 1024;

type BoxedWriter = Box<dyn Write + Send>;

/// Config for FilterProcessor
struct FilterProcessorConfig {
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    deplete: bool,
    rename: bool,
    debug: bool,
}

/// Check input file path(s) exist
fn check_input_paths(config: &FilterConfig) -> Result<()> {
    if !config.minimizers_path.exists() {
        return Err(anyhow::anyhow!(
            "Index file does not exist: {}",
            config.minimizers_path.display()
        ));
    }

    // Skip stdin case
    if config.input_path != "-" && !std::path::Path::new(config.input_path).exists() {
        return Err(anyhow::anyhow!(
            "Input file does not exist: {}",
            config.input_path
        ));
    }

    if let Some(input2_path) = config.input2_path {
        if input2_path != "-" && !std::path::Path::new(input2_path).exists() {
            return Err(anyhow::anyhow!(
                "Second input file does not exist: {}",
                input2_path
            ));
        }
    }

    Ok(())
}

/// Check if file metadata len < 5 (catches empty uncompressed files only)
fn is_empty_file(path: &str) -> Result<bool> {
    if path == "-" {
        return Ok(false); // Can't check stdin
    }
    let metadata = std::fs::metadata(path)
        .map_err(|e| anyhow::anyhow!("Failed to read file metadata {}: {}", path, e))?;
    Ok(metadata.len() < 5)
}

/// Error thrown when paraseq reads empty compressed file
fn is_empty_input_error(err: &anyhow::Error) -> bool {
    err.to_string().contains("failed to fill whole buffer")
}

/// Create a paraseq reader from optional path (stdin if None or "-")
fn create_paraseq_reader(path: Option<&str>) -> Result<Reader<Box<dyn std::io::Read + Send>>> {
    match path {
        None | Some("-") => {
            let stdin_reader = Box::new(std::io::stdin()) as Box<dyn std::io::Read + Send>;
            Reader::new(stdin_reader)
                .map_err(|e| anyhow::anyhow!("Failed to create stdin reader: {}", e))
        }
        Some(p) => {
            // Use paraseq's from_path for files (internally uses niffler for compression detection)
            Reader::from_path(p).map_err(|e| anyhow::anyhow!("Failed to open file {}: {}", p, e))
        }
    }
}

/// Format a single record into a buffer (FASTA/FASTQ format)
/// `seq` is the newline-stripped sequence corresponding to the record from `record.seq()`.
fn format_record_to_buffer<R: Record>(
    record: &R,
    seq: &[u8],
    counter: u64,
    rename: bool,
    buffer: &mut Vec<u8>,
) -> Result<()> {
    let is_fasta = record.qual().is_none();

    // Header
    buffer.write_all(if is_fasta { b">" } else { b"@" })?;
    if rename {
        buffer.extend_from_slice(counter.to_string().as_bytes());
    } else {
        buffer.extend_from_slice(record.id());
    }
    buffer.write_all(b"\n")?;

    // Sequence
    buffer.extend_from_slice(seq);

    if is_fasta {
        buffer.write_all(b"\n")?;
    } else {
        // FASTQ: plus and qual lines
        buffer.write_all(b"\n+\n")?;
        if let Some(qual) = record.qual() {
            buffer.extend_from_slice(qual);
        }
        buffer.write_all(b"\n")?;
    }
    Ok(())
}

/// Validate compression level for the given format
#[cfg(feature = "compression")]
fn validate_compression_level(level: u8, min: u8, max: u8, format: &str) -> Result<()> {
    if level < min || level > max {
        Err(anyhow::anyhow!(
            "Invalid {} compression level {}. Must be between {} and {}.",
            format,
            level,
            min,
            max
        ))
    } else {
        Ok(())
    }
}

/// Check if a path requires gzip compression
fn is_compressed_output(path: Option<&std::path::Path>) -> bool {
    path.map(|p| p.to_string_lossy().ends_with(".gz"))
        .unwrap_or(false)
}

/// Number of compressed outputs from config (0, 1, or 2)
fn count_compressed_outputs(config: &FilterConfig) -> u8 {
    let mut count = 0;
    if is_compressed_output(config.output_path) {
        count += 1;
    }
    if let Some(output2) = config.output2_path {
        if is_compressed_output(Some(std::path::Path::new(output2))) {
            count += 1;
        }
    }
    count
}

// Return a suitable writer for the output path extension
#[cfg_attr(not(feature = "compression"), allow(unused_variables))]
fn get_writer(
    output_path: Option<&std::path::Path>,
    compression_level: u8,
    compression_threads: usize,
) -> Result<BoxedWriter> {
    let Some(path) = output_path else {
        return Ok(Box::new(BufWriter::with_capacity(
            OUTPUT_BUFFER_SIZE,
            io::stdout(),
        )));
    };

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)
        .context(format!("Failed to create output file: {}", path.display()))?;

    let buffered_file = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

    match path.to_string_lossy().as_ref() {
        #[cfg(feature = "compression")]
        p if p.ends_with(".gz") => {
            validate_compression_level(compression_level, 1, 9, "gzip")?;
            use gzp::deflate::Gzip;
            use gzp::par::compress::ParCompressBuilder;

            // Use the calculated number of threads for gzip compression
            let writer = ParCompressBuilder::<Gzip>::new()
                .compression_level(gzp::Compression::new(compression_level as u32))
                .buffer_size(1024 * 1024) // 1MB buf
                .unwrap()
                .num_threads(compression_threads)
                .unwrap()
                .from_writer(buffered_file);
            Ok(Box::new(writer))
        }
        #[cfg(feature = "compression")]
        p if p.ends_with(".zst") => {
            validate_compression_level(compression_level, 1, 22, "zstd")?;
            Ok(Box::new(zstd::stream::write::Encoder::new(
                buffered_file,
                compression_level as i32,
            )?))
        }
        #[cfg(feature = "compression")]
        p if p.ends_with(".xz") => {
            validate_compression_level(compression_level, 0, 9, "xz")?;
            Ok(Box::new(liblzma::write::XzEncoder::new(
                buffered_file,
                compression_level as u32,
            )))
        }
        _ => Ok(Box::new(buffered_file)),
    }
}

// JSON summary struct
#[derive(Serialize, Deserialize)]
pub struct FilterSummary {
    version: String,
    index: String,
    input: String,
    input2: Option<String>,
    output: String,
    output2: Option<String>,
    k: u8,
    w: u8,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    deplete: bool,
    rename: bool,
    seqs_in: u64,
    seqs_out: u64,
    seqs_out_proportion: f64,
    seqs_removed: u64,
    seqs_removed_proportion: f64,
    bp_in: u64,
    bp_out: u64,
    bp_out_proportion: f64,
    bp_removed: u64,
    bp_removed_proportion: f64,
    time: f64,
    seqs_per_second: u64,
    bp_per_second: u64,
    seqs_per_second_total: u64,
    bp_per_second_total: u64,
}

#[derive(Clone)]
struct FilterProcessor {
    // Minimizer matching parameters
    minimizers: &'static MinimizerSet,
    kmer_length: u8,
    window_size: u8,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    deplete: bool,
    rename: bool,
    debug: bool,

    hasher: KmerHasher,

    // Local buffers
    local_buffer: Vec<u8>,
    local_buffer2: Vec<u8>, // Second buffer for paired output
    local_stats: ProcessingStats,
    buffers: Buffers,

    // Global state
    global_writer: Arc<Mutex<BoxedWriter>>,
    global_writer2: Option<Arc<Mutex<BoxedWriter>>>,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    filtering_start_time: Instant,
}

#[derive(Clone, Default, Debug)]
pub(crate) struct ProcessingStats {
    pub total_seqs: u64,
    filtered_seqs: u64,
    pub total_bp: u64,
    output_bp: u64,
    filtered_bp: u64,
    output_seq_counter: u64,
    pub last_reported: u64,
}

#[derive(Clone)]
pub(crate) struct Buffers {
    pub packed_nseq: PackedNSeqVec,
    pub positions: Vec<u32>,
    pub minimizers: crate::MinimizerVec,
    pub cache: (simd_minimizers::Cache, Vec<u32x8>, Vec<u32x8>),
}

impl Buffers {
    pub fn new_u64() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            minimizers: crate::MinimizerVec::U64(Vec::new()),
            cache: Default::default(),
        }
    }

    pub fn new_u128() -> Self {
        Self {
            packed_nseq: PackedNSeqVec {
                seq: Default::default(),
                ambiguous: Default::default(),
            },
            positions: Default::default(),
            minimizers: crate::MinimizerVec::U128(Vec::new()),
            cache: Default::default(),
        }
    }
}

impl FilterProcessor {
    /// Calculate required hits based on absolute and relative thresholds
    fn calculate_required_hits(&self, total_minimizers: usize) -> usize {
        let abs_required = self.abs_threshold;
        let rel_required = if total_minimizers == 0 {
            0
        } else {
            ((self.rel_threshold * total_minimizers as f64).round() as usize).max(1)
        };
        abs_required.max(rel_required)
    }

    /// Check if sequence meets filtering criteria
    fn meets_filtering_criteria(&self, hit_count: usize, total_minimizers: usize) -> bool {
        let required = self.calculate_required_hits(total_minimizers);
        if self.deplete {
            hit_count < required
        } else {
            hit_count >= required
        }
    }
    fn new(
        minimizers: &'static MinimizerSet,
        kmer_length: u8,
        window_size: u8,
        config: &FilterProcessorConfig,
        writer: BoxedWriter,
        writer2: Option<BoxedWriter>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        filtering_start_time: Instant,
    ) -> Self {
        Self {
            minimizers,
            kmer_length,
            window_size,
            abs_threshold: config.abs_threshold,
            rel_threshold: config.rel_threshold,
            prefix_length: config.prefix_length,
            deplete: config.deplete,
            rename: config.rename,
            debug: config.debug,
            hasher: KmerHasher::new(kmer_length as usize),
            local_buffer: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
            local_buffer2: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
            local_stats: ProcessingStats::default(),
            buffers: if kmer_length <= 32 {
                Buffers::new_u64()
            } else {
                Buffers::new_u128()
            },
            global_writer: Arc::new(Mutex::new(writer)),
            global_writer2: writer2.map(|w| Arc::new(Mutex::new(w))),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            spinner,
            filtering_start_time,
        }
    }

    fn should_keep_sequence(&mut self, seq: &[u8]) -> (bool, usize, usize, Vec<String>) {
        if seq.len() < self.kmer_length as usize {
            return (self.deplete, 0, 0, Vec::new()); // If too short, keep if in deplete mode
        }

        self.get_minimizer_positions_and_values(seq);
        let minimizers = &self.buffers.minimizers;
        let num_minimizers = minimizers.len();

        // Count distinct minimizer hits based on variant
        let (hit_count, hit_kmers) = match (minimizers, self.minimizers) {
            (crate::MinimizerVec::U64(vec), MinimizerSet::U64(set)) => {
                let mut seen_hits = crate::RapidHashSet::default();
                let mut hit_kmers = Vec::new();
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen_hits.insert(minimizer) {
                        if self.debug {
                            let kmer = decode_u64(minimizer, self.kmer_length);
                            hit_kmers.push(String::from_utf8_lossy(&kmer).to_string());
                        }
                    }
                }
                (seen_hits.len(), hit_kmers)
            }
            (crate::MinimizerVec::U128(vec), MinimizerSet::U128(set)) => {
                let mut seen_hits = crate::RapidHashSet::default();
                let mut hit_kmers = Vec::new();
                for &minimizer in vec {
                    if set.contains(&minimizer) && seen_hits.insert(minimizer) {
                        if self.debug {
                            let kmer = decode_u128(minimizer, self.kmer_length);
                            hit_kmers.push(String::from_utf8_lossy(&kmer).to_string());
                        }
                    }
                }
                (seen_hits.len(), hit_kmers)
            }
            _ => panic!("Mismatch between MinimizerVec and MinimizerSet types"),
        };

        (
            self.meets_filtering_criteria(hit_count, num_minimizers),
            hit_count,
            num_minimizers,
            hit_kmers,
        )
    }

    fn get_minimizer_positions_and_values<'s>(&mut self, seq: &'s [u8]) -> &'s [u8] {
        // Apply prefix length limit if specified.
        let seq = if self.prefix_length > 0 && seq.len() > self.prefix_length {
            &seq[..self.prefix_length]
        } else {
            seq
        };
        // Drop trailing newline if present.
        let seq = seq.strip_suffix(b"\n").unwrap_or(seq);

        let Buffers {
            packed_nseq,
            positions,
            minimizers,
            cache,
        } = &mut self.buffers;

        packed_nseq.seq.clear();
        packed_nseq.ambiguous.clear();
        minimizers.clear();
        positions.clear();

        // Pack the sequence into 2-bit representation.
        packed_nseq.seq.push_ascii(seq);
        packed_nseq.ambiguous.push_ascii(seq);

        // let mut positions = Vec::new();
        let k = self.kmer_length as usize;
        let w = self.window_size as usize;
        let m = simd_minimizers::canonical_minimizers(k, w)
            .hasher(&self.hasher)
            .run_skip_ambiguous_windows_with_buf(packed_nseq.as_slice(), positions, cache);

        // Store k-mer values directly based on variant
        match minimizers {
            crate::MinimizerVec::U64(vec) => {
                vec.extend(m.pos_and_values_u64().map(|(_pos, val)| val));
            }
            crate::MinimizerVec::U128(vec) => {
                vec.extend(m.pos_and_values_u128().map(|(_pos, val)| val));
            }
        }

        seq
    }

    fn should_keep_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> (bool, usize, usize, Vec<String>) {
        // Process both sequences and count distinct hits
        let (hit_count, num_minimizers, hit_kmers) = match self.minimizers {
            MinimizerSet::U64(set) => {
                let mut seen_hits = crate::RapidHashSet::default();
                let mut hit_kmers = Vec::new();
                let mut num_minimizers = 0;

                for seq in [seq1, seq2] {
                    self.get_minimizer_positions_and_values(seq);
                    if let crate::MinimizerVec::U64(vec) = &self.buffers.minimizers {
                        num_minimizers += vec.len();
                        for &minimizer in vec {
                            if set.contains(&minimizer) && seen_hits.insert(minimizer) {
                                if self.debug {
                                    let kmer = decode_u64(minimizer, self.kmer_length);
                                    hit_kmers.push(String::from_utf8_lossy(&kmer).to_string());
                                }
                            }
                        }
                    }
                }
                (seen_hits.len(), num_minimizers, hit_kmers)
            }
            MinimizerSet::U128(set) => {
                let mut seen_hits = crate::RapidHashSet::default();
                let mut hit_kmers = Vec::new();
                let mut num_minimizers = 0;

                for seq in [seq1, seq2] {
                    self.get_minimizer_positions_and_values(seq);
                    if let crate::MinimizerVec::U128(vec) = &self.buffers.minimizers {
                        num_minimizers += vec.len();
                        for &minimizer in vec {
                            if set.contains(&minimizer) && seen_hits.insert(minimizer) {
                                if self.debug {
                                    let kmer = decode_u128(minimizer, self.kmer_length);
                                    hit_kmers.push(String::from_utf8_lossy(&kmer).to_string());
                                }
                            }
                        }
                    }
                }
                (seen_hits.len(), num_minimizers, hit_kmers)
            }
        };

        (
            self.meets_filtering_criteria(hit_count, num_minimizers),
            hit_count,
            num_minimizers,
            hit_kmers,
        )
    }

    fn write_record<Rf: Record>(&mut self, record: &Rf, seq: &[u8]) -> Result<()> {
        self.local_stats.output_seq_counter += 1;
        format_record_to_buffer(
            record,
            seq,
            self.local_stats.output_seq_counter,
            self.rename,
            &mut self.local_buffer,
        )
    }

    fn write_record_to_buffer2<Rf: Record>(&mut self, record: &Rf, seq: &[u8]) -> Result<()> {
        self.local_stats.output_seq_counter += 1;
        format_record_to_buffer(
            record,
            seq,
            self.local_stats.output_seq_counter,
            self.rename,
            &mut self.local_buffer2,
        )
    }

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.filtering_start_time.elapsed();
            let seqs_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
            let mbp_per_sec = bp_per_sec / 1_000_000.0;

            let output_seqs = stats.total_seqs - stats.filtered_seqs;
            let output_proportion = if stats.total_seqs > 0 {
                output_seqs as f64 / stats.total_seqs as f64
            } else {
                0.0
            };

            let output_bp_proportion = if stats.total_bp > 0 {
                stats.output_bp as f64 / stats.total_bp as f64
            } else {
                0.0
            };

            spinner.lock().set_message(format!(
                "Retained {}/{} sequences ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
                output_seqs,
                stats.total_seqs,
                output_proportion * 100.0,
                stats.output_bp,
                stats.total_bp,
                output_bp_proportion * 100.0,
                seqs_per_sec,
                mbp_per_sec
            ));
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for FilterProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        let (should_keep, hit_count, total_minimizers, hit_kmers) = self.should_keep_sequence(&seq);

        // Show debug info for sequences with hits
        if self.debug {
            eprintln!(
                "DEBUG: {} hits={}/{} keep={} kmers=[{}]",
                String::from_utf8_lossy(record.id()),
                hit_count,
                total_minimizers,
                should_keep,
                hit_kmers.join(",")
            );
        }

        if should_keep {
            self.local_stats.output_bp += seq.len() as u64;
            self.write_record(&record, &seq)?;
        } else {
            self.local_stats.filtered_seqs += 1;
            self.local_stats.filtered_bp += seq.len() as u64;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Write buffer to output
        if !self.local_buffer.is_empty() {
            let mut global_writer = self.global_writer.lock();
            global_writer.write_all(&self.local_buffer)?;
            global_writer.flush()?;
        }

        // Clear buffer after releasing the lock
        self.local_buffer.clear();

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.filtered_seqs += self.local_stats.filtered_seqs;
            stats.total_bp += self.local_stats.total_bp;
            stats.output_bp += self.local_stats.output_bp;
            stats.filtered_bp += self.local_stats.filtered_bp;
            stats.output_seq_counter += self.local_stats.output_seq_counter;
        }

        // Update spinner
        self.update_spinner();

        // Reset local stats
        self.local_stats = ProcessingStats::default();

        Ok(())
    }
}

impl<Rf: Record> PairedParallelProcessor<Rf> for FilterProcessor {
    fn process_record_pair(&mut self, record1: Rf, record2: Rf) -> paraseq::parallel::Result<()> {
        let seq1 = record1.seq();
        let seq2 = record2.seq();

        self.local_stats.total_seqs += 2;
        self.local_stats.total_bp += (seq1.len() + seq2.len()) as u64;

        let (should_keep, hit_count, total_minimizers, hit_kmers) =
            self.should_keep_pair(&seq1, &seq2);

        // Debug info for interleaved pairs
        if self.debug && hit_count > 0 {
            eprintln!(
                "DEBUG: {}/{} hits={}/{} keep={} kmers=[{}]",
                String::from_utf8_lossy(record1.id()),
                String::from_utf8_lossy(record2.id()),
                hit_count,
                total_minimizers,
                should_keep,
                hit_kmers.join(",")
            );
        }

        if should_keep {
            self.local_stats.output_bp += (seq1.len() + seq2.len()) as u64;

            // Write to appropriate writers
            if self.global_writer2.is_some() {
                // Separate outputs
                self.write_record(&record1, &seq1)?;
                self.write_record_to_buffer2(&record2, &seq2)?;
            } else {
                // Interleaved output
                self.write_record(&record1, &seq1)?;
                self.write_record(&record2, &seq2)?;
            }
        } else {
            self.local_stats.filtered_seqs += 2;
            self.local_stats.filtered_bp += (seq1.len() + seq2.len()) as u64;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        if let Some(ref writer2) = self.global_writer2 {
            // Atomic paired batch writing
            if !self.local_buffer.is_empty() || !self.local_buffer2.is_empty() {
                let mut writer1 = self.global_writer.lock();
                let mut writer2 = writer2.lock();

                writer1.write_all(&self.local_buffer)?;
                writer1.flush()?;
                writer2.write_all(&self.local_buffer2)?;
                writer2.flush()?;
            }
        } else {
            // Interleaved output
            if !self.local_buffer.is_empty() {
                let mut writer = self.global_writer.lock();
                writer.write_all(&self.local_buffer)?;
                writer.flush()?;
            }
        }

        // Clear buffer after releasing the lock for better performance
        self.local_buffer.clear();
        self.local_buffer2.clear();

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.filtered_seqs += self.local_stats.filtered_seqs;
            stats.total_bp += self.local_stats.total_bp;
            stats.output_bp += self.local_stats.output_bp;
            stats.filtered_bp += self.local_stats.filtered_bp;
            stats.output_seq_counter += self.local_stats.output_seq_counter;
        }

        // Update spinner
        self.update_spinner();

        // Reset local stats
        self.local_stats = ProcessingStats::default();

        Ok(())
    }
}

pub fn run(config: &FilterConfig) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    let tool_version = format!("deacon {}", version);

    // Enable quiet mode when debug enabled
    let quiet = config.quiet || config.debug;

    let total_threads = if config.threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        config.threads as usize
    };

    let compressed_output_count = count_compressed_outputs(config);

    // Allocate threads between filtering (rayon) and compression (gzp).
    // Rayon pool can only be initialized once, so calculate before build_global().
    let (filtering_threads, compression_threads_per_output) = if compressed_output_count > 0 {
        let compression_threads_total = if config.compression_threads > 0 {
            config.compression_threads as usize
        } else {
            (total_threads + 1) / 2 // Auto: ceil(total_threads / 2)
        };
        let filtering_threads = total_threads
            .saturating_sub(compression_threads_total)
            .max(1);
        let output_count = compressed_output_count as usize;
        let threads_per_output =
            ((compression_threads_total + output_count - 1) / output_count).max(1);
        (filtering_threads, threads_per_output)
    } else {
        (total_threads, 0)
    };

    if filtering_threads > 0 {
        // error is OK here when we initialize a 2nd time in server mode.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(filtering_threads)
            .build_global()
            .context("Failed to initialize thread pool");
    }

    let mode = if config.deplete { "deplete" } else { "search" };

    let mut input_type = String::new();
    let mut options = Vec::<String>::new();
    let paired_stdin = config.input_path == "-"
        && config.input2_path.is_some()
        && config.input2_path.unwrap() == "-";
    if paired_stdin {
        input_type.push_str("interleaved");
    } else if config.input2_path.is_some() {
        input_type.push_str("paired");
    } else {
        input_type.push_str("single");
    }
    options.push(format!(
        "abs_threshold={}, rel_threshold={}",
        config.abs_threshold, config.rel_threshold
    ));
    if config.prefix_length > 0 {
        options.push(format!("prefix_length={}", config.prefix_length));
    }
    if config.rename {
        options.push("rename".to_string());
    }
    if config.threads > 0 {
        let threads_str = if compressed_output_count > 0 {
            let compression_total =
                compressed_output_count as usize * compression_threads_per_output;
            format!(
                "threads={}({}f+{}c)",
                config.threads, filtering_threads, compression_total
            )
        } else {
            format!("threads={}", config.threads)
        };
        options.push(threads_str);
    }

    if !quiet {
        eprintln!(
            "Deacon v{}; mode: {}; input: {}; options: {}",
            version,
            mode,
            input_type,
            options.join(", ")
        );
    }

    // Check input files exist before making user wait for index loading
    check_input_paths(config)?;

    // Load minimizers and parse header
    let (minimizers, header) = load_minimizers_cached(&config.minimizers_path)?;

    let kmer_length = header.kmer_length();
    let window_size = header.window_size();

    let load_time = start_time.elapsed();
    if !quiet {
        eprintln!(
            "Loaded index (k={}, w={}) in {:.2?}",
            kmer_length, window_size, load_time
        );
    }

    // Create appropriate writer(s) based on output path(s)
    let writer = get_writer(
        config.output_path,
        config.compression_level,
        compression_threads_per_output,
    )?;
    let writer2 = if let (Some(output2), Some(_)) = (config.output2_path, config.input2_path) {
        Some(get_writer(
            Some(std::path::Path::new(output2)),
            config.compression_level,
            compression_threads_per_output,
        )?)
    } else {
        None
    };

    // Progress bar setup if not quiet
    let spinner = if !quiet {
        let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{msg}")?,
        );
        Some(Arc::new(Mutex::new(pb)))
    } else {
        None
    };

    // Start timer for rate calculation
    let filtering_start_time = Instant::now();

    // Create processor
    let processor_config = FilterProcessorConfig {
        abs_threshold: config.abs_threshold,
        rel_threshold: config.rel_threshold,
        prefix_length: config.prefix_length,
        deplete: config.deplete,
        rename: config.rename,
        debug: config.debug,
    };
    let mut processor = FilterProcessor::new(
        minimizers,
        kmer_length,
        window_size,
        &processor_config,
        writer,
        writer2,
        spinner.clone(),
        filtering_start_time,
    );

    // Check for empty files via metadata (fast path for uncompressed files <5 bytes)
    let input1_empty = is_empty_file(config.input_path)?;
    let input2_empty = config
        .input2_path
        .map(|p| is_empty_file(p))
        .transpose()?
        .unwrap_or(false);

    // Process based on input type - use filtering threads (already calculated above)
    let num_threads = filtering_threads;

    if paired_stdin {
        // Interleaved paired stdin - use interleaved processor
        let reader = create_paraseq_reader(Some("-"))?;
        reader.process_parallel_interleaved(&mut processor, num_threads)?;
    } else if let Some(input2_path) = config.input2_path {
        // Paired files - both must be empty or both non-empty
        if input1_empty && input2_empty {
            if !quiet {
                eprintln!("Empty input file(s) detected");
            }
        } else if input1_empty || input2_empty {
            return Err(anyhow::anyhow!(
                "One paired file is empty but the other is not"
            ));
        } else {
            // Try to create readers, catching empty compressed files
            let r1_result = create_paraseq_reader(Some(config.input_path));
            let r2_result = create_paraseq_reader(Some(input2_path));

            match (r1_result, r2_result) {
                (Ok(r1), Ok(r2)) => {
                    r1.process_parallel_paired(r2, &mut processor, num_threads)?;
                }
                (Err(e1), Err(e2)) if is_empty_input_error(&e1) && is_empty_input_error(&e2) => {
                    if !quiet {
                        eprintln!("Empty input file(s) detected");
                    }
                }
                (Err(e), _) if is_empty_input_error(&e) => {
                    return Err(anyhow::anyhow!(
                        "First paired file appears empty while second is not"
                    ));
                }
                (_, Err(e)) if is_empty_input_error(&e) => {
                    return Err(anyhow::anyhow!(
                        "Second paired file appears empty while first is not"
                    ));
                }
                (Err(e), _) => return Err(e),
                (_, Err(e)) => return Err(e),
            }
        }
    } else {
        // Single file or stdin
        if input1_empty {
            if !quiet {
                eprintln!("Empty input file(s) detected");
            }
        } else {
            // Try to create reader, catching empty compressed files
            match create_paraseq_reader(Some(config.input_path)) {
                Ok(reader) => {
                    reader.process_parallel(&mut processor, num_threads)?;
                }
                Err(e) if is_empty_input_error(&e) => {
                    if !quiet {
                        eprintln!("Empty input file(s) detected");
                    }
                }
                Err(e) => return Err(e),
            }
        }
    }

    let final_stats = processor.global_stats.lock();
    let total_seqs = final_stats.total_seqs;
    let filtered_seqs = final_stats.filtered_seqs;
    let total_bp = final_stats.total_bp;
    let output_bp = final_stats.output_bp;
    let filtered_bp = final_stats.filtered_bp;

    drop(final_stats); // Release lock

    // Flush writers - they should auto-flush on drop
    drop(processor.global_writer);
    if let Some(w2) = processor.global_writer2 {
        drop(w2);
    }

    let total_time = start_time.elapsed();
    let filtering_time = filtering_start_time.elapsed();

    // Based on filtering time excluding index loading
    let seqs_per_sec = total_seqs as f64 / filtering_time.as_secs_f64();
    let bp_per_sec = total_bp as f64 / filtering_time.as_secs_f64();
    let mbp_per_sec = bp_per_sec / 1_000_000.0;

    // Based on total time, including index loading
    let seqs_per_sec_total = total_seqs as f64 / total_time.as_secs_f64();
    let bp_per_sec_total = total_bp as f64 / total_time.as_secs_f64();

    // Calculate proportions
    let filtered_proportion = if total_seqs > 0 {
        filtered_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    let filtered_bp_proportion = if total_bp > 0 {
        filtered_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    let output_seqs = total_seqs - filtered_seqs;
    let output_seq_proportion = if total_seqs > 0 {
        output_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    let output_bp_proportion = if total_bp > 0 {
        output_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    // Finish and clear spinner - disable it completely
    if let Some(ref spinner) = spinner {
        let pb = spinner.lock();
        pb.finish_with_message("");
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }

    if !quiet {
        eprintln!(
            "Retained {}/{} sequences ({:.3}%), {}/{} bp ({:.3}%) in {:.2?}. {:.0} seqs/s ({:.1} Mbp/s)",
            output_seqs,
            total_seqs,
            output_seq_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0,
            total_time,
            seqs_per_sec,
            mbp_per_sec
        );
    }

    // Build and write JSON summary if path provided
    if let Some(summary_file) = config.summary_path {
        let seqs_out = total_seqs - filtered_seqs;

        let summary = FilterSummary {
            version: tool_version,
            index: config.minimizers_path.to_string_lossy().to_string(),
            input: config.input_path.to_string(),
            input2: config.input2_path.map(|s| s.to_string()),
            output: config
                .output_path
                .map_or("-".to_string(), |p| p.display().to_string()),
            output2: config.output2_path.map(|s| s.to_string()),
            k: kmer_length,
            w: window_size,
            abs_threshold: config.abs_threshold,
            rel_threshold: config.rel_threshold,
            prefix_length: config.prefix_length,
            deplete: config.deplete,
            rename: config.rename,
            seqs_in: total_seqs as u64,
            seqs_out: seqs_out as u64,
            seqs_out_proportion: output_seq_proportion,
            seqs_removed: filtered_seqs as u64,
            seqs_removed_proportion: filtered_proportion,
            bp_in: total_bp as u64,
            bp_out: output_bp as u64,
            bp_out_proportion: output_bp_proportion,
            bp_removed: filtered_bp as u64,
            bp_removed_proportion: filtered_bp_proportion,
            time: total_time.as_secs_f64(),
            seqs_per_second: seqs_per_sec as u64,
            bp_per_second: bp_per_sec as u64,
            seqs_per_second_total: seqs_per_sec_total as u64,
            bp_per_second_total: bp_per_sec_total as u64,
        };

        let file = File::create(summary_file)
            .context(format!("Failed to create summary: {:?}", summary_file))?;
        let writer = BufWriter::new(file);

        serde_json::to_writer_pretty(writer, &summary).context("Failed to write summary")?;

        if !quiet {
            eprintln!("Filter summary saved to {:?}", summary_file);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_summary() {
        let summary = FilterSummary {
            version: "deacon 0.1.0".to_string(),
            index: "test.idx".to_string(),
            input: "test.fastq".to_string(),
            input2: Some("test2.fastq".to_string()),
            output: "output.fastq".to_string(),
            output2: Some("output2.fastq".to_string()),
            k: 31,
            w: 21,
            abs_threshold: 1,
            rel_threshold: 0.01,
            prefix_length: 0,
            deplete: false,
            rename: false,
            seqs_in: 100,
            seqs_out: 90,
            seqs_out_proportion: 0.9,
            seqs_removed: 10,
            seqs_removed_proportion: 0.1,
            bp_in: 10000,
            bp_out: 9000,
            bp_out_proportion: 0.9,
            bp_removed: 1000,
            bp_removed_proportion: 0.1,
            time: 1.5,
            seqs_per_second: 66,
            bp_per_second: 6666,
            seqs_per_second_total: 60,
            bp_per_second_total: 6000,
        };

        let json = serde_json::to_string(&summary).unwrap();
        let parsed: FilterSummary = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed.version, "deacon 0.1.0");
        assert_eq!(parsed.seqs_in, 100);
        assert_eq!(parsed.seqs_removed_proportion, 0.1);
        assert_eq!(parsed.seqs_out_proportion, 0.9);
        assert_eq!(parsed.bp_out_proportion, 0.9);
        assert_eq!(parsed.input, "test.fastq");
        assert_eq!(parsed.input2, Some("test2.fastq".to_string()));
        assert_eq!(parsed.output, "output.fastq");
        assert_eq!(parsed.output2, Some("output2.fastq".to_string()));
    }
}
