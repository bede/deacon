use crate::FilterConfig;
use crate::index::load_minimizer_hashes_cached;
use crate::minimizers::KmerHasher;
use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use liblzma::write::XzEncoder;
use packed_seq::SeqVec;
use paraseq::Record;
use paraseq::fastx::Reader;
use paraseq::parallel::{PairedParallelProcessor, ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::iter::zip;
use std::sync::Arc;
use std::time::Instant;
use crate::HashSet;
use xxhash_rust;
use zstd::stream::write::Encoder as ZstdEncoder;

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

/// Create a paraseq reader from optional path (stdin if None or "-")
fn create_paraseq_reader(path: Option<&str>) -> Result<Reader<Box<dyn std::io::Read + Send>>> {
    match path {
        None | Some("-") => {
            let stdin_reader = Box::new(std::io::stdin()) as Box<dyn std::io::Read + Send>;
            Reader::new(stdin_reader)
                .map_err(|e| anyhow::anyhow!("Failed to create stdin reader: {}", e))
        }
        Some(p) => {
            let (reader, _format) = niffler::send::from_path(p)
                .map_err(|e| anyhow::anyhow!("Failed to open file {}: {}", p, e))?;
            Reader::new(reader)
                .map_err(|e| anyhow::anyhow!("Failed to create reader for {}: {}", p, e))
        }
    }
}

/// Format a single record into a buffer (FASTA/FASTQ format)
///
/// `seq` is the newline-free sequence corresponding to the record, obtained from `record.seq()`.
fn format_record_to_buffer<R: Record>(
    record: &R,
    seq: &[u8],
    counter: u64,
    rename: bool,
    buffer: &mut Vec<u8>,
) -> Result<()> {
    let is_fasta = record.qual().is_none();

    // Header line
    buffer.write_all(if is_fasta { b">" } else { b"@" })?;
    if rename {
        buffer.extend_from_slice(counter.to_string().as_bytes());
    } else {
        buffer.extend_from_slice(record.id());
    }
    buffer.write_all(b"\n")?;

    // Sequence line
    buffer.extend_from_slice(seq);

    if is_fasta {
        buffer.write_all(b"\n")?;
    } else {
        // FASTQ: plus line and quality
        buffer.write_all(b"\n+\n")?;
        if let Some(qual) = record.qual() {
            buffer.extend_from_slice(qual);
        }
        buffer.write_all(b"\n")?;
    }
    Ok(())
}

/// Validate compression level for the given format
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

// Return a file writer appropriate for the output path extension
fn get_writer(output_path: Option<&std::path::Path>, compression_level: u8) -> Result<BoxedWriter> {
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
        p if p.ends_with(".gz") => {
            validate_compression_level(compression_level, 1, 9, "gzip")?;
            Ok(Box::new(GzEncoder::new(
                buffered_file,
                flate2::Compression::new(compression_level as u32),
            )))
        }
        p if p.ends_with(".zst") => {
            validate_compression_level(compression_level, 1, 22, "zstd")?;
            Ok(Box::new(ZstdEncoder::new(
                buffered_file,
                compression_level as i32,
            )?))
        }
        p if p.ends_with(".xz") => {
            validate_compression_level(compression_level, 0, 9, "xz")?;
            Ok(Box::new(XzEncoder::new(
                buffered_file,
                compression_level as u32,
            )))
        }
        _ => Ok(Box::new(buffered_file)),
    }
}

// JSON summary structure
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
    minimizer_hashes: &'static HashSet,
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

#[derive(Clone, Default)]
struct ProcessingStats {
    total_seqs: u64,
    filtered_seqs: u64,
    total_bp: u64,
    output_bp: u64,
    filtered_bp: u64,
    output_seq_counter: u64,
}

#[derive(Default, Clone)]
pub(crate) struct Buffers {
    pub packed_seq: packed_seq::PackedSeqVec,
    pub ambiguous: packed_seq::BitSeqVec,
    pub positions: Vec<u32>,
    pub hashes: Vec<u64>,
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
        minimizer_hashes: &'static HashSet,
        kmer_length: u8,
        window_size: u8,
        config: &FilterProcessorConfig,
        writer: BoxedWriter,
        writer2: Option<BoxedWriter>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        filtering_start_time: Instant,
    ) -> Self {
        Self {
            minimizer_hashes,
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
            buffers: Buffers::default(),
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

        let seq = self.get_minimizer_positions_and_hashes(seq);
        let positions = &self.buffers.positions;
        let hashes = &self.buffers.hashes;

        let num_minimizers = positions.len();

        // Count distinct minimizer hits and collect matching k-mers
        let mut seen_hits = FxHashSet::default();
        let mut hit_kmers = Vec::new();

        for (&pos, &hash) in zip(positions, hashes) {
            if self.minimizer_hashes.contains(hash) && seen_hits.insert(hash) {
                // Extract the k-mer sequence at this position
                if self.debug {
                    let pos = pos as usize;
                    let kmer = &seq[pos..pos + self.kmer_length as usize];
                    hit_kmers.push(String::from_utf8_lossy(kmer).to_string());
                }
            }
        }

        (
            self.meets_filtering_criteria(seen_hits.len(), num_minimizers),
            seen_hits.len(),
            num_minimizers,
            hit_kmers,
        )
    }

    fn get_minimizer_positions_and_hashes<'s>(&mut self, seq: &'s [u8]) -> &'s [u8] {
        // Apply prefix length limit if specified.
        let seq = if self.prefix_length > 0 && seq.len() > self.prefix_length {
            &seq[..self.prefix_length]
        } else {
            seq
        };
        // Drop trailing newline if present.
        let seq = seq.strip_suffix(b"\n").unwrap_or(seq);

        let Buffers {
            packed_seq,
            ambiguous,
            positions,
            hashes,
        } = &mut self.buffers;

        packed_seq.clear();
        hashes.clear();
        positions.clear();
        ambiguous.clear();

        // Pack the sequence into 2-bit representation.
        packed_seq.push_ascii(seq);
        ambiguous.push_ascii(seq);

        // let mut positions = Vec::new();
        let k = self.kmer_length as usize;
        let w = self.window_size as usize;
        let m = simd_minimizers::canonical_minimizers(k, w)
            .hasher(&self.hasher)
            .run_skip_ambiguous_windows(packed_seq.as_slice(), ambiguous.as_slice(), positions);

        // Hash valid positions
        if self.kmer_length <= 32 {
            hashes.extend(
                m.pos_and_values_u64()
                    .map(|(_pos, val)| xxhash_rust::xxh3::xxh3_64(&val.to_le_bytes())),
            );
        } else {
            hashes.extend(
                m.pos_and_values_u128()
                    .map(|(_pos, val)| xxhash_rust::xxh3::xxh3_64(&val.to_le_bytes())),
            );
        }

        seq
    }

    fn should_keep_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> (bool, usize, usize, Vec<String>) {
        // Count distinct minimizer hits and collect matching k-mers
        let mut seen_hits = FxHashSet::default();
        let mut hit_kmers = Vec::new();
        let mut num_minimizers = 0;

        for seq in [seq1, seq2] {
            let seq = self.get_minimizer_positions_and_hashes(seq);
            let positions = &self.buffers.positions;
            let hashes = &self.buffers.hashes;

            num_minimizers += positions.len();

            for (&pos, &hash) in zip(positions, hashes) {
                if self.minimizer_hashes.contains(hash) && seen_hits.insert(hash) {
                    // Extract the k-mer sequence at this position
                    if self.debug {
                        let pos = pos as usize;
                        let kmer = &seq[pos..pos + self.kmer_length as usize];
                        hit_kmers.push(String::from_utf8_lossy(kmer).to_string());
                    }
                }
            }
        }

        (
            self.meets_filtering_criteria(seen_hits.len(), num_minimizers),
            seen_hits.len(),
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

    // Configure thread pool if nonzero
    if config.threads > 0 {
        // error is OK here when we initialize a 2nd time in server mode.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
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
        options.push(format!("threads={}", config.threads));
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

    // Load minimizer hashes and parse header
    let (minimizer_hashes, header) = load_minimizer_hashes_cached(&config.minimizers_path)?;

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
    let writer = get_writer(config.output_path, config.compression_level)?;
    let writer2 = if let (Some(output2), Some(_)) = (config.output2_path, config.input2_path) {
        Some(get_writer(
            Some(std::path::Path::new(output2)),
            config.compression_level,
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
        pb.set_message("Filtering");
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
        minimizer_hashes,
        kmer_length,
        window_size,
        &processor_config,
        writer,
        writer2,
        spinner.clone(),
        filtering_start_time,
    );

    // Process based on input type
    let num_threads = if config.threads == 0 {
        rayon::current_num_threads()
    } else {
        config.threads
    };

    if paired_stdin {
        // Interleaved paired from stdin - use native interleaved processor
        let reader = create_paraseq_reader(Some("-"))?;
        reader.process_parallel_interleaved(&mut processor, num_threads)?;
    } else if let Some(input2_path) = config.input2_path {
        // Paired files
        let r1_reader = create_paraseq_reader(Some(config.input_path))?;
        let r2_reader = create_paraseq_reader(Some(input2_path))?;
        r1_reader.process_parallel_paired(r2_reader, &mut processor, num_threads)?;
    } else {
        // Single file or stdin
        let reader = create_paraseq_reader(Some(config.input_path))?;
        reader.process_parallel(&mut processor, num_threads)?;
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
            eprintln!("Summary saved to {:?}", summary_file);
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
