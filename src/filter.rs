use crate::index::load_minimizer_hashes;
use crate::minimizers::fill_minimizer_hashes;
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use needletail::parse_fastx_stdin;
use needletail::parser::{Format, SequenceRecord};
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use zstd::stream::write::Encoder as ZstdEncoder;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer

trait FastxWriter: Write {
    fn flush_all(&mut self) -> io::Result<()>;
}

struct StandardWriter<W: Write>(W);

impl<W: Write> Write for StandardWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

impl<W: Write> FastxWriter for StandardWriter<W> {
    fn flush_all(&mut self) -> io::Result<()> {
        self.flush()
    }
}

struct GzipWriter<W: Write>(GzEncoder<W>);

impl<W: Write> Write for GzipWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

impl<W: Write> FastxWriter for GzipWriter<W> {
    fn flush_all(&mut self) -> io::Result<()> {
        self.flush()?;
        self.0.try_finish()?;
        Ok(())
    }
}

struct ZstdWriter<W: Write> {
    encoder: Option<ZstdEncoder<'static, W>>,
}

impl<W: Write> ZstdWriter<W> {
    fn new(writer: W, compression_level: i32) -> Result<Self> {
        let encoder =
            ZstdEncoder::new(writer, compression_level).context("Failed to create zstd encoder")?;
        Ok(Self {
            encoder: Some(encoder),
        })
    }
}

impl<W: Write> Write for ZstdWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if let Some(encoder) = &mut self.encoder {
            encoder.write(buf)
        } else {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been closed",
            ))
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Some(encoder) = &mut self.encoder {
            encoder.flush()
        } else {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been closed",
            ))
        }
    }
}

impl<W: Write> FastxWriter for ZstdWriter<W> {
    fn flush_all(&mut self) -> io::Result<()> {
        if let Some(encoder) = self.encoder.take() {
            // Take ownership of the encoder
            encoder.finish()?;
        }
        Ok(())
    }
}

// Return a file writer appropriate for the output path extension
fn get_writer(output_path: &str) -> Result<Box<dyn FastxWriter>> {
    if output_path == "-" {
        // Write to stdout
        let stdout = io::stdout();
        let writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, stdout);
        Ok(Box::new(StandardWriter(writer)))
    } else {
        // Write to file with extension-appropriate encoder
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(output_path)
            .context(format!("Failed to create output file: {}", output_path))?;

        let buffered_file = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

        if output_path.ends_with(".gz") {
            // Use gzip
            let encoder = GzEncoder::new(buffered_file, Compression::fast());
            Ok(Box::new(GzipWriter(encoder)))
        } else if output_path.ends_with(".zst") {
            // Use zstd
            let writer =
                ZstdWriter::new(buffered_file, 1).context("Failed to create zstd encoder")?;
            Ok(Box::new(writer))
        } else {
            // Vanilla FASTQ
            Ok(Box::new(StandardWriter(buffered_file)))
        }
    }
}

// JSON log structure
#[derive(Serialize, Deserialize)]
pub struct FilterLog {
    version: String,
    index: String,
    input1: String,
    input2: Option<String>,
    output: String,
    k: usize,
    w: usize,
    m: usize,
    n: usize,
    invert: bool,
    rename: bool,
    seqs_in: u64,
    seqs_out: u64,
    seqs_removed: u64,
    seqs_removed_proportion: f64,
    bp_in: u64,
    bp_out: u64,
    bp_removed: u64,
    bp_removed_proportion: f64,
    time: f64,
    seqs_per_second: u64,
    bp_per_second: u64,
}

pub fn run<P: AsRef<Path>>(
    minimizers_path: P,
    input_path: &str,
    input2_path: Option<&str>,
    output_path: &str,
    min_matches: usize,
    prefix_length: usize,
    log_path: Option<&PathBuf>,
    invert: bool,
    rename: bool,
) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION").to_string();

    // Build preamble message
    let mut input_mode = String::new();
    let mut options = Vec::<String>::new();
    let paired_stdin = input_path == "-" && input2_path.is_some() && input2_path.unwrap() == "-";
    if paired_stdin {
        input_mode.push_str("interleaved pairs from stdin");
    } else if let Some(_) = input2_path {
        input_mode.push_str("pairs from separate files");
    } else {
        input_mode.push_str("single");
    }
    options.push(format!("min_matches={}", min_matches));
    if prefix_length > 0 {
        options.push(format!("prefix_length={}", prefix_length));
    }
    if invert {
        options.push("invert".to_string());
    }
    if rename {
        options.push("rename".to_string());
    }

    eprintln!(
        "Deacon v{}; mode: {}; options: {}",
        version,
        input_mode,
        options.join(", ")
    );

    // Load minimizers hashes and parse header
    let (minimizer_hashes, header) = load_minimizer_hashes(&minimizers_path)?;

    let kmer_length = header.kmer_length();
    let window_length = header.window_length();

    let load_time = start_time.elapsed();
    eprintln!(
        "Loaded index (k={}, w={}) in {:.2?}",
        kmer_length, window_length, load_time
    );

    // Create the appropriate writer based on the output path
    let mut writer = get_writer(output_path)?;

    // A progress bar would require a denominator, so let's spin
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[".  ", ".. ", "...", " ..", "  .", "   "])
            .template("{msg}{spinner} ")?,
    );
    spinner.set_message("Filtering");
    spinner.enable_steady_tick(Duration::from_millis(100));

    // Init counters
    let mut total_seqs = 0;
    let mut filtered_seqs = 0;
    let mut total_bp = 0;
    let mut output_bp = 0;
    let mut filtered_bp = 0;
    let mut output_seq_counter = 0;

    if paired_stdin {
        process_interleaved_paired_seqs(
            &minimizer_hashes,
            &mut writer,
            min_matches,
            prefix_length,
            kmer_length,
            window_length,
            invert,
            rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            start_time,
        )?;
    } else if let Some(input2_path) = input2_path {
        process_paired_seqs(
            &minimizer_hashes,
            input_path,
            input2_path,
            &mut writer,
            min_matches,
            prefix_length,
            kmer_length,
            window_length,
            invert,
            rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            start_time,
        )?;
    } else {
        process_single_seqs(
            &minimizer_hashes,
            input_path,
            &mut writer,
            min_matches,
            prefix_length,
            kmer_length,
            window_length,
            invert,
            rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            start_time,
        )?;
    }

    // Remainder of the function unchanged
    // Ensure everything is flushed
    writer.flush_all()?;

    let total_time = start_time.elapsed();
    let seqs_per_sec = total_seqs as f64 / total_time.as_secs_f64();
    let bp_per_sec = total_bp as f64 / total_time.as_secs_f64();
    // Convert to Mbp/s for display only
    let mbp_per_sec = bp_per_sec / 1_000_000.0;

    // Calculate filtered proportion directly
    let filtered_proportion = if total_seqs > 0 {
        filtered_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    // Calculate filtered base pair proportion
    let filtered_bp_proportion = if total_bp > 0 {
        filtered_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    // Finish spinner with final message
    spinner.finish_with_message(format!(
        "Removed {}/{} sequences ({:.3}%) , {}/{} bp ({:.3}%)",
        filtered_seqs,
        total_seqs,
        filtered_proportion * 100.0,
        filtered_bp,
        total_bp,
        filtered_bp_proportion * 100.0
    ));

    // Print completion message without spinner
    eprintln!(
        "Completed in {:.2?}. Speed: {:.0} seqs/s ({:.1} Mbp/s)",
        total_time, seqs_per_sec, mbp_per_sec
    );

    // Build and write a JSON log if path provided
    if let Some(log_file) = log_path {
        // Get number of sequences passing filter
        let seqs_out = total_seqs - filtered_seqs;

        let log = FilterLog {
            version: version,
            index: minimizers_path.as_ref().to_string_lossy().to_string(),
            input1: input_path.to_string(),
            input2: input2_path.map(|s| s.to_string()),
            output: output_path.to_string(),
            k: kmer_length,
            w: window_length,
            m: min_matches,
            n: prefix_length,
            invert,
            rename,
            seqs_in: total_seqs as u64,
            seqs_out: seqs_out as u64,
            seqs_removed: filtered_seqs as u64,
            seqs_removed_proportion: filtered_proportion,
            bp_in: total_bp as u64,
            bp_out: output_bp as u64,
            bp_removed: filtered_bp as u64,
            bp_removed_proportion: filtered_bp_proportion,
            time: total_time.as_secs_f64(),
            seqs_per_second: seqs_per_sec as u64,
            bp_per_second: bp_per_sec as u64,
        };

        // Write log file
        let file =
            File::create(log_file).context(format!("Failed to create log file: {:?}", log_file))?;
        let writer = BufWriter::new(file);

        // Serialize and write the log as JSON
        serde_json::to_writer_pretty(writer, &log).context("Failed to write JSON log")?;

        eprintln!("JSON log saved to: {:?}", log_file);
    }

    Ok(())
}

fn process_single_seqs(
    minimizer_hashes: &FxHashSet<u64>,
    input_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    min_matches: usize,
    prefix_length: usize,
    kmer_length: usize,
    window_length: usize,
    invert: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    start_time: Instant,
) -> Result<()> {
    // Opt: pre-allocate buffers for reuse
    let mut minimizer_buffer = Vec::with_capacity(64);
    let mut output_record_buffer = Vec::with_capacity(1024);
    let mut seen_hits = FxHashSet::default();

    // Parse FASTX from input (file or stdin)
    let mut reader = if input_path == "-" {
        parse_fastx_stdin()?
    } else {
        parse_fastx_file(input_path)?
    };

    while let Some(record) = reader.next() {
        let record = record?;

        let seq = record.seq();
        *total_seqs += 1;
        *total_bp += seq.len() as u64;

        // Check for minimizer hits
        let mut hit_count = 0;

        if seq.len() >= kmer_length {
            // Reuse pre-allocated buffers
            minimizer_buffer.clear();
            seen_hits.clear();

            // Apply prefix length limit if specified (if > 0)
            let effective_seq = if prefix_length > 0 && seq.len() > prefix_length {
                &seq.as_ref()[..prefix_length]
            } else {
                seq.as_ref()
            };

            // Get minimizer hash values using parameters from header
            fill_minimizer_hashes(
                effective_seq,
                kmer_length,
                window_length,
                &mut minimizer_buffer,
            );

            // Count distinct minimizer hits
            for &hash in &minimizer_buffer {
                if minimizer_hashes.contains(&hash) && seen_hits.insert(hash) {
                    hit_count += 1;
                    if hit_count >= min_matches {
                        break;
                    }
                }
            }
        }

        // Determine if we should output this sequence based on hit count and invert flag
        let should_output = if !invert {
            // When not inverted, keep sequences with fewer than min_matches
            hit_count < min_matches
        } else {
            // When inverted, keep sequences with greater than or equal to min_matches
            hit_count >= min_matches
        };

        if should_output {
            // Track output base pairs
            *output_bp += seq.len() as u64;

            // Format as FASTX directly to a byte buffer
            output_record_buffer.clear();

            // Increment output sequence counter
            *output_seq_counter += 1;

            // Write the record based on format, with optional renaming
            output_fastx_record(
                &record,
                &mut output_record_buffer,
                rename,
                *output_seq_counter,
            );

            // Write the whole record at once
            writer.write_all(&output_record_buffer)?;
        } else {
            *filtered_seqs += 1;
            *filtered_bp += seq.len() as u64; // Track filtered base pairs
        }

        // Periodically flush and update spinner message with detailed stats
        if *total_seqs % (100_000) == 0 {
            writer.flush()?;

            let elapsed = start_time.elapsed();
            let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
            let mbp_per_sec = bp_per_sec / 1_000_000.0;

            // Calculate filtered proportion directly
            let filtered_proportion = if *total_seqs > 0 {
                *filtered_seqs as f64 / *total_seqs as f64
            } else {
                0.0
            };

            // Calculate filtered base pair proportion
            let filtered_bp_proportion = if *total_bp > 0 {
                *filtered_bp as f64 / *total_bp as f64
            } else {
                0.0
            };

            // Update spinner message with detailed stats
            spinner.set_message(format!(
                "Removed {}/{} seqs ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
                filtered_seqs,
                total_seqs,
                filtered_proportion * 100.0, // Convert to percentage for display
                filtered_bp,
                total_bp,
                filtered_bp_proportion * 100.0,
                seqs_per_sec,
                mbp_per_sec
            ));
        }
    }

    Ok(())
}

fn process_paired_seqs(
    minimizer_hashes: &FxHashSet<u64>,
    input1_path: &str,
    input2_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    min_matches: usize,
    prefix_length: usize,
    kmer_length: usize,
    window_length: usize,
    invert: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    start_time: Instant,
) -> Result<()> {
    // Opt: pre-allocate buffers for reuse
    let mut minimizer_buffer1 = Vec::with_capacity(64);
    let mut minimizer_buffer2 = Vec::with_capacity(64);
    let mut output_record_buffer = Vec::with_capacity(1024);
    let mut seen_hits1 = FxHashSet::default();
    let mut seen_hits2 = FxHashSet::default();

    // Open both input files
    let mut reader1 = if input1_path == "-" {
        parse_fastx_stdin()?
    } else {
        parse_fastx_file(input1_path)?
    };

    let mut reader2 = parse_fastx_file(input2_path)?;

    // Process read pairs
    while let (Some(record1_res), Some(record2_res)) = (reader1.next(), reader2.next()) {
        let record1 = record1_res?;
        let record2 = record2_res?;

        let seq1 = record1.seq();
        let seq2 = record2.seq();

        // Count as two sequences
        *total_seqs += 2;
        *total_bp += (seq1.len() + seq2.len()) as u64;

        // Check for minimizer hits in read 1
        let mut hit_count1 = 0;
        if seq1.len() >= kmer_length {
            minimizer_buffer1.clear();
            seen_hits1.clear();

            // Apply prefix length limit if specified (if > 0)
            let effective_seq = if prefix_length > 0 && seq1.len() > prefix_length {
                &seq1.as_ref()[..prefix_length]
            } else {
                seq1.as_ref()
            };

            // Get minimizer hash values using parameters from header
            fill_minimizer_hashes(
                effective_seq,
                kmer_length,
                window_length,
                &mut minimizer_buffer1,
            );

            // Count distinct minimizer hits
            for &hash in &minimizer_buffer1 {
                if minimizer_hashes.contains(&hash) && seen_hits1.insert(hash) {
                    hit_count1 += 1;
                }
            }
        }

        // Check for minimizer hits in read 2
        let mut hit_count2 = 0;
        if seq2.len() >= kmer_length {
            minimizer_buffer2.clear();
            seen_hits2.clear();

            // Apply prefix length limit if specified (if > 0)
            let effective_seq = if prefix_length > 0 && seq2.len() > prefix_length {
                &seq2.as_ref()[..prefix_length]
            } else {
                seq2.as_ref()
            };

            // Get minimizer hash values using parameters from header
            fill_minimizer_hashes(
                effective_seq,
                kmer_length,
                window_length,
                &mut minimizer_buffer2,
            );

            // Count distinct minimizer hits
            for &hash in &minimizer_buffer2 {
                if minimizer_hashes.contains(&hash) && seen_hits2.insert(hash) {
                    hit_count2 += 1;
                }
            }
        }

        // Total hit count for the read pair
        let total_hit_count = hit_count1 + hit_count2;

        // Determine if we should output this read pair based on combined hit count and invert flag
        let should_output = if !invert {
            // When not inverted, keep pairs with fewer than min_matches combined hits
            total_hit_count < min_matches
        } else {
            // When inverted, keep pairs with greater than or equal to min_matches combined hits
            total_hit_count >= min_matches
        };

        if should_output {
            // Track output base pairs
            *output_bp += (seq1.len() + seq2.len()) as u64;

            // Increment output sequence counter (twice, once for each read)
            *output_seq_counter += 2;

            // Write interleaved seqs
            // Format s1 as FASTX to byte buffer
            output_record_buffer.clear();
            output_fastx_record(
                &record1,
                &mut output_record_buffer,
                rename,
                *output_seq_counter - 1, // Use the correct counter for s1
            );
            writer.write_all(&output_record_buffer)?;

            // Format s2 as FASTX to byte buffer
            output_record_buffer.clear();
            output_fastx_record(
                &record2,
                &mut output_record_buffer,
                rename,
                *output_seq_counter, // Use the correct counter for s2
            );
            writer.write_all(&output_record_buffer)?;
        } else {
            *filtered_seqs += 2; // Both seqs filtered out
            *filtered_bp += (seq1.len() + seq2.len()) as u64; // Track filtered base pairs
        }

        // Periodically flush and update spinner message with detailed stats
        if *total_seqs % (100_000) == 0 {
            writer.flush()?;

            let elapsed = start_time.elapsed();
            let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
            let mbp_per_sec = bp_per_sec / 1_000_000.0;

            // Calculate filtered proportion directly
            let filtered_proportion = if *total_seqs > 0 {
                *filtered_seqs as f64 / *total_seqs as f64
            } else {
                0.0
            };

            // Calculate filtered base pair proportion
            let filtered_bp_proportion = if *total_bp > 0 {
                *filtered_bp as f64 / *total_bp as f64
            } else {
                0.0
            };

            // Update spinner message with detailed stats
            spinner.set_message(format!(
                "Removed {}/{} seqs ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
                filtered_seqs,
                total_seqs,
                filtered_proportion * 100.0, // Convert to percentage
                filtered_bp,
                total_bp,
                filtered_bp_proportion * 100.0,
                seqs_per_sec,
                mbp_per_sec
            ));
        }
    }

    Ok(())
}

fn process_interleaved_paired_seqs(
    minimizer_hashes: &FxHashSet<u64>,
    writer: &mut Box<dyn FastxWriter>,
    min_matches: usize,
    prefix_length: usize,
    kmer_length: usize,
    window_length: usize,
    invert: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    start_time: Instant,
) -> Result<()> {
    // Pre-allocate buffers for reuse
    let mut minimizer_buffer1 = Vec::with_capacity(64);
    let mut minimizer_buffer2 = Vec::with_capacity(64);
    let mut output_record_buffer = Vec::with_capacity(1024);
    let mut seen_hits1 = FxHashSet::default();
    let mut seen_hits2 = FxHashSet::default();

    // Parse FASTX from stdin
    let mut reader = parse_fastx_stdin()?;
    let mut record_counter = 0;

    // Process records in pairs
    loop {
        // Read the first record and extract all needed data
        let (record1_id, record1_seq, record1_qual, record1_format) = match reader.next() {
            Some(result) => {
                record_counter += 1;
                let record = result?;
                // Extract all data we need from the record
                let id = record.id().to_vec();
                let seq = record.seq().to_vec();
                let qual = record.qual().map(|q| q.to_vec());
                let format = record.format();
                (id, seq, qual, format)
            }
            None => break, // End of input
        };

        // Read the second record and extract all needed data
        let (record2_id, record2_seq, record2_qual, record2_format) = match reader.next() {
            Some(result) => {
                record_counter += 1;
                let record = result?;
                let id = record.id().to_vec();
                let seq = record.seq().to_vec();
                let qual = record.qual().map(|q| q.to_vec());
                let format = record.format();
                (id, seq, qual, format)
            }
            None => {
                // Check if we have record1 but no record2 (mispaired)
                return Err(anyhow::anyhow!(
                    "Uneven number of interleaved sequence pairs. Found {} records.",
                    record_counter
                ));
            }
        };

        // Count as two sequences
        *total_seqs += 2;
        *total_bp += (record1_seq.len() + record2_seq.len()) as u64;

        // Check for minimizer hits in seq 1
        let mut hit_count1 = 0;
        if record1_seq.len() >= kmer_length {
            minimizer_buffer1.clear();
            seen_hits1.clear();

            // Apply prefix length limit if specified
            let effective_seq = if prefix_length > 0 && record1_seq.len() > prefix_length {
                &record1_seq[..prefix_length]
            } else {
                &record1_seq
            };

            // Get minimizer hash values
            fill_minimizer_hashes(
                effective_seq,
                kmer_length,
                window_length,
                &mut minimizer_buffer1,
            );

            // Count distinct minimizer hits
            for &hash in &minimizer_buffer1 {
                if minimizer_hashes.contains(&hash) && seen_hits1.insert(hash) {
                    hit_count1 += 1;
                }
            }
        }

        // Check for minimizer hits in seq 2
        let mut hit_count2 = 0;
        if record2_seq.len() >= kmer_length {
            minimizer_buffer2.clear();
            seen_hits2.clear();

            // Apply prefix length limit if specified
            let effective_seq = if prefix_length > 0 && record2_seq.len() > prefix_length {
                &record2_seq[..prefix_length]
            } else {
                &record2_seq
            };

            // Get minimizer hash values
            fill_minimizer_hashes(
                effective_seq,
                kmer_length,
                window_length,
                &mut minimizer_buffer2,
            );

            // Count distinct minimizer hits
            for &hash in &minimizer_buffer2 {
                if minimizer_hashes.contains(&hash) && seen_hits2.insert(hash) {
                    hit_count2 += 1;
                }
            }
        }

        // Total hit count for the seq pair
        let total_hit_count = hit_count1 + hit_count2;

        // Determine if we should output this sequence pair based on combined hit count and invert flag
        let should_output = if !invert {
            // When not inverted, keep pairs with fewer than min_matches combined hits
            total_hit_count < min_matches
        } else {
            // When inverted, keep pairs with greater than or equal to min_matches combined hits
            total_hit_count >= min_matches
        };

        if should_output {
            // Track output base pairs
            *output_bp += (record1_seq.len() + record2_seq.len()) as u64;

            // Increment output sequence counter (twice, once for each sequence)
            *output_seq_counter += 2;

            // Format and write record 1
            output_record_buffer.clear();
            output_fastx_record_from_parts(
                &record1_id,
                &record1_seq,
                record1_qual.as_deref(),
                record1_format,
                &mut output_record_buffer,
                rename,
                *output_seq_counter - 1,
            );
            writer.write_all(&output_record_buffer)?;

            // Format and write record 2
            output_record_buffer.clear();
            output_fastx_record_from_parts(
                &record2_id,
                &record2_seq,
                record2_qual.as_deref(),
                record2_format,
                &mut output_record_buffer,
                rename,
                *output_seq_counter,
            );
            writer.write_all(&output_record_buffer)?;
        } else {
            *filtered_seqs += 2; // Both seqs filtered out
            *filtered_bp += (record1_seq.len() + record2_seq.len()) as u64; // Track filtered base pairs
        }

        // Periodically flush and update spinner message with detailed stats
        if *total_seqs % (100_000) == 0 {
            writer.flush()?;

            let elapsed = start_time.elapsed();
            let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
            let mbp_per_sec = bp_per_sec / 1_000_000.0;

            // Calculate filtered proportion directly
            let filtered_proportion = if *total_seqs > 0 {
                *filtered_seqs as f64 / *total_seqs as f64
            } else {
                0.0
            };

            // Calculate filtered base pair proportion
            let filtered_bp_proportion = if *total_bp > 0 {
                *filtered_bp as f64 / *total_bp as f64
            } else {
                0.0
            };

            // Update spinner message with detailed stats
            spinner.set_message(format!(
                "Removed {}/{} seqs ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
                filtered_seqs,
                total_seqs,
                filtered_proportion * 100.0, // Convert to percentage for display
                filtered_bp,
                total_bp,
                filtered_bp_proportion * 100.0,
                seqs_per_sec,
                mbp_per_sec
            ));
        }
    }

    Ok(())
}

/// Push FASTA or FASTQ record to output buffer
fn output_fastx_record(
    record: &SequenceRecord,
    buffer: &mut Vec<u8>,
    rename: bool,
    seq_number: u64,
) {
    match record.format() {
        Format::Fasta => {
            buffer.push(b'>');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(record.id());
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(record.seq().as_ref());
            buffer.push(b'\n');
        }
        Format::Fastq => {
            buffer.push(b'@');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(record.id());
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(record.seq().as_ref());
            buffer.extend_from_slice(b"\n+\n");
            if let Some(qual) = record.qual() {
                buffer.extend_from_slice(qual);
            }
            buffer.push(b'\n');
        }
    }
}

/// Push FASTA or FASTQ record to output buffer from component parts
/// Workaround for borrowing misery with interleaved pairs from stdin
fn output_fastx_record_from_parts(
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
    format: Format,
    buffer: &mut Vec<u8>,
    rename: bool,
    seq_number: u64,
) {
    match format {
        Format::Fasta => {
            buffer.push(b'>');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(id);
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(seq);
            buffer.push(b'\n');
        }
        Format::Fastq => {
            buffer.push(b'@');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(id);
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(seq);
            buffer.extend_from_slice(b"\n+\n");
            if let Some(qual_data) = qual {
                buffer.extend_from_slice(qual_data);
            }
            buffer.push(b'\n');
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::write_minimizers;
    use crate::index_format::IndexHeader;
    use tempfile::TempDir;

    #[allow(dead_code)] // Suppress unused warnings
    fn create_test_index() -> (PathBuf, IndexHeader, TempDir) {
        // Create a temporary directory
        let temp_dir = TempDir::new().unwrap();
        let index_path = temp_dir.path().join("test.idx");

        // Create dummy minimizers
        let minimizers: FxHashSet<u64> = [1, 2, 3, 4, 5].iter().cloned().collect();
        let header = IndexHeader::new(5, 3);

        // Write to file
        write_minimizers(&minimizers, &header, Some(&index_path)).unwrap();

        // Return the TempDir along with the other values to keep it in scope
        (index_path, header, temp_dir)
    }

    #[test]
    fn test_filter_log() {
        // Create a sample log
        let log = FilterLog {
            version: "0.1.0".to_string(),
            index: "test.idx".to_string(),
            input1: "test.fastq".to_string(),
            input2: Some("test2.fastq".to_string()), // Added for paired-end support
            output: "output.fastq".to_string(),
            k: 31,
            w: 21,
            m: 1,
            n: 0,
            invert: false,
            rename: false,
            seqs_in: 100,
            seqs_out: 90,
            seqs_removed: 10,
            seqs_removed_proportion: 0.1,
            bp_in: 10000,
            bp_out: 9000,
            bp_removed: 1000,
            bp_removed_proportion: 0.1,
            time: 1.5,
            seqs_per_second: 66,
            bp_per_second: 6666,
        };

        // Test JSON ser+de
        let json = serde_json::to_string(&log).unwrap();
        let parsed: FilterLog = serde_json::from_str(&json).unwrap();

        // Check values
        assert_eq!(parsed.version, "0.1.0");
        assert_eq!(parsed.seqs_in, 100);
        assert_eq!(parsed.seqs_removed_proportion, 0.1);
        assert_eq!(parsed.input1, "test.fastq");
        assert_eq!(parsed.input2, Some("test2.fastq".to_string()));
        assert_eq!(parsed.output, "output.fastq");
    }
}
