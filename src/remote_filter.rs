use crate::FilterConfig;
use crate::FilterSummary;
use crate::filter_common::get_minimizer_hashes_and_positions;
use crate::filter_common::get_paired_minimizer_hashes_and_positions;
use crate::filter_common::pair_matches;
use crate::filter_common::sequence_matches;
use crate::filter_common::{get_summary_index, meets_filtering_criteria};
use crate::index::load_minimizer_hashes;
#[cfg(feature = "server")]
use crate::server_common::{FilterResponse, PairedFilterRequest, UnpairedFilterRequest};
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use liblzma::write::XzEncoder;
use needletail::parse_fastx_file;
use needletail::parse_fastx_stdin;
use needletail::parser::Format;
use rayon::prelude::*;
#[cfg(feature = "server")]
use reqwest::blocking::Client;
use rustc_hash::FxHashSet;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::time::Instant;
use zstd::stream::write::Encoder as ZstdEncoder;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer

/// Data structure to hold a fastq record
struct RecordData {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
    format: Format,
}
trait FastxWriter: Write {
    fn flush_all(&mut self) -> io::Result<()>;
}

trait CompressionEncoder: Write {
    fn finish(self: Box<Self>) -> io::Result<()>;
}

#[derive(Debug, Clone, Copy)]
enum CompressionFormat {
    None,
    Gzip,
    Zstd,
    Xz,
}

impl CompressionFormat {
    fn from_extension(path: &str) -> Self {
        if path.ends_with(".gz") {
            Self::Gzip
        } else if path.ends_with(".zst") {
            Self::Zstd
        } else if path.ends_with(".xz") {
            Self::Xz
        } else {
            Self::None
        }
    }

    fn validate_compression_level(&self, level: u8) -> Result<()> {
        match self {
            Self::None => Ok(()),
            Self::Gzip => {
                if !(1..=9).contains(&level) {
                    Err(anyhow::anyhow!(
                        "Invalid gzip compression level {}. Must be between 1 and 9.",
                        level
                    ))
                } else {
                    Ok(())
                }
            }
            Self::Zstd => {
                if !(1..=22).contains(&level) {
                    Err(anyhow::anyhow!(
                        "Invalid zstd compression level {}. Must be between 1 and 22.",
                        level
                    ))
                } else {
                    Ok(())
                }
            }
            Self::Xz => {
                if level > 9 {
                    Err(anyhow::anyhow!(
                        "Invalid xz compression level {}. Must be between 0 and 9.",
                        level
                    ))
                } else {
                    Ok(())
                }
            }
        }
    }
}

impl<W: Write> CompressionEncoder for GzEncoder<W> {
    fn finish(mut self: Box<Self>) -> io::Result<()> {
        self.try_finish()
    }
}

impl<W: Write> CompressionEncoder for ZstdEncoder<'static, W> {
    fn finish(self: Box<Self>) -> io::Result<()> {
        (*self).finish().map(|_| ())
    }
}

impl<W: Write> CompressionEncoder for XzEncoder<W> {
    fn finish(self: Box<Self>) -> io::Result<()> {
        (*self).finish().map(|_| ())
    }
}

struct CompressedWriter {
    encoder: Option<Box<dyn CompressionEncoder>>,
}

impl CompressedWriter {
    fn new(encoder: Box<dyn CompressionEncoder>) -> Self {
        Self {
            encoder: Some(encoder),
        }
    }

    fn uncompressed<W: Write>(writer: W) -> StandardWriter<W> {
        StandardWriter(writer)
    }
}

impl Write for CompressedWriter {
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

impl FastxWriter for CompressedWriter {
    fn flush_all(&mut self) -> io::Result<()> {
        if let Some(encoder) = self.encoder.take() {
            encoder.finish()?;
        }
        Ok(())
    }
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

// Return a file writer appropriate for the output path extension
fn get_writer(output_path: &str, compression_level: u8) -> Result<Box<dyn FastxWriter>> {
    if output_path == "-" {
        // Write to stdout
        let stdout = io::stdout();
        let writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, stdout);
        Ok(Box::new(CompressedWriter::uncompressed(writer)))
    } else {
        // Write to file with extension-appropriate encoder
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(output_path)
            .context(format!("Failed to create output file: {output_path}"))?;

        let buffered_file = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);
        let format = CompressionFormat::from_extension(output_path);

        // Validate compression level for the format
        format.validate_compression_level(compression_level)?;

        match format {
            CompressionFormat::None => Ok(Box::new(CompressedWriter::uncompressed(buffered_file))),
            CompressionFormat::Gzip => {
                let encoder =
                    GzEncoder::new(buffered_file, Compression::new(compression_level as u32));
                Ok(Box::new(CompressedWriter::new(Box::new(encoder))))
            }
            CompressionFormat::Zstd => {
                let encoder = ZstdEncoder::new(buffered_file, compression_level as i32)
                    .context("Failed to create zstd encoder")?;
                Ok(Box::new(CompressedWriter::new(Box::new(encoder))))
            }
            CompressionFormat::Xz => {
                let encoder = XzEncoder::new(buffered_file, compression_level as u32);
                Ok(Box::new(CompressedWriter::new(Box::new(encoder))))
            }
        }
    }
}

pub fn unpaired_should_keep(
    input_minimizers_and_positions: &Vec<(Vec<u64>, Vec<u32>, Vec<u8>)>,
    kmer_length: u8,
    index_minimizers: &FxHashSet<u64>,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
    debug: bool,
) -> Vec<(bool, usize, usize, Vec<String>)> {
    input_minimizers_and_positions
        .par_iter()
        .map(|(minimizers, positions, seq)| {
            let (hit_count, hit_kmers) = sequence_matches(
                index_minimizers,
                minimizers,
                positions,
                seq,
                kmer_length,
                debug,
            );
            (
                meets_filtering_criteria(
                    hit_count,
                    minimizers.len(),
                    abs_threshold,
                    rel_threshold,
                    deplete,
                ),
                hit_count,
                minimizers.len(),
                hit_kmers,
            )
        })
        .collect()
}

pub fn paired_should_keep(
    input_minimizers_and_positions: &Vec<(Vec<u64>, Vec<u32>, Vec<&[u8]>)>,
    kmer_length: u8,
    index_minimizers: &FxHashSet<u64>,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
    debug: bool,
) -> Vec<(bool, usize, usize, Vec<String>)> {
    input_minimizers_and_positions
        .par_iter()
        .map(|(minimizers, positions, seq)| {
            let (pair_hit_count, hit_kmers) = pair_matches(
                minimizers,
                positions,
                seq,
                index_minimizers,
                kmer_length,
                debug,
            );

            (
                meets_filtering_criteria(
                    pair_hit_count,
                    minimizers.len(),
                    abs_threshold,
                    rel_threshold,
                    deplete,
                ),
                pair_hit_count,
                minimizers.len(),
                hit_kmers,
            )
        })
        .collect()
}

/// Given a set of input minimizers from unpaired reads, check if they should be output
/// If index minimizers are provided, check locally.
/// If not, send to server for checking. Requires the `server` feature to be enabled.
pub fn check_single_inputs_should_be_output(
    index_minimizers: &Option<FxHashSet<u64>>,
    input_minimizers_and_positions: &Vec<(Vec<u64>, Vec<u32>, Vec<u8>)>,
    _server_address: &Option<String>,
    deplete: bool,
    kmer_length: u8,
    debug: bool,
    abs_threshold: usize,
    rel_threshold: f64,
) -> Vec<(bool, usize, usize, Vec<String>)> {
    // If index minimizers are provided, check if input matches locally
    if let Some(index_minimizers) = index_minimizers {
        unpaired_should_keep(
            input_minimizers_and_positions,
            kmer_length,
            index_minimizers,
            abs_threshold,
            rel_threshold,
            deplete,
            debug,
        )
    } else {
        // Else, send the input minimizers to the server for checking
        #[cfg(feature = "server")]
        {
            if _server_address.is_none() {
                panic!("Server address is required when using the server feature.");
            }
            let server_address = _server_address.as_ref().map(String::as_str).unwrap();
            // Create a client to send the minimizers to the server
            let client = Client::new();

            // Send the minimizers as a POST request
            let response = client
                .post(server_address.to_owned() + "/should_output_unpaired")
                .json(&UnpairedFilterRequest {
                    input: input_minimizers_and_positions.to_vec(),
                    abs_threshold,
                    rel_threshold,
                    deplete,
                    kmer_length,
                    debug,
                })
                .send()
                .unwrap();

            // Check if the response indicates a match
            if response.status().is_success() {
                response.json::<FilterResponse>().unwrap().should_output
            } else {
                panic!("Server returned an error: {}", response.status())
            }
        }
        #[cfg(not(feature = "server"))]
        {
            panic!("Server feature is not enabled. Cannot check input against index.");
        }
    }
}

/// Given a set of input minimizers from paired reads, check if they should be output
/// If index minimizers are provided, check locally.
/// If not, send to server for checking. Requires the `server` feature to be enabled.
pub fn check_paired_inputs_should_be_output(
    index_minimizers: &Option<FxHashSet<u64>>,
    input_minimizers_and_positions: &Vec<(Vec<u64>, Vec<u32>, Vec<&[u8]>)>,
    _server_address: &Option<String>,
    deplete: bool,
    kmer_length: u8,
    debug: bool,
    abs_threshold: usize,
    rel_threshold: f64,
) -> Vec<(bool, usize, usize, Vec<String>)> {
    // If index minimizers are provided, check if input matches locally
    if let Some(index_minimizers) = index_minimizers {
        paired_should_keep(
            input_minimizers_and_positions,
            kmer_length,
            index_minimizers,
            abs_threshold,
            rel_threshold,
            deplete,
            debug,
        )
    } else {
        // Else, send the input minimizers to the server for checking
        #[cfg(feature = "server")]
        {
            if _server_address.is_none() {
                panic!("Server address is required when using the server feature.");
            }
            let server_address = _server_address.as_ref().map(String::as_str).unwrap();

            // Quickly wrangle the seqs into vecs instead of slices so serde can cope
            // Not perfect, but if it has to happen anywhere, here is the best
            let input_minimizers_and_positions: Vec<(Vec<u64>, Vec<u32>, Vec<Vec<u8>>)> =
                input_minimizers_and_positions
                    .iter()
                    .map(|(minimizers, positions, seqs)| {
                        (
                            minimizers.to_vec(),
                            positions.to_vec(),
                            seqs.iter().map(|s| s.to_vec()).collect(),
                        )
                    })
                    .collect();

            // Create a client to send the minimizers to the server
            let client = Client::new();

            // Send the minimizers as a POST request
            let response = client
                .post(server_address.to_owned() + "/should_output_paired")
                .json(&PairedFilterRequest {
                    input: input_minimizers_and_positions.to_vec(),
                    abs_threshold,
                    rel_threshold,
                    deplete,
                    kmer_length,
                    debug,
                })
                .send()
                .unwrap();

            // Check if the response indicates a match
            if response.status().is_success() {
                response.json::<FilterResponse>().unwrap().should_output
            } else {
                panic!("Server returned an error: {}", response.status())
            }
        }
        #[cfg(not(feature = "server"))]
        {
            panic!("Server feature is not enabled. Cannot check input against index.");
        }
    }
}

/// Run deacon filter with the provided parameters.
pub fn run(config: &FilterConfig) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    let tool_version = format!("deacon {version}");

    // Configure thread pool if nonzero
    if config.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build_global()
            .context("Failed to initialize thread pool")?;
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

    eprintln!(
        "Deacon v{}; mode: {}; input: {}; options: {}",
        version,
        mode,
        input_type,
        options.join(", ")
    );

    // Load minimizers hashes and parse header
    let (minimizer_hashes, header) =
        load_minimizer_hashes(&config.minimizers_path, &config.server_address)?;

    let kmer_length = header.kmer_length();
    let window_size = header.window_size();

    let load_time = start_time.elapsed();
    eprintln!("Loaded index (k={kmer_length}, w={window_size}) in {load_time:.2?}");

    // Create the appropriate writer(s) based on the output path(s)
    let mut writer = get_writer(config.output_path, config.compression_level)?;
    let mut writer2 = if let (Some(output2), Some(_)) = (config.output2_path, config.input2_path) {
        // Only create second writer if both output2 and input2 are specified
        Some(get_writer(output2, config.compression_level)?)
    } else {
        None
    };

    // A progress bar would require a denominator, so let's spin
    let spinner = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[".  ", ".. ", "...", " ..", "  .", "   "])
            .template("{msg}{spinner} ")?,
    );
    spinner.set_message("Filtering");

    // Init counters
    let mut total_seqs = 0;
    let mut filtered_seqs = 0;
    let mut total_bp = 0;
    let mut output_bp = 0;
    let mut filtered_bp = 0;
    let mut output_seq_counter = 0;

    // Start timer for filtering rate calculation (excludes index loading time)
    let filtering_start_time = Instant::now();

    if paired_stdin {
        process_interleaved_paired_seqs(
            &minimizer_hashes,
            &mut writer,
            writer2.as_mut(),
            config.abs_threshold,
            config.rel_threshold,
            config.prefix_length as u8,
            kmer_length,
            window_size,
            config.deplete,
            config.rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            filtering_start_time,
            &config.server_address,
            config.debug,
        )?;
    } else if let Some(input2_path) = config.input2_path {
        process_paired_seqs(
            &minimizer_hashes,
            config.input_path,
            input2_path,
            &mut writer,
            writer2.as_mut(),
            config.abs_threshold,
            config.rel_threshold,
            config.prefix_length as u8,
            kmer_length,
            window_size,
            config.deplete,
            config.rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            filtering_start_time,
            &config.server_address,
            config.debug,
        )?;
    } else {
        process_single_seqs(
            &minimizer_hashes,
            config.input_path,
            &mut writer,
            config.abs_threshold,
            config.rel_threshold,
            config.prefix_length as u8,
            kmer_length,
            window_size,
            config.deplete,
            config.rename,
            &mut total_seqs,
            &mut filtered_seqs,
            &mut total_bp,
            &mut output_bp,
            &mut filtered_bp,
            &mut output_seq_counter,
            &spinner,
            filtering_start_time,
            &config.server_address,
            config.debug,
        )?;
    }

    writer.flush_all()?;
    if let Some(ref mut w2) = writer2 {
        w2.flush_all()?;
    }

    let total_time = start_time.elapsed();
    let seqs_per_sec = total_seqs as f64 / total_time.as_secs_f64();
    let bp_per_sec = total_bp as f64 / total_time.as_secs_f64();
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

    // Calculate output proportions
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

    // Finish and clear spinner, print final message
    spinner.finish_and_clear();
    eprintln!(
        "Retained {}/{} sequences ({:.3}%), {}/{} bp ({:.3}%)",
        output_seqs,
        total_seqs,
        output_seq_proportion * 100.0,
        output_bp,
        total_bp,
        output_bp_proportion * 100.0
    );

    // Print completion message with speed
    eprintln!(
        "Completed in {total_time:.2?}. Speed: {seqs_per_sec:.0} seqs/s ({mbp_per_sec:.1} Mbp/s)"
    );

    // Build and write a JSON summary if path provided
    if let Some(summary_file) = config.summary_path {
        // Get number of sequences passing filter
        let seqs_out = total_seqs - filtered_seqs;

        let summary = FilterSummary {
            version: tool_version,
            index: get_summary_index(&config.minimizers_path, &config.server_address),
            input: config.input_path.to_string(),
            input2: config.input2_path.map(|s| s.to_string()),
            output: config.output_path.to_string(),
            output2: config.output2_path.map(|s| s.to_string()),
            k: kmer_length,
            w: window_size,
            abs_threshold: config.abs_threshold,
            rel_threshold: config.rel_threshold,
            prefix_length: config.prefix_length,
            deplete: config.deplete,
            rename: config.rename,
            seqs_in: total_seqs,
            seqs_out,
            seqs_out_proportion: output_seq_proportion,
            seqs_removed: filtered_seqs,
            seqs_removed_proportion: filtered_proportion,
            bp_in: total_bp,
            bp_out: output_bp,
            bp_out_proportion: output_bp_proportion,
            bp_removed: filtered_bp,
            bp_removed_proportion: filtered_bp_proportion,
            time: total_time.as_secs_f64(),
            seqs_per_second: seqs_per_sec as u64,
            bp_per_second: bp_per_sec as u64,
        };

        // Write summary file
        let file = File::create(summary_file)
            .context(format!("Failed to create summary: {summary_file:?}"))?;
        let writer = BufWriter::new(file);

        // Serialise and write the summary JSON
        serde_json::to_writer_pretty(writer, &summary).context("Failed to write summary")?;

        eprintln!("Summary saved to {summary_file:?}");
    }

    Ok(())
}

/// Filter a single (unpaired) sequence.
#[allow(clippy::too_many_arguments)]
fn process_single_seqs(
    minimizer_hashes: &Option<FxHashSet<u64>>,
    input_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: u8,
    kmer_length: u8,
    window_size: u8,
    deplete: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    filtering_start_time: Instant,
    server_address: &Option<String>,
    debug: bool,
) -> Result<()> {
    // Create a reader based on the input source
    let mut reader = if input_path == "-" {
        parse_fastx_stdin()?
    } else {
        parse_fastx_file(input_path)?
    };

    // Process in batches
    let batch_size = 10000;
    let mut output_record_buffer = Vec::with_capacity(1024);

    // Process batches
    loop {
        // Collect a batch of records with owned data
        let mut batch: Vec<RecordData> = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        // Fill the batch (sequential read from either stdin or file)
        for _ in 0..batch_size {
            if let Some(record_result) = reader.next() {
                match record_result {
                    Ok(record) => {
                        let record_data = RecordData {
                            id: record.id().to_vec(),
                            seq: record.seq().to_vec(),
                            qual: record.qual().map(|q| q.to_vec()),
                            format: record.format(),
                        };
                        batch.push(record_data);
                    }
                    Err(e) => return Err(e.into()),
                }
            } else {
                reached_end = true;
                break;
            }
        }

        if batch.is_empty() {
            break;
        }

        // Get batch minimizers in parallel
        let batch_result: Vec<(Vec<u64>, Vec<u32>, Vec<u8>)> = batch
            .par_iter()
            .map(|record_data| {
                let (hashes, positions, seqs) = get_minimizer_hashes_and_positions(
                    &record_data.seq,
                    prefix_length as usize,
                    kmer_length,
                    window_size,
                );
                (hashes, positions, seqs.to_vec())
                // get_hashes_from_record(record_data, kmer_length, prefix_length, window_size)
            })
            .collect();

        // let (batch_minimizers, batch_positons, effective_seqs): (Vec<Vec<u64>>, Vec<Vec<u32>>, Vec<&[u8]>) =
        //     batch_result.into_iter().multiunzip();

        // Check if minimizers match the index
        // Separated from initial par_iter to allow flexibility with local/server processing
        let batch_should_outputs = check_single_inputs_should_be_output(
            minimizer_hashes,
            &batch_result,
            server_address,
            deplete,
            kmer_length,
            debug,
            abs_threshold,
            rel_threshold,
        );

        // Process results sequentially to maintain order
        for (i, (should_output, hit_count, total_minimizers, hit_kmers)) in
            batch_should_outputs.into_iter().enumerate()
        {
            let record_data = &batch[i];
            let seq_len = record_data.seq.len();
            *total_seqs += 1;
            *total_bp += seq_len as u64;

            if debug {
                eprintln!(
                    "DEBUG: {} hits={}/{} keep={} kmers=[{}]",
                    String::from_utf8_lossy(&record_data.id),
                    hit_count,
                    total_minimizers,
                    should_output,
                    hit_kmers.join(",")
                );
            }

            if should_output {
                // Track output base pairs
                *output_bp += seq_len as u64;

                // Increment output sequence counter
                *output_seq_counter += 1;

                // Format as FASTX and write
                output_record_buffer.clear();
                output_fastx_record_from_parts(
                    &record_data.id,
                    &record_data.seq,
                    record_data.qual.as_deref(),
                    record_data.format,
                    &mut output_record_buffer,
                    rename,
                    *output_seq_counter,
                );
                writer.write_all(&output_record_buffer)?;
            } else {
                *filtered_seqs += 1;
                *filtered_bp += seq_len as u64;
            }
        }

        // Update spinner and flush periodically
        let elapsed = filtering_start_time.elapsed();
        let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
        let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
        let mbp_per_sec = bp_per_sec / 1_000_000.0;

        // Calculate output proportion
        let output_seqs = *total_seqs - *filtered_seqs;
        let output_proportion = if *total_seqs > 0 {
            output_seqs as f64 / *total_seqs as f64
        } else {
            0.0
        };

        // Calculate output base pair proportion
        let output_bp_proportion = if *total_bp > 0 {
            *output_bp as f64 / *total_bp as f64
        } else {
            0.0
        };

        // Update spinner message
        spinner.set_message(format!(
            "Retained {}/{} sequences ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
            output_seqs,
            total_seqs,
            output_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0,
            seqs_per_sec,
            mbp_per_sec
        ));

        // Flush writer periodically
        writer.flush()?;

        // Check if we've reached the end of the file/stdin
        if reached_end {
            break;
        }
    }

    Ok(())
}

/// Filter a pair of sequences
#[allow(clippy::too_many_arguments)]
fn process_paired_seqs(
    minimizer_hashes: &Option<FxHashSet<u64>>,
    input1_path: &str,
    input2_path: &str,
    writer: &mut Box<dyn FastxWriter>,
    mut writer2: Option<&mut Box<dyn FastxWriter>>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: u8,
    kmer_length: u8,
    window_size: u8,
    deplete: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    filtering_start_time: Instant,
    server_address: &Option<String>,
    debug: bool,
) -> Result<()> {
    // Open both input files
    let mut reader1 = if input1_path == "-" {
        parse_fastx_stdin()?
    } else {
        parse_fastx_file(input1_path)?
    };

    let mut reader2 = parse_fastx_file(input2_path)?;

    // Process in batches
    let batch_size = 10000;
    let mut output_record_buffer = Vec::with_capacity(1024);

    // Process batches
    loop {
        // Collect a batch of read pairs with owned data
        let mut batch1: Vec<RecordData> = Vec::with_capacity(batch_size);
        let mut batch2: Vec<RecordData> = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        // Fill the batch (sequential read from files)
        for _ in 0..batch_size {
            if let (Some(record1_res), Some(record2_res)) = (reader1.next(), reader2.next()) {
                match (record1_res, record2_res) {
                    (Ok(record1), Ok(record2)) => {
                        let record_data1 = RecordData {
                            id: record1.id().to_vec(),
                            seq: record1.seq().to_vec(),
                            qual: record1.qual().map(|q| q.to_vec()),
                            format: record1.format(),
                        };
                        let record_data2 = RecordData {
                            id: record2.id().to_vec(),
                            seq: record2.seq().to_vec(),
                            qual: record2.qual().map(|q| q.to_vec()),
                            format: record2.format(),
                        };
                        batch1.push(record_data1);
                        batch2.push(record_data2);
                    }
                    (Err(e), _) => return Err(e.into()),
                    (_, Err(e)) => return Err(e.into()),
                }
            } else {
                reached_end = true;
                break;
            }
        }

        if batch1.is_empty() {
            break;
        }

        // Get batch minimizers in parallel
        let batch_result: Vec<(Vec<u64>, Vec<u32>, Vec<&[u8]>)> = batch1
            .par_iter()
            .zip(batch2.par_iter())
            .map(|(record_data1, record_data2)| {
                get_paired_minimizer_hashes_and_positions(
                    &record_data1.seq,
                    &record_data2.seq,
                    prefix_length.into(),
                    kmer_length,
                    window_size,
                )
            })
            .collect();

        // let batch_result: Vec<(Vec<u64>, u8, u8)> = batch1
        //     .par_iter()
        //     .zip(batch2.par_iter())
        //     .map(|(record_data1, record_data2)| {
        //         get_hashes_from_record_pair(
        //             record_data1,
        //             record_data2,
        //             kmer_length,
        //             prefix_length,
        //             window_size,
        //         )
        //     })
        //     .collect();

        // let (batch_minimizers, seq_lens1, seq_lens2): (Vec<Vec<u64>>, Vec<u8>, Vec<u8>) =
        //     batch_result.into_iter().multiunzip();

        // Check if minimizers match the index
        // Separated from initial par_iter to allow flexibility with local/server processing
        // let batch_should_outputs = check_inputs_should_be_output(
        //     minimizer_hashes,
        //     &batch_minimizers,
        //     abs_threshold,
        //     rel_threshold,
        //     server_address,
        //     deplete,
        // );

        let batch_should_outputs = check_paired_inputs_should_be_output(
            minimizer_hashes,
            &batch_result,
            server_address,
            deplete,
            kmer_length,
            debug,
            abs_threshold,
            rel_threshold,
        );

        // Process results sequentially to maintain order
        for (i, (should_output, hit_count, total_minimizers, hit_kmers)) in
            batch_should_outputs.into_iter().enumerate()
        {
            let record_data1 = &batch1[i];
            let record_data2 = &batch2[i];
            let seq1_len = record_data1.seq.len();
            let seq2_len = record_data2.seq.len();

            *total_seqs += 2;
            *total_bp += (seq1_len + seq2_len) as u64;

            if debug && hit_count > 0 {
                eprintln!(
                    "DEBUG: {}/{} hits={}/{} keep={} kmers=[{}]",
                    String::from_utf8_lossy(&record_data1.id),
                    String::from_utf8_lossy(&record_data2.id),
                    hit_count,
                    total_minimizers,
                    should_output,
                    hit_kmers.join(",")
                );
            }

            if should_output {
                // Track output base pairs
                *output_bp += (seq1_len + seq2_len) as u64;

                // Increment output sequence counter (twice, once for each read)
                *output_seq_counter += 2;

                // Format s1 as FASTX to byte buffer and write to appropriate writer
                output_record_buffer.clear();
                output_fastx_record_from_parts(
                    &record_data1.id,
                    &record_data1.seq,
                    record_data1.qual.as_deref(),
                    record_data1.format,
                    &mut output_record_buffer,
                    rename,
                    *output_seq_counter - 1,
                );

                if let Some(ref mut w2) = writer2 {
                    // Write read 1 to primary writer
                    writer.write_all(&output_record_buffer)?;

                    // Format s2 as FASTX to byte buffer and write to second writer
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        &record_data2.id,
                        &record_data2.seq,
                        record_data2.qual.as_deref(),
                        record_data2.format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    w2.write_all(&output_record_buffer)?;
                } else {
                    // Interleaved output
                    writer.write_all(&output_record_buffer)?;

                    // Format s2 as FASTX to byte buffer
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        &record_data2.id,
                        &record_data2.seq,
                        record_data2.qual.as_deref(),
                        record_data2.format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    writer.write_all(&output_record_buffer)?;
                }
            } else {
                *filtered_seqs += 2; // Both seqs filtered out
                *filtered_bp += (seq1_len + seq2_len) as u64; // Track filtered base pairs
            }
        }

        // Update spinner and flush periodically
        let elapsed = filtering_start_time.elapsed();
        let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
        let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
        let mbp_per_sec = bp_per_sec / 1_000_000.0;

        // Calculate output proportion directly
        let output_seqs = *total_seqs - *filtered_seqs;
        let output_proportion = if *total_seqs > 0 {
            output_seqs as f64 / *total_seqs as f64
        } else {
            0.0
        };

        // Calculate output base pair proportion
        let output_bp_proportion = if *total_bp > 0 {
            *output_bp as f64 / *total_bp as f64
        } else {
            0.0
        };

        // Update spinner message with detailed stats
        spinner.set_message(format!(
            "Retained {}/{} sequences ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
            output_seqs,
            total_seqs,
            output_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0,
            seqs_per_sec,
            mbp_per_sec
        ));

        // Flush writer periodically
        writer.flush()?;
        if let Some(ref mut w2) = writer2 {
            w2.flush()?;
        }

        // Check if we've reached the end of the files
        if reached_end {
            break;
        }
    }

    Ok(())
}

/// Filter a pair of interleaved sequences
/// Functionally very similar to `process_paired_seqs`, but handles interleaved input
#[allow(clippy::too_many_arguments)]
fn process_interleaved_paired_seqs(
    minimizer_hashes: &Option<FxHashSet<u64>>,
    writer: &mut Box<dyn FastxWriter>,
    mut writer2: Option<&mut Box<dyn FastxWriter>>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: u8,
    kmer_length: u8,
    window_size: u8,
    deplete: bool,
    rename: bool,
    total_seqs: &mut u64,
    filtered_seqs: &mut u64,
    total_bp: &mut u64,
    output_bp: &mut u64,
    filtered_bp: &mut u64,
    output_seq_counter: &mut u64,
    spinner: &ProgressBar,
    filtering_start_time: Instant,
    server_address: &Option<String>,
    debug: bool,
) -> Result<()> {
    // Parse FASTX from stdin
    let mut reader = parse_fastx_stdin()?;
    let mut output_record_buffer = Vec::with_capacity(1024);
    let mut record_counter = 0;

    // Process in batches
    let batch_size = 10000;

    loop {
        // Collect a batch of read pairs with owned data
        let mut batch_pairs = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        // Fill the batch with interleaved pairs
        for _ in 0..batch_size {
            // Read the first record of the pair
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
                None => {
                    reached_end = true;
                    break; // End of input
                }
            };

            // Read the second record of the pair
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

            // Store the pair in the batch
            batch_pairs.push((
                RecordData {
                    id: record1_id,
                    seq: record1_seq,
                    qual: record1_qual,
                    format: record1_format,
                },
                RecordData {
                    id: record2_id,
                    seq: record2_seq,
                    qual: record2_qual,
                    format: record2_format,
                },
            ));
        }

        if batch_pairs.is_empty() {
            break;
        }

        // Get batch minimizers in parallel
        let batch_result: Vec<(Vec<u64>, Vec<u32>, Vec<&[u8]>)> = batch_pairs
            .par_iter()
            .map(|(record_data1, record_data2)| {
                get_paired_minimizer_hashes_and_positions(
                    &record_data1.seq,
                    &record_data2.seq,
                    prefix_length.into(),
                    kmer_length,
                    window_size,
                )
            })
            .collect();

        let batch_should_outputs = check_paired_inputs_should_be_output(
            minimizer_hashes,
            &batch_result,
            server_address,
            deplete,
            kmer_length,
            debug,
            abs_threshold,
            rel_threshold,
        );
        // let batch_result: Vec<(Vec<u64>, u8, u8)> = batch_pairs
        //     .par_iter()
        //     .map(|(record_data1, record_data2)| {
        //         get_hashes_from_record_pair(
        //             record_data1,
        //             record_data2,
        //             kmer_length,
        //             prefix_length,
        //             window_size,
        //         )
        //     })
        //     .collect();

        // let (batch_minimizers, seq_lens1, seq_lens2): (Vec<Vec<u64>>, Vec<u8>, Vec<u8>) =
        //     batch_result.into_iter().multiunzip();

        // // Check if minimizers match the index
        // // Separated from initial par_iter to allow flexibility with local/server processing
        // let batch_should_outputs = check_inputs_should_be_output(
        //     minimizer_hashes,
        //     &batch_minimizers,
        //     abs_threshold,
        //     rel_threshold,
        //     server_address,
        //     deplete,
        // );

        // Process results sequentially to maintain order
        for (i, (should_output, hit_count, total_minimizers, hit_kmers)) in
            batch_should_outputs.into_iter().enumerate()
        {
            // for (i, result) in batch_results.into_iter().enumerate() {
            let (record1, record2) = &batch_pairs[i];
            let seq1_len = record1.seq.len();
            let seq2_len = record2.seq.len();

            *total_seqs += 2;
            *total_bp += (seq1_len + seq2_len) as u64;

            if debug && hit_count > 0 {
                eprintln!(
                    "DEBUG: {}/{} hits={}/{} keep={} kmers=[{}]",
                    String::from_utf8_lossy(&record1.id),
                    String::from_utf8_lossy(&record2.id),
                    hit_count,
                    total_minimizers,
                    should_output,
                    hit_kmers.join(",")
                );
            }

            if should_output {
                // Track output base pairs
                *output_bp += (seq1_len + seq2_len) as u64;

                // Increment output sequence counter (twice, once for each seq)
                *output_seq_counter += 2;

                // Format and write record 1
                output_record_buffer.clear();
                output_fastx_record_from_parts(
                    &record1.id,
                    &record1.seq,
                    record1.qual.as_deref(),
                    record1.format,
                    &mut output_record_buffer,
                    rename,
                    *output_seq_counter - 1,
                );

                if let Some(ref mut w2) = writer2 {
                    // Write read 1 to primary writer
                    writer.write_all(&output_record_buffer)?;

                    // Format and write record 2 to second writer
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        &record2.id,
                        &record2.seq,
                        record2.qual.as_deref(),
                        record2.format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    w2.write_all(&output_record_buffer)?;
                } else {
                    // Interleaved output (existing behavior)
                    writer.write_all(&output_record_buffer)?;

                    // Format and write record 2
                    output_record_buffer.clear();
                    output_fastx_record_from_parts(
                        &record2.id,
                        &record2.seq,
                        record2.qual.as_deref(),
                        record2.format,
                        &mut output_record_buffer,
                        rename,
                        *output_seq_counter,
                    );
                    writer.write_all(&output_record_buffer)?;
                }
            } else {
                *filtered_seqs += 2; // Both seqs filtered out
                *filtered_bp += (seq1_len + seq2_len) as u64; // Track filtered base pairs
            }
        }

        // Update spinner and flush periodically
        let elapsed = filtering_start_time.elapsed();
        let seqs_per_sec = *total_seqs as f64 / elapsed.as_secs_f64();
        let bp_per_sec = *total_bp as f64 / elapsed.as_secs_f64();
        let mbp_per_sec = bp_per_sec / 1_000_000.0;

        // Calculate output proportion directly
        let output_seqs = *total_seqs - *filtered_seqs;
        let output_proportion = if *total_seqs > 0 {
            output_seqs as f64 / *total_seqs as f64
        } else {
            0.0
        };

        // Calculate output base pair proportion
        let output_bp_proportion = if *total_bp > 0 {
            *output_bp as f64 / *total_bp as f64
        } else {
            0.0
        };

        // Update spinner message with detailed stats
        spinner.set_message(format!(
            "Retained {}/{} seqs ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
            output_seqs,
            total_seqs,
            output_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0,
            seqs_per_sec,
            mbp_per_sec
        ));

        // Flush writer periodically
        writer.flush()?;
        if let Some(ref mut w2) = writer2 {
            w2.flush()?;
        }

        // Check if we've reached the end of input
        if reached_end {
            break;
        }
    }

    Ok(())
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
    use crate::index::IndexHeader;
    use crate::index::write_minimizers;
    use std::path::PathBuf;
    use tempfile::TempDir;

    #[allow(dead_code)] // Suppress unused warnings
    fn create_test_index() -> (PathBuf, IndexHeader, TempDir) {
        // Create a temporary directory
        let temp_dir = TempDir::new().unwrap();
        let index_path = temp_dir.path().join("test.idx");

        // Create dummy minimizers
        let minimizers: FxHashSet<u64> = [1, 2, 3, 4, 5].iter().cloned().collect();
        let header = IndexHeader::new(5, 3);

        write_minimizers(&minimizers, &header, Some(&index_path)).unwrap();

        // Return the TempDir along with the other values to keep it in scope
        (index_path, header, temp_dir)
    }

    #[test]
    fn test_filter_summary() {
        // Create a sample summary
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
        };

        // Test JSON ser+de
        let json = serde_json::to_string(&summary).unwrap();
        let parsed: FilterSummary = serde_json::from_str(&json).unwrap();

        // Check values
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
