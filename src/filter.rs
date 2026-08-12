use crate::index::{load_index_from_path_auto, load_minimizers_cached};
use crate::{
    ComplexityAlgorithm, FilterConfig, FilterDecision, FilterKernel, FilterParams, IndexHeader,
    MinimizerSet, validate_unit_interval,
};
use anyhow::{Context, Result};
use binseq::cbq;
use binseq::write::{BinseqWriterBuilder, Format as BinseqFormat};
use binseq::{BinseqRecord, ParallelReader as BinseqParallelReader, SequencingRecordBuilder};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use paraseq::Record;
use paraseq::fastx::Reader;
use paraseq::parallel::{PairedParallelProcessor, ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer
const DEFAULT_BUFFER_SIZE: usize = 64 * 1024;
/// Default CBQ block size in MiB; a block must hold the largest single record
pub const DEFAULT_CBQ_BLOCK_SIZE_MIB: u16 = 16;

type BoxedWriter = Box<dyn Write + Send>;
/// CBQ output is always a named file, no stdout
type CbqWriter = binseq::BinseqWriter<BufWriter<File>>;

/// Sequence file format
#[derive(Clone, Copy, PartialEq, Eq)]
enum Format {
    Fastx,
    Cbq,
}

/// Input metadata, resolved before any output is opened
struct InputLayout {
    format: Format,
    paired: bool,
    qualities: bool,
    headers: bool,
    /// CBQ per-record flags present (always false for FASTX)
    flags: bool,
    /// CBQ input block size, a lower bound for CBQ output (None for FASTX)
    block_size: Option<usize>,
}

#[allow(clippy::large_enum_variant)]
enum Input {
    Cbq(cbq::MmapReader),
    FastxSingle(Reader<Box<dyn std::io::Read + Send>>),
    FastxInterleaved(Reader<Box<dyn std::io::Read + Send>>),
    FastxPaired(
        Reader<Box<dyn std::io::Read + Send>>,
        Reader<Box<dyn std::io::Read + Send>>,
    ),
    Empty,
}

/// Borrowed view of one read, shared by every reader and writer
struct ReadView<'a> {
    id: &'a [u8],
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
    /// CBQ per-record flag, preserved on CBQ output (None elsewhere)
    flag: Option<u64>,
}

/// Format-specific output state: thread-local buffers plus the shared
/// writers each processor clone merges into
#[allow(clippy::large_enum_variant)] // mirrors binseq's BinseqWriter, which is not boxed upstream
#[derive(Clone)]
enum Output {
    Fastx {
        local: Vec<u8>,
        local2: Vec<u8>,
        pending: Vec<PendingRename>,
        pending2: Vec<PendingRename>,
        batch_kept: u64,
        shared: Arc<Mutex<BoxedWriter>>,
        shared2: Option<Arc<Mutex<BoxedWriter>>>,
    },
    Cbq {
        /// Thread-local block writer; completed blocks merge into the shared writer
        local: binseq::BinseqWriter<Vec<u8>>,
        shared: Arc<Mutex<CbqWriter>>,
    },
}

impl Output {
    fn fastx(writer: BoxedWriter, writer2: Option<BoxedWriter>) -> Self {
        Output::Fastx {
            local: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
            local2: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
            pending: Vec::new(),
            pending2: Vec::new(),
            batch_kept: 0,
            shared: Arc::new(Mutex::new(writer)),
            shared2: writer2.map(|w| Arc::new(Mutex::new(w))),
        }
    }

    fn cbq(writer: CbqWriter) -> Result<Self> {
        let shared = Arc::new(Mutex::new(writer));
        let local = shared.lock().new_headless_buffer()?;
        Ok(Output::Cbq { local, shared })
    }

    /// Format one read into a FASTX local buffer, deferring the header to
    /// [`write_renamed`] when renaming
    fn push_fastx(
        buffer: &mut Vec<u8>,
        pending: &mut Vec<PendingRename>,
        batch_kept: u64,
        read: &ReadView,
        suffix: &'static [u8],
        rename: bool,
        output_fasta: bool,
    ) -> Result<()> {
        let offset = buffer.len();
        let marker = format_record_to_buffer(read, rename, output_fasta, buffer)?;
        if rename {
            pending.push(PendingRename {
                offset,
                ordinal: batch_kept,
                marker,
                suffix,
            });
        }
        Ok(())
    }

    /// Write one unpaired read
    fn push_read(
        &mut self,
        read: &ReadView,
        rename: bool,
        output_fasta: bool,
        rename_counter: &AtomicU64,
    ) -> Result<()> {
        match self {
            Output::Fastx {
                local,
                pending,
                batch_kept,
                ..
            } => {
                Self::push_fastx(local, pending, *batch_kept, read, b"", rename, output_fasta)?;
                if rename {
                    *batch_kept += 1;
                }
            }
            Output::Cbq { local, .. } => {
                let number = rename.then(|| rename_counter.fetch_add(1, Ordering::Relaxed) + 1);
                let header = cbq_header(read.id, number, b"");
                let seq = sanitise_for_cbq(read.seq);
                let seq_record = SequencingRecordBuilder::default()
                    .s_seq(&seq)
                    .s_header(&header)
                    .opt_s_qual(read.qual)
                    .opt_flag(read.flag)
                    .build()?;
                local.push(seq_record)?;
            }
        }
        Ok(())
    }

    /// Write a kept pair: one native CBQ record, two FASTX files, or
    /// interleaved FASTX. Both mates take the pair's number when renaming.
    fn push_pair(
        &mut self,
        read1: &ReadView,
        read2: &ReadView,
        rename: bool,
        output_fasta: bool,
        rename_counter: &AtomicU64,
    ) -> Result<()> {
        match self {
            Output::Fastx {
                local,
                local2,
                pending,
                pending2,
                batch_kept,
                shared2,
                ..
            } => {
                Self::push_fastx(
                    local,
                    pending,
                    *batch_kept,
                    read1,
                    b"/1",
                    rename,
                    output_fasta,
                )?;
                if shared2.is_some() {
                    // Separate outputs
                    Self::push_fastx(
                        local2,
                        pending2,
                        *batch_kept,
                        read2,
                        b"/2",
                        rename,
                        output_fasta,
                    )?;
                } else {
                    // Interleaved output
                    Self::push_fastx(
                        local,
                        pending,
                        *batch_kept,
                        read2,
                        b"/2",
                        rename,
                        output_fasta,
                    )?;
                }
                if rename {
                    *batch_kept += 1;
                }
            }
            Output::Cbq { local, .. } => {
                let number = rename.then(|| rename_counter.fetch_add(1, Ordering::Relaxed) + 1);
                let header1 = cbq_header(read1.id, number, b"/1");
                let header2 = cbq_header(read2.id, number, b"/2");
                let seq1 = sanitise_for_cbq(read1.seq);
                let seq2 = sanitise_for_cbq(read2.seq);
                let seq_record = SequencingRecordBuilder::default()
                    .s_seq(&seq1)
                    .s_header(&header1)
                    .opt_s_qual(read1.qual)
                    .x_seq(&seq2)
                    .x_header(&header2)
                    .opt_x_qual(read2.qual)
                    .opt_flag(read1.flag)
                    .build()?;
                local.push(seq_record)?;
            }
        }
        Ok(())
    }

    /// Merge thread-local buffers into the shared writer(s) (per batch). Rename
    /// numbers are claimed under the writer lock: monotone in output order, but
    /// deterministic only with `--ordered` or `-t 1`.
    fn flush_batch(
        &mut self,
        rename: bool,
        rename_counter: &AtomicU64,
        ordered: bool,
    ) -> Result<()> {
        match self {
            Output::Cbq { local, shared } => {
                if ordered {
                    // Drain the incomplete block too, or records linger across batches
                    shared.lock().ingest(local)?;
                } else {
                    shared.lock().ingest_completed(local)?;
                }
            }
            Output::Fastx {
                local,
                local2,
                pending,
                pending2,
                batch_kept,
                shared,
                shared2,
            } => {
                if let Some(shared2) = shared2 {
                    // Atomic paired batch writing
                    if !local.is_empty() || !local2.is_empty() {
                        let mut writer1 = shared.lock();
                        let mut writer2 = shared2.lock();

                        if rename {
                            // Both mates share one block of numbers
                            let base = rename_counter.fetch_add(*batch_kept, Ordering::Relaxed) + 1;
                            write_renamed(&mut writer1, local, pending, base)?;
                            write_renamed(&mut writer2, local2, pending2, base)?;
                        } else {
                            writer1.write_all(local)?;
                            writer2.write_all(local2)?;
                        }
                        writer1.flush()?;
                        writer2.flush()?;
                    }
                } else if !local.is_empty() {
                    let mut writer = shared.lock();
                    if rename {
                        let base = rename_counter.fetch_add(*batch_kept, Ordering::Relaxed) + 1;
                        write_renamed(&mut writer, local, pending, base)?;
                    } else {
                        writer.write_all(local)?;
                    }
                    writer.flush()?;
                }
                local.clear();
                local2.clear();
                pending.clear();
                pending2.clear();
                *batch_kept = 0;
            }
        }
        Ok(())
    }

    /// Merge any remaining local CBQ blocks into the shared writer (per thread)
    fn flush_thread(&mut self) -> Result<()> {
        if let Output::Cbq { local, shared } = self {
            shared.lock().ingest(local)?;
        }
        Ok(())
    }

    /// Finish a CBQ stream: flush remaining blocks and write the embedded index
    fn finish(&self) -> Result<()> {
        if let Output::Cbq { shared, .. } = self {
            shared.lock().finish()?;
        }
        Ok(())
    }
}

/// Header for a CBQ record: the original id, or `number[ suffix]` when
/// renaming. Numbers are claimed per record: unique, but sequential only
/// with `--ordered` or `-t 1`.
fn cbq_header<'a>(id: &'a [u8], number: Option<u64>, suffix: &[u8]) -> Cow<'a, [u8]> {
    match number {
        Some(n) => {
            let mut header = n.to_string().into_bytes();
            if !suffix.is_empty() {
                header.push(b' ');
                header.extend_from_slice(suffix);
            }
            Cow::Owned(header)
        }
        None => Cow::Borrowed(id),
    }
}

/// Map unsupported bases to N.
fn sanitise_for_cbq(seq: &[u8]) -> Cow<'_, [u8]> {
    #[inline]
    fn supported(b: u8) -> bool {
        matches!(
            b,
            b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't'
        )
    }

    let Some(first) = seq.iter().position(|&b| !supported(b)) else {
        return Cow::Borrowed(seq);
    };
    let mut sanitised = seq.to_vec();
    sanitised[first..].iter_mut().for_each(|b| {
        if !supported(*b) {
            *b = b'N';
        }
    });
    Cow::Owned(sanitised)
}

/// Filtering config for an already-loaded index (no index path; see [`FilterConfig`]).
pub struct FilterRunConfig {
    /// Path to input fastx file (or - for stdin)
    pub input_path: String,
    /// Path to optional second paired fastx file (or - for interleaved stdin)
    pub input2_path: Option<String>,
    /// Treat input_path as an interleaved paired stream
    pub interleaved: bool,
    /// Validate paired record names (Illumina CASAVA or /1 /2 suffixes)
    pub check_pairs: bool,
    /// Path to output fastx file (None for stdout; detects .gz/.zst/.xz)
    pub output_path: Option<PathBuf>,
    /// Path to optional second output fastx file for paired reads
    pub output2_path: Option<String>,
    /// Absolute threshold for filtering sequences
    pub abs_threshold: usize,
    /// Relative threshold for filtering sequences (0.0-1.0)
    pub rel_threshold: f64,
    /// Consider only the first N nucleotides per sequence (0 = entire sequence)
    pub prefix_length: usize,
    /// Path to JSON summary file (None to skip writing one; stats are always returned)
    pub summary_path: Option<PathBuf>,
    /// Deplete mode (remove sequences WITH matches)
    pub deplete: bool,
    /// Replace sequence headers with incrementing numbers
    pub rename: bool,
    /// Force FASTA output (discards quality scores)
    pub output_fasta: bool,
    /// Preserve input record ordering (deterministic, slightly slower)
    pub ordered: bool,
    /// Number of execution threads (0 = auto)
    pub threads: u16,
    /// Compression level for output files (1-22 for zst, 1-9 for gz)
    pub compression_level: u8,
    /// cbq output block size in MiB (raised to the cbq input's block size if larger)
    pub cbq_block_size: u16,
    /// Number of threads for compression (0 = auto)
    pub compression_threads: u16,
    /// Debug mode: output sequences with minimizer hits to stderr
    pub debug: bool,
    /// Suppress progress reporting
    pub quiet: bool,
    /// Label recorded in the summary's `index` field (no filesystem check)
    pub index_label: String,
}

/// Config for FilterProcessor
struct FilterProcessorConfig {
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    deplete: bool,
    rename: bool,
    output_fasta: bool,
    debug: bool,
    check_pairs: bool,
    ordered: bool,
}

/// Split a FASTA/Q header into its first whitespace-delimited token and description
#[inline]
fn split_record_id(id: &[u8]) -> (&[u8], &[u8]) {
    match id.iter().position(|byte| byte.is_ascii_whitespace()) {
        Some(index) => (&id[..index], &id[index..]),
        None => (id, &[]),
    }
}

#[inline]
fn casava_mate_number(description: &[u8]) -> Option<u8> {
    let description = description
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .map_or(&[][..], |index| &description[index..]);

    match description {
        [mate @ (b'1' | b'2'), b':', ..] => Some(*mate),
        _ => None,
    }
}

#[inline]
fn paired_record_names_match(id1: &[u8], id2: &[u8]) -> bool {
    let (name1, description1) = split_record_id(id1);
    let (name2, description2) = split_record_id(id2);

    if let (Some(core1), Some(core2)) = (name1.strip_suffix(b"/1"), name2.strip_suffix(b"/2")) {
        return !core1.is_empty() && core1 == core2;
    }

    !name1.is_empty()
        && name1 == name2
        && casava_mate_number(description1) == Some(b'1')
        && casava_mate_number(description2) == Some(b'2')
}

fn validate_check_pairs_mode(check_pairs: bool, paired_input: bool) -> Result<()> {
    if check_pairs && !paired_input {
        anyhow::bail!("--check-pairs requires paired input (INPUT2 or --interleaved)");
    }
    Ok(())
}

/// Check if path is a named pipe or process substitution / /dev/fd/*
#[cfg(unix)]
fn is_special_input_path(path: &str) -> bool {
    use std::os::unix::fs::FileTypeExt;
    path.starts_with("/dev/fd/")
        || path.starts_with("/proc/self/fd/")
        || std::path::Path::new(path)
            .metadata()
            .map(|m| m.file_type().is_fifo())
            .unwrap_or(false)
}

#[cfg(not(unix))]
fn is_special_input_path(_path: &str) -> bool {
    false
}

/// Check input fastx file path(s) exist (the index is already loaded, so not checked here)
fn check_input_paths(config: &FilterRunConfig) -> Result<()> {
    if config.input_path != "-"
        && !is_special_input_path(&config.input_path)
        && !std::path::Path::new(&config.input_path).exists()
    {
        return Err(anyhow::anyhow!(
            "Input file does not exist: {}",
            config.input_path
        ));
    }

    if let Some(input2_path) = &config.input2_path
        && input2_path != "-"
        && !is_special_input_path(input2_path)
        && !std::path::Path::new(input2_path).exists()
    {
        return Err(anyhow::anyhow!(
            "Second input file does not exist: {}",
            input2_path
        ));
    }

    Ok(())
}

/// Check if file metadata len < 5 (catches empty uncompressed files only)
fn is_empty_file(path: &str) -> Result<bool> {
    if path == "-" || is_special_input_path(path) {
        return Ok(false);
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

/// Resolve the output format from the output path suffix
fn resolve_output_format(config: &FilterRunConfig) -> Format {
    match config.output_path.as_deref() {
        Some(path) if path.to_string_lossy().ends_with(".cbq") => Format::Cbq,
        _ => Format::Fastx,
    }
}

/// Resolve the input format and layout, opening all readers up front so input
/// errors surface before any output file is created.
fn open_input(config: &FilterRunConfig, interleaved_input: bool) -> Result<(InputLayout, Input)> {
    // CBQ input is file-only (mmap reader); stdin would need binseq's streaming reader
    if config.input_path == "-" || is_special_input_path(&config.input_path) {
        return open_fastx(config, interleaved_input);
    }

    // Regular files: sniff the CBQ magic; everything else is FASTX
    let file = File::open(&config.input_path)
        .map_err(|e| anyhow::anyhow!("Failed to open file {}: {}", config.input_path, e))?;
    let file_len = file.metadata()?.len();
    let mut magic = Vec::with_capacity(64);
    file.take(64).read_to_end(&mut magic)?;
    match BinseqFormat::sniff(&magic) {
        Some(BinseqFormat::Cbq) => {
            // The mmap reader slices header and footer unchecked, so guard the length
            let min_len = (std::mem::size_of::<cbq::FileHeader>()
                + std::mem::size_of::<cbq::IndexFooter>()) as u64;
            if file_len < min_len {
                anyhow::bail!("Truncated or corrupt CBQ input: {}", config.input_path);
            }
            let reader =
                cbq::MmapReader::new(&config.input_path).context("Failed to open CBQ input")?;
            let header = reader.header();
            let layout = InputLayout {
                format: Format::Cbq,
                paired: header.is_paired(),
                qualities: header.has_qualities(),
                headers: header.has_headers(),
                flags: header.has_flags(),
                block_size: Some(header.block_size as usize),
            };
            // binseq's parallel reader rejects an empty record range
            let input = if reader.num_records() == 0 {
                Input::Empty
            } else {
                Input::Cbq(reader)
            };
            Ok((layout, input))
        }
        Some(f) => anyhow::bail!("{f:?} input is not supported"),
        None => open_fastx(config, interleaved_input),
    }
}

/// Open a FASTX input and resolve its layout
fn open_fastx(config: &FilterRunConfig, interleaved_input: bool) -> Result<(InputLayout, Input)> {
    let layout = InputLayout {
        format: Format::Fastx,
        paired: interleaved_input || config.input2_path.is_some(),
        qualities: false,
        headers: true,
        flags: false,
        block_size: None,
    };

    let input1_empty = is_empty_file(&config.input_path)?;
    let input2_empty = config
        .input2_path
        .as_deref()
        .map(is_empty_file)
        .transpose()?
        .unwrap_or(false);

    if interleaved_input {
        if input1_empty {
            return Ok((layout, Input::Empty));
        }
        return match create_paraseq_reader(Some(config.input_path.as_str())) {
            Ok(reader) => {
                let qualities = reader.format() == paraseq::fastx::Format::Fastq;
                Ok((
                    InputLayout {
                        qualities,
                        ..layout
                    },
                    Input::FastxInterleaved(reader),
                ))
            }
            Err(e) if is_empty_input_error(&e) => Ok((layout, Input::Empty)),
            Err(e) => Err(e),
        };
    }

    if let Some(input2_path) = config.input2_path.as_deref() {
        if input1_empty && input2_empty {
            return Ok((layout, Input::Empty));
        }
        if input1_empty || input2_empty {
            return Err(anyhow::anyhow!(
                "One paired file is empty but the other is not"
            ));
        }
        let r1 = create_paraseq_reader(Some(config.input_path.as_str()));
        let r2 = create_paraseq_reader(Some(input2_path));
        return match (r1, r2) {
            (Ok(reader1), Ok(reader2)) => {
                let qualities = reader1.format() == paraseq::fastx::Format::Fastq;
                Ok((
                    InputLayout {
                        qualities,
                        ..layout
                    },
                    Input::FastxPaired(reader1, reader2),
                ))
            }
            (Err(e1), Err(e2)) if is_empty_input_error(&e1) && is_empty_input_error(&e2) => {
                Ok((layout, Input::Empty))
            }
            (Err(e), _) if is_empty_input_error(&e) => Err(anyhow::anyhow!(
                "First paired file appears empty while second is not"
            )),
            (_, Err(e)) if is_empty_input_error(&e) => Err(anyhow::anyhow!(
                "Second paired file appears empty while first is not"
            )),
            (Err(e), _) => Err(e),
            (_, Err(e)) => Err(e),
        };
    }

    if input1_empty {
        return Ok((layout, Input::Empty));
    }
    match create_paraseq_reader(Some(config.input_path.as_str())) {
        Ok(reader) => {
            let qualities = reader.format() == paraseq::fastx::Format::Fastq;
            Ok((
                InputLayout {
                    qualities,
                    ..layout
                },
                Input::FastxSingle(reader),
            ))
        }
        Err(e) if is_empty_input_error(&e) => Ok((layout, Input::Empty)),
        Err(e) => Err(e),
    }
}

/// Validate format combinations before any output file is opened or truncated
fn validate_input_output(
    layout: &InputLayout,
    output_format: Format,
    config: &FilterRunConfig,
) -> Result<()> {
    if layout.format == Format::Cbq && config.input2_path.is_some() {
        anyhow::bail!("CBQ input does not support INPUT2");
    }
    if layout.format == Format::Cbq && config.interleaved {
        anyhow::bail!("CBQ input does not support --interleaved");
    }
    if output_format == Format::Cbq && config.output2_path.is_some() {
        anyhow::bail!("CBQ output does not support OUTPUT2; CBQ pairing is native");
    }
    if config
        .output2_path
        .as_deref()
        .is_some_and(|output2| output2.ends_with(".cbq"))
    {
        anyhow::bail!("OUTPUT2 cannot be CBQ; CBQ pairing is native (use a single --output)");
    }
    if output_format == Format::Cbq && !(1..=22).contains(&config.compression_level) {
        anyhow::bail!(
            "Invalid CBQ compression level {} (must be 1-22)",
            config.compression_level
        );
    }
    if config.check_pairs && layout.format == Format::Cbq && !layout.headers {
        anyhow::bail!("--check-pairs requires CBQ input with headers");
    }
    validate_check_pairs_mode(config.check_pairs, layout.paired)?;
    Ok(())
}

/// A record buffered without its header, which is written later by [`write_renamed`]
/// once its final number is known (numbering depends on earlier batches finishing).
#[derive(Clone)]
struct PendingRename {
    /// Start of the headerless body in the local buffer
    offset: usize,
    /// Number relative to the first of its batch
    ordinal: u64,
    /// `>` for FASTA, `@` for FASTQ
    marker: u8,
    /// Mate suffix after the number, e.g. `/1`; empty when unpaired
    suffix: &'static [u8],
}

/// Format a record into a buffer (FASTA/FASTQ), returning its line prefix
/// `seq` is the newline-stripped sequence from `record.seq()`.
/// When renaming, the header is left to [`write_renamed`].
fn format_record_to_buffer(
    read: &ReadView,
    rename: bool,
    output_fasta: bool,
    buffer: &mut Vec<u8>,
) -> Result<u8> {
    let is_fasta = output_fasta || read.qual.is_none();
    let marker = if is_fasta { b'>' } else { b'@' };

    // Header (omitted when renaming)
    if !rename {
        buffer.push(marker);
        buffer.extend_from_slice(read.id);
        buffer.write_all(b"\n")?;
    }

    // Sequence
    buffer.extend_from_slice(read.seq);

    if is_fasta {
        buffer.write_all(b"\n")?;
    } else {
        // FASTQ: plus and qual lines
        buffer.write_all(b"\n+\n")?;
        if let Some(qual) = read.qual {
            buffer.extend_from_slice(qual);
        }
        buffer.write_all(b"\n")?;
    }
    Ok(marker)
}

/// Write buffered records to `writer`, numbering them from `base`
///
/// Each body runs from its own offset to the next one's, or to the end of `buffer`.
fn write_renamed(
    writer: &mut BoxedWriter,
    buffer: &[u8],
    pending: &[PendingRename],
    base: u64,
) -> Result<()> {
    for (i, record) in pending.iter().enumerate() {
        let end = pending.get(i + 1).map_or(buffer.len(), |next| next.offset);
        writer.write_all(&[record.marker])?;
        writer.write_all((base + record.ordinal).to_string().as_bytes())?;
        if !record.suffix.is_empty() {
            writer.write_all(b" ")?;
            writer.write_all(record.suffix)?;
        }
        writer.write_all(b"\n")?;
        writer.write_all(&buffer[record.offset..end])?;
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
fn count_compressed_outputs(config: &FilterRunConfig) -> u8 {
    let mut count = 0;
    if is_compressed_output(config.output_path.as_deref()) {
        count += 1;
    }
    if let Some(output2) = &config.output2_path
        && is_compressed_output(Some(std::path::Path::new(output2)))
    {
        count += 1;
    }
    count
}

/// Open a CBQ output file directly, bypassing the compressed writer stack
fn open_cbq_output(path: &std::path::Path) -> Result<BufWriter<File>> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create output file: {}", path.display()))?;
    Ok(BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file))
}

/// Return a suitable writer for the output path extension
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
            // `auto_finish()` yields a writer that writes the zstd frame
            // epilogue on drop. Without it, dropping a bare `Encoder` closes
            // the file without finalizing the frame, producing a truncated
            // `.zst` stream (`zstd -t` reports "premature end").
            Ok(Box::new(
                zstd::stream::write::Encoder::new(buffered_file, compression_level as i32)?
                    .auto_finish(),
            ))
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
    ordered: bool,
    check_pairs: bool,
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
    minimizers: Arc<MinimizerSet>,
    rename: bool,
    output_fasta: bool,
    debug: bool,
    check_pairs: bool,
    /// Write batches in input order, not completion order
    ordered: bool,
    kernel: FilterKernel,

    output: Output,
    local_stats: ProcessingStats,

    /// Shared across workers, handing out rename numbers
    rename_counter: Arc<AtomicU64>,

    // Global state
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
    pub last_reported: u64,
}

impl FilterProcessor {
    fn new(
        minimizers: Arc<MinimizerSet>,
        kmer_length: u8,
        window_size: u8,
        config: &FilterProcessorConfig,
        output: Output,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        filtering_start_time: Instant,
    ) -> Result<Self> {
        Ok(Self {
            minimizers,
            rename: config.rename,
            output_fasta: config.output_fasta,
            debug: config.debug,
            check_pairs: config.check_pairs,
            ordered: config.ordered,
            kernel: FilterKernel::new(
                kmer_length,
                window_size,
                FilterParams {
                    deplete: config.deplete,
                    abs_threshold: config.abs_threshold,
                    rel_threshold: config.rel_threshold,
                    prefix_length: config.prefix_length,
                },
            )?,
            output,
            local_stats: ProcessingStats::default(),
            rename_counter: Arc::new(AtomicU64::new(0)),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            spinner,
            filtering_start_time,
        })
    }

    fn should_keep_sequence(&mut self, seq: &[u8]) -> FilterDecision {
        self.kernel.classify_read(&self.minimizers, seq, self.debug)
    }

    fn should_keep_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> FilterDecision {
        self.kernel
            .classify_pair(&self.minimizers, seq1, seq2, self.debug)
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

    /// Shared per-read logic for every reader (paraseq FASTX, CBQ mmap)
    fn handle_read(&mut self, read: &ReadView) -> Result<()> {
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += read.seq.len() as u64;

        let decision = self.should_keep_sequence(read.seq);

        // Show debug info for sequences with hits
        if self.debug {
            eprintln!(
                "DEBUG: {} hits={}/{} keep={} kmers=[{}]",
                String::from_utf8_lossy(read.id),
                decision.hit_count,
                decision.total_minimizers,
                decision.keep,
                decision.hit_kmers.join(",")
            );
        }

        if decision.keep {
            self.local_stats.output_bp += read.seq.len() as u64;
            self.output
                .push_read(read, self.rename, self.output_fasta, &self.rename_counter)?;
        } else {
            self.local_stats.filtered_seqs += 1;
            self.local_stats.filtered_bp += read.seq.len() as u64;
        }

        Ok(())
    }

    /// Shared per-pair logic for every reader (paraseq FASTX, CBQ mmap)
    fn handle_pair(&mut self, read1: &ReadView, read2: &ReadView) -> Result<()> {
        if self.check_pairs && !paired_record_names_match(read1.id, read2.id) {
            return Err(anyhow::anyhow!(
                "Paired record name mismatch: R1='{}', R2='{}'. Expected matching Illumina CASAVA 1: and 2: fields or names suffixed with /1 and /2",
                String::from_utf8_lossy(read1.id),
                String::from_utf8_lossy(read2.id)
            ));
        }

        self.local_stats.total_seqs += 2;
        self.local_stats.total_bp += (read1.seq.len() + read2.seq.len()) as u64;

        let decision = self.should_keep_pair(read1.seq, read2.seq);

        // Debug info for interleaved pairs
        if self.debug && decision.hit_count > 0 {
            eprintln!(
                "DEBUG: {}/{} hits={}/{} keep={} kmers=[{}]",
                String::from_utf8_lossy(read1.id),
                String::from_utf8_lossy(read2.id),
                decision.hit_count,
                decision.total_minimizers,
                decision.keep,
                decision.hit_kmers.join(",")
            );
        }

        if decision.keep {
            self.local_stats.output_bp += (read1.seq.len() + read2.seq.len()) as u64;
            self.output.push_pair(
                read1,
                read2,
                self.rename,
                self.output_fasta,
                &self.rename_counter,
            )?;
        } else {
            self.local_stats.filtered_seqs += 2;
            self.local_stats.filtered_bp += (read1.seq.len() + read2.seq.len()) as u64;
        }

        Ok(())
    }

    fn flush_batch(&mut self) -> Result<()> {
        self.output
            .flush_batch(self.rename, &self.rename_counter, self.ordered)?;

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.filtered_seqs += self.local_stats.filtered_seqs;
            stats.total_bp += self.local_stats.total_bp;
            stats.output_bp += self.local_stats.output_bp;
            stats.filtered_bp += self.local_stats.filtered_bp;
        }

        // Update spinner
        self.update_spinner();

        // Reset local stats
        self.local_stats = ProcessingStats::default();

        Ok(())
    }

    fn flush_thread(&mut self) -> Result<()> {
        self.output.flush_thread()
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for FilterProcessor {
    fn requires_ordering(&self) -> bool {
        self.ordered
    }

    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq = record.seq();
        self.handle_read(&ReadView {
            id: record.id(),
            seq: &seq,
            qual: record.qual(),
            flag: None,
        })?;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        self.flush_batch()?;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::parallel::Result<()> {
        self.flush_thread()?;
        Ok(())
    }
}

impl<Rf: Record> PairedParallelProcessor<Rf> for FilterProcessor {
    fn requires_ordering(&self) -> bool {
        self.ordered
    }

    fn process_record_pair(&mut self, record1: Rf, record2: Rf) -> paraseq::parallel::Result<()> {
        let seq1 = record1.seq();
        let seq2 = record2.seq();
        self.handle_pair(
            &ReadView {
                id: record1.id(),
                seq: &seq1,
                qual: record1.qual(),
                flag: None,
            },
            &ReadView {
                id: record2.id(),
                seq: &seq2,
                qual: record2.qual(),
                flag: None,
            },
        )?;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        self.flush_batch()?;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::parallel::Result<()> {
        self.flush_thread()?;
        Ok(())
    }
}

impl binseq::ParallelProcessor for FilterProcessor {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        let read1 = ReadView {
            id: record.sheader(),
            seq: record.sseq(),
            qual: record.has_quality().then(|| record.squal()),
            flag: record.flag(),
        };
        if record.is_paired() {
            let read2 = ReadView {
                id: record.xheader(),
                seq: record.xseq(),
                qual: record.has_quality().then(|| record.xqual()),
                flag: record.flag(),
            };
            self.handle_pair(&read1, &read2)?;
        } else {
            self.handle_read(&read1)?;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        self.flush_batch()?;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.flush_thread()?;
        Ok(())
    }
}

pub fn run(config: &FilterConfig) -> Result<FilterSummary> {
    validate_unit_interval("relative threshold", config.rel_threshold)?;
    if let Some(threshold) = config.complexity_threshold {
        validate_unit_interval("complexity threshold", threshold)?;
    }

    // Fail --check-pairs before the index loads; CBQ may be natively paired,
    // so skip it there and let validate_input_output recheck once the header is read
    if config.check_pairs
        && !config.interleaved
        && config.input2_path.is_none()
        && !config.input_path.ends_with(".cbq")
    {
        validate_check_pairs_mode(true, false)?;
    }

    // Validate the index path once here; run_with_index never touches it again.
    if !config.minimizers_path.exists() {
        return Err(anyhow::anyhow!(
            "Index file does not exist: {}",
            config.minimizers_path.display()
        ));
    }

    let quiet = config.quiet || config.debug;
    let load_start = Instant::now();

    let run_config = FilterRunConfig {
        input_path: config.input_path.to_string(),
        input2_path: config.input2_path.map(str::to_string),
        interleaved: config.interleaved,
        check_pairs: config.check_pairs,
        output_path: config.output_path.map(|p| p.to_path_buf()),
        output2_path: config.output2_path.map(str::to_string),
        abs_threshold: config.abs_threshold,
        rel_threshold: config.rel_threshold,
        prefix_length: config.prefix_length,
        summary_path: config.summary_path.cloned(),
        deplete: config.deplete,
        rename: config.rename,
        output_fasta: config.output_fasta,
        ordered: config.ordered,
        threads: config.threads,
        compression_level: config.compression_level,
        cbq_block_size: config.cbq_block_size,
        compression_threads: config.compression_threads,
        debug: config.debug,
        quiet: config.quiet,
        index_label: config.minimizers_path.to_string_lossy().into_owned(),
    };

    // Discard low-complexity (kdust) index minimizers once at load
    if let Some(threshold) = config.complexity_threshold {
        let (mut minimizers, header) = load_index_from_path_auto(config.minimizers_path)?;
        if matches!(minimizers, MinimizerSet::Fuse(_)) {
            return Err(anyhow::anyhow!(
                "Complexity filtering is not supported on BFF indexes; use an exact index"
            ));
        }
        let before = minimizers.len();
        minimizers.retain_complexity(
            header.kmer_length(),
            ComplexityAlgorithm::Kdust,
            threshold,
            false,
        )?;
        if !quiet {
            eprintln!(
                "Loaded index (k={}, w={}) in {:.2?}; kept {} of {} minimizers (kdust >= {})",
                header.kmer_length(),
                header.window_size(),
                load_start.elapsed(),
                minimizers.len(),
                before,
                threshold
            );
        }
        return run_with_index(Arc::new(minimizers), &header, &run_config);
    }

    let (minimizers, header) = load_minimizers_cached(config.minimizers_path)?;
    if !quiet {
        eprintln!(
            "Loaded index (k={}, w={}) in {:.2?}",
            header.kmer_length(),
            header.window_size(),
            load_start.elapsed()
        );
    }

    run_with_index(minimizers, header, &run_config)
}

/// Filter an already-loaded index against input fastx file(s), returning summary stats.
///
/// Reusable entry point behind the Python bindings... load the index once, call repeatedly.
/// Does no index-path validation; the index is already in memory.
pub fn run_with_index(
    minimizers: Arc<MinimizerSet>,
    header: &IndexHeader,
    config: &FilterRunConfig,
) -> Result<FilterSummary> {
    validate_unit_interval("relative threshold", config.rel_threshold)?;

    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    let tool_version = format!("deacon {}", version);

    let quiet = config.quiet || config.debug;

    let kmer_length = header.kmer_length();
    let window_size = header.window_size();

    let total_threads = if config.threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        config.threads as usize
    };

    let compressed_output_count = count_compressed_outputs(config);

    // Allocate threads between filtering (rayon) and compression (gzp).
    // Rayon pool can only be initialised once, so calculate before build_global().
    let (filtering_threads, compression_threads_per_output) = if compressed_output_count > 0 {
        let compression_threads_total = if config.compression_threads > 0 {
            config.compression_threads as usize
        } else {
            total_threads.div_ceil(2) // Auto: ceil(total_threads / 2)
        };
        let filtering_threads = total_threads
            .saturating_sub(compression_threads_total)
            .max(1);
        let output_count = compressed_output_count as usize;
        let threads_per_output = compression_threads_total.div_ceil(output_count).max(1);
        (filtering_threads, threads_per_output)
    } else {
        (total_threads, 0)
    };

    if filtering_threads > 0 {
        // error is OK here when we initialise a 2nd time in server mode.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(filtering_threads)
            .build_global()
            .context("Failed to initialise thread pool");
    }

    check_input_paths(config)?;

    // Resolve formats and input metadata before opening or truncating outputs
    let interleaved_stdin = config.input_path == "-" && config.input2_path.as_deref() == Some("-");
    let interleaved_input = config.interleaved || interleaved_stdin;
    let output_format = resolve_output_format(config);
    let (layout, input) = open_input(config, interleaved_input)?;
    validate_input_output(&layout, output_format, config)?;

    // The binseq reader and CBQ rename numbering cannot honour --ordered
    // multithreaded
    let ordered_cbq = config.ordered
        && (layout.format == Format::Cbq || (output_format == Format::Cbq && config.rename));
    let filtering_threads = if ordered_cbq && filtering_threads > 1 {
        if !quiet {
            eprintln!("Using 1 filtering thread: --ordered with CBQ input or renamed CBQ output");
        }
        1
    } else {
        filtering_threads
    };

    let mode = if config.deplete { "deplete" } else { "search" };

    let mut input_type = String::new();
    let mut options = Vec::<String>::new();
    if interleaved_input {
        input_type.push_str("interleaved");
    } else if layout.paired {
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
    if config.ordered {
        options.push("ordered".to_string());
    }
    if config.check_pairs {
        options.push("check-pairs".to_string());
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

    let output = match output_format {
        Format::Cbq => {
            let path = config
                .output_path
                .as_deref()
                .expect("CBQ output implies a named path");
            let cbq_writer = BinseqWriterBuilder::new(BinseqFormat::Cbq)
                .paired(layout.paired)
                .quality(layout.qualities && !config.output_fasta)
                .headers(layout.headers || config.rename)
                .flags(layout.flags)
                .block_size(
                    (config.cbq_block_size as usize * 1024 * 1024)
                        .max(layout.block_size.unwrap_or(0)),
                )
                .compression_level(i32::from(config.compression_level))
                .build(open_cbq_output(path)?)
                .context("Failed to create CBQ writer")?;
            Output::cbq(cbq_writer)?
        }
        Format::Fastx => {
            let writer = get_writer(
                config.output_path.as_deref(),
                config.compression_level,
                compression_threads_per_output,
            )?;
            let writer2 = if let Some(output2) = config.output2_path.as_deref() {
                if layout.paired {
                    Some(get_writer(
                        Some(std::path::Path::new(output2)),
                        config.compression_level,
                        compression_threads_per_output,
                    )?)
                } else {
                    None
                }
            } else {
                None
            };
            Output::fastx(writer, writer2)
        }
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
        output_fasta: config.output_fasta,
        debug: config.debug,
        check_pairs: config.check_pairs,
        ordered: config.ordered,
    };
    let mut processor = FilterProcessor::new(
        minimizers,
        kmer_length,
        window_size,
        &processor_config,
        output,
        spinner.clone(),
        filtering_start_time,
    )?;

    // Process based on input type - use filtering threads (already calculated above)
    let num_threads = filtering_threads;

    match input {
        Input::Cbq(reader) => {
            reader.process_parallel(processor.clone(), num_threads)?;
        }
        Input::FastxSingle(reader) => {
            reader.process_parallel(&mut processor, num_threads)?;
        }
        Input::FastxInterleaved(reader) => {
            reader.process_parallel_interleaved(&mut processor, num_threads)?;
        }
        Input::FastxPaired(reader1, reader2) => {
            reader1.process_parallel_paired(reader2, &mut processor, num_threads)?;
        }
        Input::Empty => {
            if !quiet {
                eprintln!("Empty input file(s) detected");
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

    // Finish any CBQ stream (writes the embedded index), then drop the
    // processor so the writers flush
    processor.output.finish()?;
    drop(processor);

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

    let summary = FilterSummary {
        version: tool_version,
        index: config.index_label.clone(),
        input: config.input_path.clone(),
        input2: config.input2_path.clone(),
        output: config
            .output_path
            .as_deref()
            .map_or("-".to_string(), |p| p.display().to_string()),
        output2: config.output2_path.clone(),
        k: kmer_length,
        w: window_size,
        abs_threshold: config.abs_threshold,
        rel_threshold: config.rel_threshold,
        prefix_length: config.prefix_length,
        deplete: config.deplete,
        rename: config.rename,
        ordered: config.ordered,
        check_pairs: config.check_pairs,
        seqs_in: total_seqs,
        seqs_out: output_seqs,
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
        seqs_per_second_total: seqs_per_sec_total as u64,
        bp_per_second_total: bp_per_sec_total as u64,
    };

    if let Some(summary_file) = &config.summary_path {
        let file = File::create(summary_file)
            .context(format!("Failed to create summary: {:?}", summary_file))?;
        serde_json::to_writer_pretty(BufWriter::new(file), &summary)
            .context("Failed to write summary")?;
        if !quiet {
            eprintln!("Filter summary saved to {:?}", summary_file);
        }
    }

    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paired_record_names_match_supported_formats() {
        assert!(paired_record_names_match(b"cluster/1", b"cluster/2"));
        assert!(paired_record_names_match(
            b"cluster/1 legacy description",
            b"cluster/2 other description"
        ));
        assert!(paired_record_names_match(
            b"A00123:1:H5J2TDSX7:1:1101:1000:1000 1:N:0:ACGT",
            b"A00123:1:H5J2TDSX7:1:1101:1000:1000 2:Y:7:TGCA"
        ));
        assert!(paired_record_names_match(
            b"A00123:1:H5J2TDSX7:1:1101:1000:1000\t1:N:0:ACGT",
            b"A00123:1:H5J2TDSX7:1:1101:1000:1000\t 2:N:0:ACGT"
        ));
    }

    #[test]
    fn test_paired_record_names_reject_invalid_pairs() {
        let invalid = [
            (&b"other/1"[..], &b"cluster/2"[..]),
            (&b"cluster/2"[..], &b"cluster/1"[..]),
            (&b"cluster/1"[..], &b"cluster/1"[..]),
            (&b"cluster"[..], &b"cluster"[..]),
            (&b"cluster/1"[..], &b"cluster 2:N:0:1"[..]),
            (&b"/1"[..], &b"/2"[..]),
            (&b"cluster 1:N:0:1"[..], &b"other 2:N:0:1"[..]),
            (&b"cluster 2:N:0:1"[..], &b"cluster 1:N:0:1"[..]),
            (&b"cluster 1:N:0:1"[..], &b"cluster 1:N:0:1"[..]),
            (&b"cluster 1"[..], &b"cluster 2"[..]),
        ];

        for (id1, id2) in invalid {
            assert!(
                !paired_record_names_match(id1, id2),
                "accepted invalid pair"
            );
        }
    }

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
            ordered: false,
            check_pairs: false,
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

    #[cfg(feature = "compression")]
    fn decode_output(path: &std::path::Path) -> std::io::Result<Vec<u8>> {
        use std::io::Read;
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
        let file = std::fs::File::open(path)?;
        let mut out = Vec::new();
        match ext {
            "gz" => flate2::read::MultiGzDecoder::new(file).read_to_end(&mut out)?,
            "zst" => zstd::stream::read::Decoder::new(file)?.read_to_end(&mut out)?,
            "xz" => liblzma::read::XzDecoder::new(file).read_to_end(&mut out)?,
            _ => std::io::BufReader::new(file).read_to_end(&mut out)?,
        };
        Ok(out)
    }

    // Test regression for #88: dropped writers must leave complete output.
    #[cfg(feature = "compression")]
    #[rstest::rstest]
    #[case("out.fq.gz")]
    #[case("out.fq.zst")]
    #[case("out.fq.xz")]
    #[case("out.fq")]
    fn test_output_stream_is_finalized(#[case] filename: &str) {
        for payload in [
            b"@read0\nACGTACGTACGT\n+\nIIIIIIIIIIII\n".to_vec(),
            Vec::new(),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join(filename);

            {
                let mut writer = get_writer(Some(&path), 2, 1).unwrap();
                writer.write_all(&payload).unwrap();
                writer.flush().unwrap();
            }

            let decoded = decode_output(&path).unwrap_or_else(|e| {
                panic!("{filename} did not round-trip (truncated stream?): {e}")
            });
            assert_eq!(decoded, payload, "{filename} content mismatch");
        }
    }
}
