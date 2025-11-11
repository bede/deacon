use crate::filter::{Buffers, ProcessingStats};
use crate::minimizers::KmerHasher;
use crate::{FixedRapidHasher, IndexConfig, RapidHashSet};
use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use paraseq::Record;
use paraseq::prelude::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};
use std::time::Instant;

#[cfg(feature = "fetch")]
use indicatif::ProgressBar;
#[cfg(feature = "fetch")]
use std::time::Duration;

/// Index format version
pub const INDEX_FORMAT_VERSION: u8 = 3;

/// Serialisable header for the index file
#[derive(Serialize, Deserialize, Debug)]
pub struct IndexHeader {
    pub format_version: u8,
    pub kmer_length: u8,
    pub window_size: u8,
}

impl IndexHeader {
    pub fn new(kmer_length: u8, window_size: u8) -> Self {
        IndexHeader {
            format_version: INDEX_FORMAT_VERSION,
            kmer_length,
            window_size,
        }
    }

    /// Validate header
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.format_version != INDEX_FORMAT_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported index format version: {}",
                self.format_version
            ));
        }

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

    /// Get k
    pub fn kmer_length(&self) -> u8 {
        self.kmer_length
    }

    /// Get w
    pub fn window_size(&self) -> u8 {
        self.window_size
    }
}

/// Load just the header and count from an index file
pub fn load_header_and_count<P: AsRef<Path>>(path: &P) -> Result<(IndexHeader, usize)> {
    let file =
        File::open(path).context(format!("Failed to open index file {:?}", path.as_ref()))?;
    let mut reader = BufReader::new(file);

    let config = bincode::config::standard().with_fixed_int_encoding();

    // Deserialise header
    let header: IndexHeader =
        decode_from_std_read(&mut reader, config).context("Failed to deserialise index header")?;
    header.validate()?;

    // Deserialise the count of minimizers
    let count: usize = decode_from_std_read(&mut reader, config)
        .context("Failed to deserialise minimizer count")?;

    Ok((header, count))
}

static INDEX: OnceLock<(PathBuf, crate::MinimizerSet, IndexHeader)> = OnceLock::new();

pub fn load_minimizers_cached(
    path: &Path,
) -> Result<(&'static crate::MinimizerSet, &'static IndexHeader)> {
    let (p, minimizers, header) = INDEX.get_or_init(|| {
        let (m, h) = load_minimizers(path).unwrap();
        (path.to_owned(), m, h)
    });
    assert_eq!(
        p, path,
        "Currently, the server can only have one index loaded."
    );

    Ok((minimizers, header))
}

/// Load minimizers from an index file
pub fn load_minimizers(path: &Path) -> Result<(crate::MinimizerSet, IndexHeader)> {
    let file = File::open(path).context(format!("Failed to open index file {:?}", path))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let config = bincode::config::standard().with_fixed_int_encoding();

    // Deserialise header
    let header: IndexHeader =
        decode_from_std_read(&mut reader, config).context("Failed to deserialise index header")?;
    header.validate()?;

    // Deserialise the count of minimizers
    let count: usize = decode_from_std_read(&mut reader, config)
        .context("Failed to deserialise minimizer count")?;

    let bytes_per_minimizer = (header.kmer_length as usize).div_ceil(4);

    let minimizers = if header.kmer_length <= 32 {
        // Read as u64 with packed byte-aligned format
        let mut set = RapidHashSet::<u64>::with_capacity_and_hasher(count, FixedRapidHasher);
        const B: usize = 16 * 1024;
        let mut buffer = vec![0u8; bytes_per_minimizer * B];
        for i in (0..count).step_by(B) {
            let batch_count = B.min(count - i);
            let batch_bytes = bytes_per_minimizer * batch_count;
            let batch = &mut buffer[..batch_bytes];
            reader.read_exact(batch).unwrap();

            // Extract minimizers from packed bytes
            for j in 0..batch_count {
                let start = j * bytes_per_minimizer;
                let end = start + bytes_per_minimizer;
                let mut minimizer_bytes = [0u8; 8];
                minimizer_bytes[..bytes_per_minimizer].copy_from_slice(&batch[start..end]);
                set.insert(u64::from_le_bytes(minimizer_bytes));
            }
        }
        crate::MinimizerSet::U64(set)
    } else {
        // Read as u128 with packed byte-aligned format
        let mut set = RapidHashSet::<u128>::with_capacity_and_hasher(count, FixedRapidHasher);
        const B: usize = 16 * 1024;
        let mut buffer = vec![0u8; bytes_per_minimizer * B];
        for i in (0..count).step_by(B) {
            let batch_count = B.min(count - i);
            let batch_bytes = bytes_per_minimizer * batch_count;
            let batch = &mut buffer[..batch_bytes];
            reader.read_exact(batch).unwrap();

            // Extract minimizers from packed bytes
            for j in 0..batch_count {
                let start = j * bytes_per_minimizer;
                let end = start + bytes_per_minimizer;
                let mut minimizer_bytes = [0u8; 16];
                minimizer_bytes[..bytes_per_minimizer].copy_from_slice(&batch[start..end]);
                set.insert(u128::from_le_bytes(minimizer_bytes));
            }
        }
        crate::MinimizerSet::U128(set)
    };

    Ok((minimizers, header))
}

/// Helper function to write minimizers to output file or stdout
pub fn dump_minimizers(
    minimizers: &crate::MinimizerSet,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    // Create writer based on output path
    let mut writer: BufWriter<Box<dyn Write>> = match output_path {
        Some(path) if path.as_os_str() != "-" => BufWriter::new(Box::new(
            File::create(path).context("Failed to create output file")?,
        )),
        _ => BufWriter::new(Box::new(io::stdout())),
    };

    // Use fixed-width little-endian encoding for all values.
    let config = bincode::config::standard().with_fixed_int_encoding();
    encode_into_std_write(header, &mut writer, config)
        .context("Failed to serialise index header")?;

    // Serialise the count of minimizers first
    let count = minimizers.len();
    encode_into_std_write(count, &mut writer, config)
        .context("Failed to serialise minimizer count")?;

    // Serialise minimizers in byte-aligned packed format
    let bytes_per_minimizer = (header.kmer_length as usize).div_ceil(4);
    match minimizers {
        crate::MinimizerSet::U64(set) => {
            for &val in set {
                // Write only the required bytes (little-endian)
                let bytes = val.to_le_bytes();
                writer
                    .write_all(&bytes[..bytes_per_minimizer])
                    .context("Failed to write minimizer")?;
            }
        }
        crate::MinimizerSet::U128(set) => {
            for &val in set {
                // Write only the required bytes (little-endian)
                let bytes = val.to_le_bytes();
                writer
                    .write_all(&bytes[..bytes_per_minimizer])
                    .context("Failed to write minimizer")?;
            }
        }
    }
    Ok(())
}

/// Dump indexed minimizers to FASTA
pub fn dump(index_path: &Path, output_path: Option<&Path>) -> Result<()> {
    // Load the index
    let (minimizers, header) = load_minimizers(index_path).context("Failed to load index")?;

    // Create writer based on output path
    let mut writer: BufWriter<Box<dyn Write>> = match output_path {
        Some(path) if path.as_os_str() != "-" => BufWriter::new(Box::new(
            File::create(path).context("Failed to create output file")?,
        )),
        _ => BufWriter::new(Box::new(io::stdout())),
    };

    // Write FASTA
    let mut counter = 0;
    match minimizers {
        crate::MinimizerSet::U64(set) => {
            for &minimizer in &set {
                counter += 1;
                let sequence = crate::minimizers::decode_u64(minimizer, header.kmer_length);
                writeln!(writer, ">{}", counter)?;
                writeln!(writer, "{}", String::from_utf8_lossy(&sequence))?;
            }
        }
        crate::MinimizerSet::U128(set) => {
            for &minimizer in &set {
                counter += 1;
                let sequence = crate::minimizers::decode_u128(minimizer, header.kmer_length);
                writeln!(writer, ">{}", counter)?;
                writeln!(writer, "{}", String::from_utf8_lossy(&sequence))?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}

fn reader_with_inferred_batch_size(
    in_path: Option<&Path>,
) -> Result<paraseq::fastx::Reader<Box<dyn Read + Send>>> {
    let mut reader = paraseq::fastx::Reader::from_optional_path(in_path).unwrap();
    reader.update_batch_size_in_bp(256 * 1024)?;
    Ok(reader)
}

#[derive(Clone)]
struct BuildIndexProcessor<'c> {
    config: &'c IndexConfig,
    hasher: KmerHasher,
    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_minimizers_u64: Option<RapidHashSet<u64>>,
    local_minimizers_u128: Option<RapidHashSet<u128>>,
    local_record_info: Vec<(Vec<u8>, usize)>, // (record_id, seq_len) for progress output
    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_minimizers_u64: Arc<Mutex<Option<RapidHashSet<u64>>>>,
    global_minimizers_u128: Arc<Mutex<Option<RapidHashSet<u128>>>>,
}

impl<Rf: Record> ParallelProcessor<Rf> for BuildIndexProcessor<'_> {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        crate::minimizers::fill_minimizers(
            &seq,
            &self.hasher,
            self.config.kmer_length,
            self.config.window_size,
            self.config.entropy_threshold,
            &mut self.buffers,
        );

        // Extend appropriate local set based on type
        match &mut self.buffers.minimizers {
            crate::MinimizerVec::U64(vec) => {
                self.local_minimizers_u64
                    .as_mut()
                    .unwrap()
                    .extend(vec.iter());
            }
            crate::MinimizerVec::U128(vec) => {
                self.local_minimizers_u128
                    .as_mut()
                    .unwrap()
                    .extend(vec.iter());
            }
        }

        // Store record info for progress output (printed sequentially in on_batch_complete)
        if !self.config.quiet {
            self.local_record_info
                .push((record.id().to_vec(), seq.len()));
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local minimizers into global set and get new total count
        let minimizer_count = if let Some(local) = &mut self.local_minimizers_u64 {
            let mut global = self.global_minimizers_u64.lock();
            global.as_mut().unwrap().extend(local.iter());
            local.clear();
            global.as_ref().unwrap().len()
        } else {
            let mut global = self.global_minimizers_u128.lock();
            global
                .as_mut()
                .unwrap()
                .extend(self.local_minimizers_u128.as_ref().unwrap().iter());
            self.local_minimizers_u128.as_mut().unwrap().clear();
            global.as_ref().unwrap().len()
        };

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;
            self.local_stats = ProcessingStats::default();
        }

        // Print progress for each record in this batch (sequentially, after merging)
        if !self.config.quiet {
            for (record_id, seq_len) in &self.local_record_info {
                let id_str = std::str::from_utf8(record_id).unwrap_or("unknown");
                eprintln!(
                    "  {} ({}bp), total minimizers: {}",
                    id_str, seq_len, minimizer_count
                );
            }
            self.local_record_info.clear();
        }

        Ok(())
    }
}

/// Build an index of minimizers from a fastx file
pub fn build(config: &IndexConfig) -> Result<()> {
    let start_time = Instant::now();
    let path = &config.input_path;

    let version: String = env!("CARGO_PKG_VERSION").to_string();

    // Build options string similar to filter
    let mut options = Vec::<String>::new();
    if config.threads > 0 {
        options.push(format!("threads={}", config.threads));
    }

    eprintln!(
        "Deacon v{}; mode: build; input: single; options: {}",
        version,
        options.join(", ")
    );

    // Validate k-mer and window size constraints
    config.validate()?;

    let in_path = if path.as_os_str() == "-" {
        None
    } else {
        Some(path.as_path())
    };
    let reader = reader_with_inferred_batch_size(in_path)?;

    eprintln!(
        "Building index (k={}, w={})",
        config.kmer_length, config.window_size
    );

    let mut processor = if config.kmer_length <= 32 {
        BuildIndexProcessor {
            config,
            hasher: KmerHasher::new(config.kmer_length as usize),
            local_stats: ProcessingStats::default(),
            buffers: Buffers::new_u64(),
            local_minimizers_u64: Some(RapidHashSet::default()),
            local_minimizers_u128: None,
            local_record_info: Vec::new(),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            global_minimizers_u64: Arc::new(Mutex::new(Some(RapidHashSet::default()))),
            global_minimizers_u128: Arc::new(Mutex::new(None)),
        }
    } else {
        BuildIndexProcessor {
            config,
            hasher: KmerHasher::new(config.kmer_length as usize),
            local_stats: ProcessingStats::default(),
            buffers: Buffers::new_u128(),
            local_minimizers_u64: None,
            local_minimizers_u128: Some(RapidHashSet::default()),
            local_record_info: Vec::new(),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            global_minimizers_u64: Arc::new(Mutex::new(None)),
            global_minimizers_u128: Arc::new(Mutex::new(Some(RapidHashSet::default()))),
        }
    };
    reader.process_parallel(&mut processor, config.threads as usize)?;

    let all_minimizers = if config.kmer_length <= 32 {
        let set = Arc::try_unwrap(processor.global_minimizers_u64)
            .unwrap()
            .into_inner()
            .unwrap();
        crate::MinimizerSet::U64(set)
    } else {
        let set = Arc::try_unwrap(processor.global_minimizers_u128)
            .unwrap()
            .into_inner()
            .unwrap();
        crate::MinimizerSet::U128(set)
    };
    let stats = Arc::try_unwrap(processor.global_stats)
        .unwrap()
        .into_inner();

    eprintln!(
        "Indexed {} minimizers from {} record(s) ({}bp)",
        all_minimizers.len(),
        stats.total_seqs,
        stats.total_bp
    );

    let header = IndexHeader::new(config.kmer_length, config.window_size);

    // Write to output path or stdout
    dump_minimizers(&all_minimizers, &header, config.output_path.as_deref())?;

    let total_time = start_time.elapsed();
    eprintln!("Completed in {:.2?}", total_time);

    Ok(())
}

#[derive(Clone)]
struct DiffIndexProcessor {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_minimizers_u64: Option<RapidHashSet<u64>>,
    local_minimizers_u128: Option<RapidHashSet<u128>>,
    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    initial_size: usize,
    global_minimizers_u64: Arc<Mutex<Option<RapidHashSet<u64>>>>,
    global_minimizers_u128: Arc<Mutex<Option<RapidHashSet<u128>>>>,
}

impl<Rf: Record> ParallelProcessor<Rf> for DiffIndexProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        crate::minimizers::fill_minimizers(
            &seq,
            &self.hasher,
            self.kmer_length,
            self.window_size,
            0.0,
            &mut self.buffers,
        );

        // Extend appropriate local set based on type
        match &mut self.buffers.minimizers {
            crate::MinimizerVec::U64(vec) => {
                self.local_minimizers_u64
                    .as_mut()
                    .unwrap()
                    .extend(vec.iter());
            }
            crate::MinimizerVec::U128(vec) => {
                self.local_minimizers_u128
                    .as_mut()
                    .unwrap()
                    .extend(vec.iter());
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Write buffer to output
        let len = if let Some(local) = &mut self.local_minimizers_u64 {
            let mut global = self.global_minimizers_u64.lock();
            let global_set = global.as_mut().unwrap();
            for &minimizer in local.iter() {
                global_set.remove(&minimizer);
            }
            let len = global_set.len();
            local.clear();
            len
        } else {
            let mut global = self.global_minimizers_u128.lock();
            let global_set = global.as_mut().unwrap();
            for &minimizer in self.local_minimizers_u128.as_ref().unwrap().iter() {
                global_set.remove(&minimizer);
            }
            let len = global_set.len();
            self.local_minimizers_u128.as_mut().unwrap().clear();
            len
        };

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            let current_gb = stats.total_bp / 1_000_000_000;
            if current_gb > stats.last_reported {
                eprintln!(
                    "  Processed {} sequences ({}bp), removed {} minimizers",
                    stats.total_seqs,
                    stats.total_bp,
                    self.initial_size - len
                );
                stats.last_reported = current_gb;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

/// Stream minimizers from a FASTX file or stdin and remove those present in first_minimizers
fn stream_diff_fastx(
    fastx_path: &Path,
    kmer_length: u8,
    window_size: u8,
    first_header: &IndexHeader,
    threads: u16,
    first_minimizers: &mut crate::MinimizerSet,
) -> Result<(usize, usize)> {
    let path = fastx_path;

    // Validate k-mer and window size constraints
    let temp_config = crate::IndexConfig {
        input_path: PathBuf::new(),
        kmer_length,
        window_size,
        output_path: None,
        threads: 0,
        quiet: false,
        entropy_threshold: 0.0,
    };
    temp_config.validate()?;

    // Validate parameters match the first index
    if kmer_length != first_header.kmer_length() || window_size != first_header.window_size() {
        return Err(anyhow::anyhow!(
            "FASTX parameters (k={}, w={}) must match first index (k={}, w={})",
            kmer_length,
            window_size,
            first_header.kmer_length(),
            first_header.window_size()
        ));
    }

    if path.as_os_str() == "-" {
        eprintln!(
            "Second index: processing FASTX from stdin (k={}, w={})…",
            kmer_length, window_size
        );
    } else {
        eprintln!(
            "Second index: processing FASTX from file (k={}, w={})…",
            kmer_length, window_size
        );
    }

    let in_path = if path.as_os_str() == "-" {
        None
    } else {
        Some(path)
    };

    let reader = reader_with_inferred_batch_size(in_path)?;

    // Move first_minimizers into Arc<Mutex<>> for parallel access
    let initial_size = first_minimizers.len();
    let (mut processor, global_minimizers_u64, global_minimizers_u128) = match first_minimizers {
        crate::MinimizerSet::U64(set) => {
            let moved_set = std::mem::take(set);
            let arc = Arc::new(Mutex::new(Some(moved_set)));
            (
                DiffIndexProcessor {
                    kmer_length,
                    window_size,
                    hasher: KmerHasher::new(kmer_length as usize),
                    local_stats: ProcessingStats::default(),
                    buffers: Buffers::new_u64(),
                    local_minimizers_u64: Some(RapidHashSet::default()),
                    local_minimizers_u128: None,
                    global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
                    initial_size,
                    global_minimizers_u64: arc.clone(),
                    global_minimizers_u128: Arc::new(Mutex::new(None)),
                },
                Some(arc),
                None,
            )
        }
        crate::MinimizerSet::U128(set) => {
            let moved_set = std::mem::take(set);
            let arc = Arc::new(Mutex::new(Some(moved_set)));
            (
                DiffIndexProcessor {
                    kmer_length,
                    window_size,
                    hasher: KmerHasher::new(kmer_length as usize),
                    local_stats: ProcessingStats::default(),
                    buffers: Buffers::new_u128(),
                    local_minimizers_u64: None,
                    local_minimizers_u128: Some(RapidHashSet::default()),
                    global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
                    initial_size,
                    global_minimizers_u64: Arc::new(Mutex::new(None)),
                    global_minimizers_u128: arc.clone(),
                },
                None,
                Some(arc),
            )
        }
    };

    reader.process_parallel(&mut processor, threads as usize)?;

    // Extract results from Arc<Mutex<>> after processing completes
    let stats = processor.global_stats.lock().clone();
    *first_minimizers = if let Some(arc) = global_minimizers_u64 {
        let set = arc.lock().take().unwrap();
        crate::MinimizerSet::U64(set)
    } else {
        let set = global_minimizers_u128.unwrap().lock().take().unwrap();
        crate::MinimizerSet::U128(set)
    };

    eprintln!(
        "Processed {} sequences ({}bp) from FASTX file",
        stats.total_seqs, stats.total_bp
    );

    Ok((stats.total_seqs as usize, stats.total_bp as usize))
}

/// Compute the set difference between two minimizer indexes (A - B)
pub fn diff(
    first: &Path,
    second: &Path,
    kmer_length: Option<u8>,
    window_size: Option<u8>,
    threads: u16,
    output: Option<&Path>,
) -> Result<()> {
    let start_time = Instant::now();

    // Load first file (always an index)
    let (mut first_minimizers, header) = load_minimizers(first)?;
    eprintln!("First index: loaded {} minimizers", first_minimizers.len());

    // Guess if second file is an index or FASTX file
    let second_minimizers = if let (Some(k), Some(w)) = (kmer_length, window_size) {
        // Second file is a FASTX file - stream diff with provided k, w
        let before_count = first_minimizers.len();
        let (_seq_count, _total_bp) =
            stream_diff_fastx(second, k, w, &header, threads, &mut first_minimizers)?;

        // Report results
        eprintln!(
            "Removed {} minimizers, {} remaining",
            before_count - first_minimizers.len(),
            first_minimizers.len()
        );

        dump_minimizers(&first_minimizers, &header, output)?;

        let total_time = start_time.elapsed();
        eprintln!("Completed difference operation in {:.2?}", total_time);

        return Ok(());
    } else {
        // Try to load as index file first
        if let Ok((second_minimizers, second_header)) = load_minimizers(&second) {
            // Second file is an index file
            eprintln!(
                "Second index: loaded {} minimizers",
                second_minimizers.len()
            );

            // Validate compatible headers
            if second_header.kmer_length() != header.kmer_length()
                || second_header.window_size() != header.window_size()
            {
                return Err(anyhow::anyhow!(
                    "Incompatible headers: second index has k={}, w={}, but first index has k={}, w={}",
                    second_header.kmer_length(),
                    second_header.window_size(),
                    header.kmer_length(),
                    header.window_size()
                ));
            }

            second_minimizers
        } else {
            // Second file is not a valid index, treat as FASTX file
            // Use k and w from first index header and do a streaming diff
            let k = header.kmer_length();
            let w = header.window_size();

            // Count minimizers before diff
            let before_count = first_minimizers.len();

            let (_seq_count, _total_bp) =
                stream_diff_fastx(&second, k, w, &header, threads, &mut first_minimizers)?;

            // Report results
            eprintln!(
                "Removed {} minimizers, {} remaining",
                before_count - first_minimizers.len(),
                first_minimizers.len()
            );

            dump_minimizers(&first_minimizers, &header, output)?;

            let total_time = start_time.elapsed();
            eprintln!("Completed difference operation in {:.2?}", total_time);

            return Ok(());
        }
    };

    // Handle straightforward index-to-index diffing
    // Count minimizers before diff
    let before_count = first_minimizers.len();

    // Remove all minimizers in second_minimizers from first_minimizers
    first_minimizers.remove_all(&second_minimizers);

    // Report results
    eprintln!(
        "Removed {} minimizers, {} remaining",
        before_count - first_minimizers.len(),
        first_minimizers.len()
    );

    dump_minimizers(&first_minimizers, &header, output)?;

    let total_time = start_time.elapsed();
    eprintln!("Completed diff operation in {:.2?}", total_time);

    Ok(())
}

/// Show info about an index
pub fn info(index_path: &Path) -> Result<()> {
    let start_time = Instant::now();

    // Load index file
    let (minimizers, header) = load_minimizers(index_path)?;

    // Show index info
    eprintln!("Index information:");
    eprintln!("  Format version: {}", header.format_version);
    eprintln!("  K-mer length (k): {}", header.kmer_length());
    eprintln!("  Window size (w): {}", header.window_size());
    eprintln!("  Distinct minimizer count: {}", minimizers.len());

    let total_time = start_time.elapsed();
    eprintln!("Loaded index info in {:.2?}", total_time);

    Ok(())
}

/// Combine minimizer indexes (set union)
pub fn union(inputs: &[PathBuf], output: Option<&Path>) -> Result<()> {
    let start_time = Instant::now();
    // Check input files
    if inputs.is_empty() {
        return Err(anyhow::anyhow!(
            "No input files provided for union operation"
        ));
    }

    // Read all headers first to determine total capacity needed
    let mut headers_and_counts = Vec::new();

    for path in inputs {
        let (header, count) = load_header_and_count(path)?;
        headers_and_counts.push((header, count));
    }

    // Get header from first file for output
    let header = &headers_and_counts[0].0;

    eprintln!(
        "Combining indexes (k={}, w={})",
        header.kmer_length(),
        header.window_size()
    );

    // Verify all headers are compatible
    for (i, (file_header, _)) in headers_and_counts.iter().enumerate() {
        if file_header.kmer_length() != header.kmer_length()
            || file_header.window_size() != header.window_size()
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index {} has k={}, w={}, but first index has k={}, w={}",
                i,
                file_header.kmer_length(),
                file_header.window_size(),
                header.kmer_length(),
                header.window_size()
            ));
        }
    }

    // Load first index to determine type (u64 vs u128)
    let (mut all_minimizers, _) = load_minimizers(&inputs[0])?;
    eprintln!("Index 1: loaded {} minimizers", all_minimizers.len());

    // Now load and merge remaining indexes
    for (i, path) in inputs.iter().enumerate().skip(1) {
        let (minimizers, _) = load_minimizers(path)?;
        let before_count = all_minimizers.len();

        // Merge minimizers (set union)
        all_minimizers.extend(minimizers);

        let expected_count = headers_and_counts[i].1;
        eprintln!(
            "Index {}: expected {} minimizers, added {}, total: {}",
            i + 1,
            expected_count,
            all_minimizers.len() - before_count,
            all_minimizers.len()
        );
    }

    dump_minimizers(&all_minimizers, header, output)?;

    let total_time = start_time.elapsed();
    eprintln!(
        "United {} indexes with {} total minimizers in {:.2?}",
        inputs.len(),
        all_minimizers.len(),
        total_time
    );

    Ok(())
}

pub fn intersect(inputs: &[PathBuf], output: Option<&Path>) -> Result<()> {
    let start_time = Instant::now();
    // Check inputs
    if inputs.len() < 2 {
        return Err(anyhow::anyhow!(
            "At least two input files are required for intersection operation"
        ));
    }

    // Read all headers first
    let mut headers_and_counts = Vec::new();

    for path in inputs {
        let (header, count) = load_header_and_count(path)?;
        headers_and_counts.push((header, count));
    }

    // Get header from first file for output
    let header = &headers_and_counts[0].0;

    eprintln!(
        "Intersecting indexes (k={}, w={})",
        header.kmer_length(),
        header.window_size()
    );

    // Check header compat
    for (i, (file_header, _)) in headers_and_counts.iter().enumerate() {
        if file_header.kmer_length() != header.kmer_length()
            || file_header.window_size() != header.window_size()
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index {} has k={}, w={}, but first index has k={}, w={}",
                i,
                file_header.kmer_length(),
                file_header.window_size(),
                header.kmer_length(),
                header.window_size()
            ));
        }
    }

    // Load first index
    let (mut result_minimizers, _) = load_minimizers(&inputs[0])?;
    eprintln!("Index 1: loaded {} minimizers", result_minimizers.len());

    // Intersect with remaining indexes
    for (i, path) in inputs.iter().enumerate().skip(1) {
        let (minimizers, _) = load_minimizers(path)?;

        // Intersect minimizers (set intersection)
        result_minimizers.intersect(&minimizers);

        let expected_count = headers_and_counts[i].1;
        eprintln!(
            "Index {}: expected {} minimizers, retained {}, total: {}",
            i + 1,
            expected_count,
            result_minimizers.len(),
            result_minimizers.len()
        );
    }

    dump_minimizers(&result_minimizers, header, output)?;

    let total_time = start_time.elapsed();
    eprintln!(
        "Intersected {} indexes with {} common minimizers in {:.2?}",
        inputs.len(),
        result_minimizers.len(),
        total_time
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_creation() {
        let header = IndexHeader::new(31, 21);

        assert_eq!(header.format_version, 3);
        assert_eq!(header.kmer_length(), 31);
        assert_eq!(header.window_size(), 21);
    }

    #[test]
    fn test_header_validation() {
        // Valid header (v3 only)
        let valid_header_v3 = IndexHeader {
            format_version: 3,
            kmer_length: 31,
            window_size: 21,
        };
        assert!(valid_header_v3.validate().is_ok());

        // Invalid format versions
        let invalid_header_v1 = IndexHeader {
            format_version: 1,
            kmer_length: 31,
            window_size: 21,
        };
        assert!(invalid_header_v1.validate().is_err());

        let invalid_header_v2 = IndexHeader {
            format_version: 2,
            kmer_length: 31,
            window_size: 21,
        };
        assert!(invalid_header_v2.validate().is_err());

        let invalid_header_v4 = IndexHeader {
            format_version: 4,
            kmer_length: 31,
            window_size: 21,
        };
        assert!(invalid_header_v4.validate().is_err());

        // Test k <= 61 constraint
        let valid_max_k = IndexHeader {
            format_version: 3,
            kmer_length: 61,
            window_size: 15,
        };
        assert!(valid_max_k.validate().is_ok()); // k=61 is valid

        let invalid_max_k = IndexHeader {
            format_version: 3,
            kmer_length: 63,
            window_size: 15,
        };
        assert!(invalid_max_k.validate().is_err()); // k=63 exceeds 61

        // Test k+w <= 96 constraint
        let valid_max_kw = IndexHeader {
            format_version: 3,
            kmer_length: 61,
            window_size: 35,
        };
        assert!(valid_max_kw.validate().is_ok()); // k+w=96 is valid

        let invalid_max_kw = IndexHeader {
            format_version: 3,
            kmer_length: 61,
            window_size: 36,
        };
        assert!(invalid_max_kw.validate().is_err()); // k+w=97 exceeds 96

        // Test k must be odd
        let invalid_even = IndexHeader {
            format_version: 3,
            kmer_length: 30,
            window_size: 15,
        };
        assert!(invalid_even.validate().is_err()); // k=30 is even
    }
}

/// Fetch a pre-built index from remote storage
#[cfg(feature = "fetch")]
pub fn fetch(
    index_name: &str,
    kmer_length: u8,
    window_size: u8,
    output: Option<&Path>,
) -> Result<()> {
    const DEFAULT_REPOSITORY_URL: &str =
        "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o";

    let base_url = std::env::var("DEACON_REPOSITORY_URL")
        .unwrap_or_else(|_| DEFAULT_REPOSITORY_URL.to_string());

    let filename = format!("{}.k{}w{}.idx", index_name, kmer_length, window_size);
    let url = format!("{}/deacon/{}/{}", base_url, INDEX_FORMAT_VERSION, filename);

    eprintln!("Fetching {}", url);

    let agent = ureq::Agent::config_builder()
        .timeout_connect(Some(Duration::from_secs(15)))
        .build()
        .new_agent();

    let response = agent.get(&url).call().context("Failed to download index")?;

    if response.status() != 200 {
        anyhow::bail!("Failed to fetch index: HTTP {}", response.status());
    }

    let content_length = response
        .headers()
        .get("Content-Length")
        .and_then(|v| v.to_str().ok())
        .and_then(|s| s.parse::<u64>().ok());

    let pb = ProgressBar::new(content_length.unwrap_or(0));

    let mut body = response.into_body();
    let output_path = output
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from(&filename));

    let mut temp_path = output_path.clone();
    temp_path.as_mut_os_string().push(".tmp");

    let mut file = File::create(&temp_path).context("Failed to create temporary file")?;
    std::io::copy(&mut pb.wrap_read(body.as_reader()), &mut file)
        .context("Failed to write index to file")?;

    std::fs::rename(&temp_path, &output_path).context("Failed to finalise index")?;

    pb.finish_and_clear();
    eprintln!("Index saved to: {}", output_path.display());

    Ok(())
}
