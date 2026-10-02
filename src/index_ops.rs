//! Index building, set operations, inspection and conversion.

use crate::dedup_vec::DedupVec;
use crate::filter_io::{is_special_input_path, open_cbq};
use crate::index::{FixedRapidHasher, IndexStorage, RapidHashSet};
use crate::index_format::{
    FuseIndexHeader, INDEX_FORMAT_VERSION, IndexHeader, load_minimizers_from_path, table_buckets,
    write_exact_index_shards_to_path, write_exact_index_to_path,
};
use crate::minimizers::{MinimizerShards, Minimizers, validate_k_w};
use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use binseq::{BinseqRecord, ParallelReader as BinSeqParalleReader};
use paraseq::Record;
use paraseq::prelude::{ParallelProcessor, ParallelReader};
use parking_lot::Mutex;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use std::fs::File;
use std::hash::BuildHasher;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use tracing::info;

/// Load an index header and minimizer count.
pub fn load_header_and_count<P: AsRef<Path>>(path: &P) -> Result<(IndexHeader, usize)> {
    let file = std::fs::File::open(path)
        .context(format!("Failed to open index file {:?}", path.as_ref()))?;
    let mut reader = BufReader::new(file);

    let config = bincode::config::standard().with_fixed_int_encoding();

    // Deserialise header
    let header: IndexHeader =
        decode_from_std_read(&mut reader, config).context("Failed to deserialise index header")?;
    header.validate()?;

    // Deserialise the count of minimizers (stored as u64 for cross-platform compatibility)
    let count: u64 = decode_from_std_read(&mut reader, config)
        .context("Failed to deserialise minimizer count")?;

    Ok((header, count as usize))
}

/// Reshard by high bucket bits, then sort shards in parallel for concatenation.
fn sort_sharded_lists<T>(shards: Vec<Vec<T>>) -> Vec<Vec<T>>
where
    T: Copy + std::hash::Hash + Ord + Send + Sync,
{
    // Match the bucket count of the single, final hash set used by
    // sort_hashset so concatenating these vectors preserves its ordering.
    let total_len: usize = shards.iter().map(Vec::len).sum();
    let num_buckets = table_buckets(total_len);
    let bucket = |x: &T| -> usize { FixedRapidHasher.hash_one(x) as usize & (num_buckets - 1) };
    assert!(num_buckets.is_power_of_two());
    assert!(SHARDS.is_power_of_two());
    // Assign contiguous ranges of target buckets to each output shard.
    // Tables with fewer buckets than shards use one shard per bucket.
    let shift = num_buckets
        .trailing_zeros()
        .saturating_sub(SHARDS.trailing_zeros());

    info!("Resharding minimizers");
    let sorted_shards: Vec<Mutex<Vec<T>>> = (0..SHARDS).map(|_| Mutex::new(Vec::new())).collect();
    shards.into_par_iter().for_each(|shard| {
        debug_assert!(shard.array_windows().all(|[v1, v2]| v1 < v2));
        let mut buffers = (0..SHARDS).map(|_| vec![]).collect::<Vec<_>>();
        for value in shard {
            let target_shard = bucket(&value) >> shift;
            buffers[target_shard].push(value);
        }
        for shard_idx in random_shard_order() {
            sorted_shards[shard_idx]
                .lock()
                .extend(buffers[shard_idx].drain(..));
        }
    });

    info!("Sorting shards");
    let mut sorted_shards = sorted_shards
        .into_iter()
        .map(|mutex| mutex.into_inner())
        .collect::<Vec<_>>();

    sorted_shards.par_iter_mut().for_each(|values| {
        values.sort_unstable_by_key(|value| (bucket(value), *value));
    });

    for [s1, s2] in sorted_shards.array_windows() {
        assert!(s1.last().is_none_or(|v1| {
            s2.first()
                .is_none_or(|v2| (bucket(v1), v1) <= (bucket(v2), v2))
        }));
    }

    sorted_shards
}

/// Read file magic, skipping pipes to avoid consuming input.
fn sniff(path: &Path, buf: &mut [u8]) -> bool {
    !is_special_input_path(path) && File::open(path).and_then(|mut f| f.read_exact(buf)).is_ok()
}

/// Detect BFF magic.
fn is_fuse_index_file(path: &Path) -> bool {
    let mut magic = [0u8; 4];
    sniff(path, &mut magic) && magic.starts_with(b"DBF")
}

/// Exact indexes start with a version byte.
fn is_exact_index_file(path: &Path) -> bool {
    let mut first = [0u8; 1];
    sniff(path, &mut first) && first[0] <= INDEX_FORMAT_VERSION
}

/// Reject BFF indexes for exact-only operations.
fn reject_fuse_index(path: &Path, operation: &str) -> Result<()> {
    if is_fuse_index_file(path) {
        return Err(anyhow::anyhow!(
            "{:?} is a binary fuse filter (.pidx) index; {} requires an exact (.idx) index",
            path,
            operation
        ));
    }
    Ok(())
}

/// Dump minimizers to FASTA.
pub fn dump(index_path: &Path, output_path: Option<&Path>) -> Result<()> {
    if is_fuse_index_file(index_path) {
        return Err(anyhow::anyhow!(
            "Cannot dump a BFF index: keys are not recoverable from a binary fuse filter"
        ));
    }

    // Load the index
    let (minimizers, header) =
        load_minimizers_from_path(index_path).context("Failed to load index")?;

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
        IndexStorage::ExactU64(set) => {
            for &minimizer in &set {
                counter += 1;
                let sequence = crate::minimizers::decode_u64(minimizer, header.kmer_length);
                writeln!(writer, ">{}", counter)?;
                writeln!(writer, "{}", String::from_utf8_lossy(&sequence))?;
            }
        }
        IndexStorage::ExactU128(set) => {
            for &minimizer in &set {
                counter += 1;
                let sequence = crate::minimizers::decode_u128(minimizer, header.kmer_length);
                writeln!(writer, ">{}", counter)?;
                writeln!(writer, "{}", String::from_utf8_lossy(&sequence))?;
            }
        }
        IndexStorage::Fuse(_) => {
            unreachable!("BFF dump is rejected by is_fuse_index_file above")
        }
    }

    writer.flush()?;
    Ok(())
}

/// Convert an exact index to BFF (k <= 32).
pub fn freeze(index_path: &Path, output_path: Option<&Path>, bits: u8) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    info!("Deacon v{}; mode: freeze", version);

    if !matches!(bits, 16 | 32) {
        return Err(anyhow::anyhow!(
            "Unsupported BFF fingerprint width: {} bits (expected 16 or 32)",
            bits
        ));
    }

    reject_fuse_index(index_path, "freeze")?;

    let (minimizers, header) =
        load_minimizers_from_path(index_path).context("Failed to load index")?;

    // Collect unique u64 keys
    let keys: Vec<u64> = match minimizers {
        IndexStorage::ExactU64(set) => set.into_iter().collect(),
        IndexStorage::ExactU128(_) => {
            return Err(anyhow::anyhow!(
                "BFF supports k <= 32 (u64 minimizers); got k={} (u128). Build a k<=32 index first.",
                header.kmer_length
            ));
        }
        IndexStorage::Fuse(_) => {
            return Err(anyhow::anyhow!("Input is already a BFF index"));
        }
    };

    let key_count = keys.len();
    info!(
        "Building {}-bit binary fuse filter from {} minimizers (k={}, w={})",
        bits, key_count, header.kmer_length, header.window_size
    );

    let fuse_header = FuseIndexHeader::new(
        bits,
        header.kmer_length,
        header.window_size,
        key_count as u64,
    );

    // Header via serde, filter via native bincode
    let mut writer: BufWriter<Box<dyn Write>> = match output_path {
        Some(path) if path.as_os_str() != "-" => BufWriter::new(Box::new(
            File::create(path).context("Failed to create output file")?,
        )),
        _ => BufWriter::new(Box::new(io::stdout())),
    };
    let config = bincode::config::standard().with_fixed_int_encoding();
    encode_into_std_write(&fuse_header, &mut writer, config)
        .context("Failed to serialise BFF header")?;

    // Construct and serialise the filter at the requested fingerprint width
    let filter_bytes = match bits {
        16 => {
            let filter = xorf::BinaryFuse16::try_from(&keys)
                .map_err(|e| anyhow::anyhow!("Failed to construct binary fuse filter: {}", e))?;
            bincode::encode_into_std_write(&filter, &mut writer, config)
                .context("Failed to serialise BFF filter")?
        }
        _ => {
            let filter = xorf::BinaryFuse32::try_from(&keys)
                .map_err(|e| anyhow::anyhow!("Failed to construct binary fuse filter: {}", e))?;
            bincode::encode_into_std_write(&filter, &mut writer, config)
                .context("Failed to serialise BFF filter")?
        }
    };
    writer.flush()?;

    let bits_per_key = if key_count > 0 {
        (filter_bytes as f64 * 8.0) / key_count as f64
    } else {
        0.0
    };
    info!(
        "Wrote BFF: {} keys, {} filter bytes ({:.2} bits/key); false-positive rate ~2^-{}",
        key_count, filter_bytes, bits_per_key, bits
    );
    info!("Completed freeze in {:.2?}", start_time.elapsed());
    Ok(())
}

fn reader_with_inferred_batch_size(
    in_path: Option<&Path>,
) -> Result<paraseq::fastx::Reader<Box<dyn Read + Send>>> {
    let mut reader = paraseq::ReaderBuilder::optional_path(in_path)
        .build()
        .map_err(|e| anyhow::anyhow!("Failed to open input: {}", e))?;
    reader.update_batch_size_in_bp(256 * 1024)?;
    Ok(reader)
}

#[derive(Clone, Default)]
struct IndexStats {
    total_seqs: u64,
    total_bp: u64,
    last_reported: u64,
}

#[derive(Clone)]
pub struct BuildConfig {
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
}

impl BuildConfig {
    /// Default k/w, automatic threads and stdout output. `-` reads stdin.
    pub fn new(input_path: impl Into<PathBuf>) -> Self {
        Self {
            input_path: input_path.into(),
            kmer_length: crate::DEFAULT_KMER_LENGTH,
            window_size: crate::DEFAULT_WINDOW_SIZE,
            output_path: None,
            threads: 0,
        }
    }

    /// Validate k-mer and window size constraints
    pub fn validate(&self) -> Result<()> {
        validate_k_w(self.kmer_length, self.window_size)
    }
}

#[derive(Clone)]
struct BuildIndexProcessor {
    /// Paired CBQ input, mates indexed as separate records
    paired: bool,
    // Local buffers
    seq: Vec<u8>,
    minimizers: Minimizers,
    local_stats: IndexStats,
    local_minimizers_u64: Option<Vec<Vec<u64>>>,
    local_minimizers_u128: Option<Vec<Vec<u128>>>,
    // Global state
    global_stats: Arc<Mutex<IndexStats>>,
    global_minimizers_u64: Arc<Vec<Mutex<DedupVec<u64>>>>,
    global_minimizers_u128: Arc<Vec<Mutex<DedupVec<u128>>>>,
}

const SHARDS: usize = 1024;

const LOCAL_BUF_SIZE: usize = 1024;

impl BuildIndexProcessor {
    /// Count and index one record
    fn add_seq(&mut self, seq: &[u8]) {
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        // Partition minimizers by value so each worker can merge its shards
        // independently at thread completion.
        match self.minimizers.compute(seq) {
            crate::MinimizerVec::U64(vec) => {
                let local = self.local_minimizers_u64.as_mut().unwrap();
                for &minimizer in vec.iter() {
                    let shard = (minimizer % SHARDS as u64) as usize;
                    local[shard].push(minimizer);

                    if local[shard].len() >= LOCAL_BUF_SIZE {
                        let mut global = self.global_minimizers_u64[shard].lock();
                        global.extend(&mut local[shard]);
                    }
                }
            }
            crate::MinimizerVec::U128(vec) => {
                let local = self.local_minimizers_u128.as_mut().unwrap();
                for &minimizer in vec.iter() {
                    let shard = (minimizer % SHARDS as u128) as usize;
                    local[shard].push(minimizer);

                    if local[shard].len() >= LOCAL_BUF_SIZE {
                        let mut global = self.global_minimizers_u128[shard].lock();
                        global.extend(&mut local[shard]);
                    }
                }
            }
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for BuildIndexProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        self.add_seq(&record.seq());
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        // Tick to stderr once every Gbp
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            let current_gb = stats.total_bp / 1_000_000_000;
            if current_gb > stats.last_reported {
                info!(
                    "Processed {} sequences ({}bp)",
                    stats.total_seqs, stats.total_bp
                );
                stats.last_reported = current_gb;
            }

            self.local_stats = IndexStats::default();
        }

        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        if let Some(local) = &mut self.local_minimizers_u64 {
            for shard in random_shard_order() {
                self.global_minimizers_u64[shard]
                    .lock()
                    .extend(&mut local[shard]);
            }
        } else {
            let local = self.local_minimizers_u128.as_mut().unwrap();
            for shard in random_shard_order() {
                self.global_minimizers_u128[shard]
                    .lock()
                    .extend(&mut local[shard]);
            }
        }
        Ok(())
    }
}

impl binseq::ParallelProcessor for BuildIndexProcessor {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        // Take the buffer so add_seq can borrow self
        let mut seq = std::mem::take(&mut self.seq);
        seq.clear();
        record.decode_s(&mut seq)?;
        self.add_seq(&seq);
        if self.paired {
            seq.clear();
            record.decode_x(&mut seq)?;
            self.add_seq(&seq);
        }
        self.seq = seq;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        // Tick to stderr once every Gbp
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            let current_gb = stats.total_bp / 1_000_000_000;
            if current_gb > stats.last_reported {
                info!(
                    "Processed {} sequences ({}bp)",
                    stats.total_seqs, stats.total_bp
                );
                stats.last_reported = current_gb;
            }

            self.local_stats = IndexStats::default();
        }

        Ok(())
    }

    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        if let Some(local) = &mut self.local_minimizers_u64 {
            for shard in random_shard_order() {
                self.global_minimizers_u64[shard]
                    .lock()
                    .extend(&mut local[shard]);
            }
        } else {
            let local = self.local_minimizers_u128.as_mut().unwrap();
            for shard in random_shard_order() {
                self.global_minimizers_u128[shard]
                    .lock()
                    .extend(&mut local[shard]);
            }
        }
        Ok(())
    }
}

/// Randomize shard order to reduce contention.
fn random_shard_order() -> Vec<usize> {
    let mut shard_order: Vec<_> = (0..SHARDS).collect();
    shard_order.shuffle(&mut rand::rng());
    shard_order
}

/// Build an index from FASTX or CBQ.
pub fn build(config: &BuildConfig) -> Result<()> {
    let start_time = Instant::now();
    let path = &config.input_path;

    let version: String = env!("CARGO_PKG_VERSION").to_string();

    // Build options string similar to filter
    let mut options = Vec::<String>::new();
    if config.threads > 0 {
        options.push(format!("threads={}", config.threads));
    }

    info!(
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

    info!(
        "Building index (k={}, w={})",
        config.kmer_length, config.window_size
    );

    let global_stats = Mutex::new(IndexStats::default());
    let global_minimizers_u64: Vec<Mutex<DedupVec<u64>>> = (0..SHARDS)
        .map(|_| Mutex::new(Default::default()))
        .collect();
    let global_minimizers_u128: Vec<Mutex<DedupVec<u128>>> = (0..SHARDS)
        .map(|_| Mutex::new(Default::default()))
        .collect();

    // Staging buffers matching minimizer width.
    fn local_shards<T>(used: bool) -> Option<Vec<Vec<T>>> {
        used.then(|| {
            (0..SHARDS)
                .map(|_| Vec::with_capacity(LOCAL_BUF_SIZE))
                .collect()
        })
    }
    let wide = config.kmer_length > 32;
    let mut processor = BuildIndexProcessor {
        minimizers: Minimizers::new(config.kmer_length, config.window_size)?,
        local_stats: IndexStats::default(),
        paired: false,
        seq: vec![],
        local_minimizers_u64: local_shards(!wide),
        local_minimizers_u128: local_shards(wide),
        global_stats: Arc::new(global_stats),
        global_minimizers_u64: Arc::new(global_minimizers_u64),
        global_minimizers_u128: Arc::new(global_minimizers_u128),
    };
    if let Some(reader) = open_cbq(path)? {
        processor.paired = reader.is_paired();
        // binseq's parallel reader rejects an empty record range
        if reader.num_records() > 0 {
            reader.process_parallel(processor.clone(), config.threads as usize)?;
        }
    } else {
        let reader = reader_with_inferred_batch_size(in_path)?;
        reader.process_parallel(&mut processor, config.threads as usize)?;
    }
    let global_stats = Arc::into_inner(processor.global_stats).expect("workers still hold stats");
    let global_minimizers_u64 =
        Arc::into_inner(processor.global_minimizers_u64).expect("workers still hold shards");
    let global_minimizers_u128 =
        Arc::into_inner(processor.global_minimizers_u128).expect("workers still hold shards");

    info!("Deduplicating shards");
    let all_minimizers = if config.kmer_length <= 32 {
        let shards: Vec<_> = global_minimizers_u64
            .into_par_iter()
            .map(|mutex| mutex.into_inner().finish())
            .collect();
        MinimizerShards::U64(sort_sharded_lists(shards))
    } else {
        let shards: Vec<_> = global_minimizers_u128
            .into_par_iter()
            .map(|mutex| mutex.into_inner().finish())
            .collect();
        MinimizerShards::U128(sort_sharded_lists(shards))
    };

    let stats = global_stats.into_inner();
    info!(
        "Indexed {} minimizers from {} record(s) ({}bp)",
        all_minimizers.len(),
        stats.total_seqs,
        stats.total_bp
    );

    let header = IndexHeader::new(config.kmer_length, config.window_size);

    // Write to output path or stdout
    write_exact_index_shards_to_path(all_minimizers, &header, config.output_path.as_deref())?;

    let total_time = start_time.elapsed();
    info!("Completed build in {:.2?}", total_time);

    Ok(())
}

/// Hits to subtract from the first index.
#[derive(Clone)]
enum HitSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl HitSet {
    /// An empty hit set matching the width of `set`
    fn empty_like(set: &IndexStorage) -> Self {
        match set {
            IndexStorage::ExactU64(_) => HitSet::U64(RapidHashSet::default()),
            IndexStorage::ExactU128(_) => HitSet::U128(RapidHashSet::default()),
            // load_minimizers rejects BFF files, so this is unreachable
            IndexStorage::Fuse(_) => {
                unreachable!("diff does not operate on BFF indexes")
            }
        }
    }

    fn len(&self) -> usize {
        match self {
            HitSet::U64(set) => set.len(),
            HitSet::U128(set) => set.len(),
        }
    }

    /// Drain `other` into `self`, keeping `other`'s allocation for reuse
    fn drain_from(&mut self, other: &mut Self) {
        match (self, other) {
            (HitSet::U64(dst), HitSet::U64(src)) => dst.extend(src.drain()),
            (HitSet::U128(dst), HitSet::U128(src)) => dst.extend(src.drain()),
            _ => unreachable!("minimizer width mismatch between hit sets"),
        }
    }

    /// Remove every hit from `set`
    fn remove_from(&self, set: &mut IndexStorage) {
        match (self, set) {
            (HitSet::U64(hits), IndexStorage::ExactU64(set)) => {
                for minimizer in hits {
                    set.remove(minimizer);
                }
            }
            (HitSet::U128(hits), IndexStorage::ExactU128(set)) => {
                for minimizer in hits {
                    set.remove(minimizer);
                }
            }
            _ => unreachable!("minimizer width mismatch between hit set and index"),
        }
    }
}

#[derive(Clone)]
struct DiffIndexProcessor<'a> {
    /// Index being subtracted from. Read-only while streaming, so every thread probes
    /// it without locking. Removal takes one pass after streaming.
    first: &'a IndexStorage,
    // Local buffers
    minimizers: Minimizers,
    local_stats: IndexStats,
    /// Hits seen by this thread since the last batch. Bounded by `first`, and usually a
    /// tiny fraction of it.
    local_hits: HitSet,
    // Global state
    global_stats: Arc<Mutex<IndexStats>>,
    global_hits: Arc<Mutex<HitSet>>,
}

impl<Rf: Record> ParallelProcessor<Rf> for DiffIndexProcessor<'_> {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        // Dispatch on width once, then probe lock-free, keeping only the hits
        match (
            self.first,
            self.minimizers.compute(&seq),
            &mut self.local_hits,
        ) {
            (IndexStorage::ExactU64(first), crate::MinimizerVec::U64(vec), HitSet::U64(hits)) => {
                for &minimizer in vec.iter() {
                    if first.contains(&minimizer) {
                        hits.insert(minimizer);
                    }
                }
            }
            (
                IndexStorage::ExactU128(first),
                crate::MinimizerVec::U128(vec),
                HitSet::U128(hits),
            ) => {
                for &minimizer in vec.iter() {
                    if first.contains(&minimizer) {
                        hits.insert(minimizer);
                    }
                }
            }
            _ => unreachable!("minimizer width mismatch between index and sequences"),
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        // Merge this batch's hits into the global set, holding the lock only for the
        // hits rather than for every minimizer seen
        let hits = {
            let mut global = self.global_hits.lock();
            global.drain_from(&mut self.local_hits);
            global.len()
        };

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            let current_10gb = stats.total_bp / 10_000_000_000;
            if current_10gb > stats.last_reported {
                info!(
                    "Processed {} sequences ({}bp), removed {} minimizers",
                    stats.total_seqs, stats.total_bp, hits
                );
                stats.last_reported = current_10gb;
            }

            self.local_stats = IndexStats::default();
        }

        Ok(())
    }
}

/// Subtract streamed FASTX minimizers from the first index.
fn stream_diff_fastx(
    fastx_path: &Path,
    window_size: u8,
    first_header: &IndexHeader,
    threads: u16,
    first_minimizers: &mut IndexStorage,
) -> Result<(usize, usize)> {
    let path = fastx_path;
    let kmer_length = first_header.kmer_length;
    let minimizers = Minimizers::new(kmer_length, window_size)?;

    // w must match or be 1 (w=1 emits every k-mer, for exact masking)
    if window_size != first_header.window_size && window_size != 1 {
        return Err(anyhow::anyhow!(
            "FASTX w={} must match first index w={} or be 1 (for exact k-mer subtraction)",
            window_size,
            first_header.window_size
        ));
    }

    if path.as_os_str() == "-" {
        info!(
            "Second index: processing FASTX from stdin (k={}, w={})",
            kmer_length, window_size
        );
    } else {
        info!(
            "Second index: processing FASTX from file (k={}, w={})",
            kmer_length, window_size
        );
    }

    let in_path = if path.as_os_str() == "-" {
        None
    } else {
        Some(path)
    };

    let reader = reader_with_inferred_batch_size(in_path)?;

    // Read only index while streaming - remove minimizers in single final op
    let global_stats = Arc::new(Mutex::new(IndexStats::default()));
    let global_hits = Arc::new(Mutex::new(HitSet::empty_like(first_minimizers)));
    let start_time = Instant::now();

    let mut processor = DiffIndexProcessor {
        first: first_minimizers,
        minimizers,
        local_stats: IndexStats::default(),
        local_hits: HitSet::empty_like(first_minimizers),
        global_stats: global_stats.clone(),
        global_hits: global_hits.clone(),
    };

    reader.process_parallel(&mut processor, threads as usize)?;
    drop(processor);

    let stats = global_stats.lock().clone();
    global_hits.lock().remove_from(first_minimizers);

    let elapsed = start_time.elapsed();
    info!(
        "Processed {} sequences ({}bp) from FASTX file in {:.2?} ({:.1} Mbp/s)",
        stats.total_seqs,
        stats.total_bp,
        elapsed,
        stats.total_bp as f64 / elapsed.as_secs_f64() / 1_000_000.0
    );

    Ok((stats.total_seqs as usize, stats.total_bp as usize))
}

/// Subtract an index or FASTX input (A - B).
pub fn diff(
    first: &Path,
    second: &Path,
    window_size: Option<u8>,
    threads: u16,
    output: Option<&Path>,
) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    info!("Deacon v{}; mode: diff", version);

    reject_fuse_index(first, "diff")?;
    reject_fuse_index(second, "diff")?;

    // Load first file (always an index)
    let (mut first_minimizers, header) = load_minimizers_from_path(first)?;
    info!("First index: loaded {} minimizers", first_minimizers.len());

    let before_count = first_minimizers.len();
    // Explicit w forces FASTX input.
    if window_size.is_none() && is_exact_index_file(second) {
        let (second_minimizers, second_header) = load_minimizers_from_path(second)?;
        info!(
            "Second index: loaded {} minimizers",
            second_minimizers.len()
        );

        // k must match. w must match or be 1 (every k-mer).
        if second_header.kmer_length != header.kmer_length
            || (second_header.window_size != header.window_size && second_header.window_size != 1)
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: second index has k={}, w={}, but first index has k={}, w={} (w must match or second index must be w=1)",
                second_header.kmer_length,
                second_header.window_size,
                header.kmer_length,
                header.window_size
            ));
        }
        first_minimizers.remove_all(&second_minimizers);
    } else {
        stream_diff_fastx(
            second,
            window_size.unwrap_or(header.window_size),
            &header,
            threads,
            &mut first_minimizers,
        )?;
    }

    // Report results
    info!(
        "Removed {} minimizers, {} remaining",
        before_count - first_minimizers.len(),
        first_minimizers.len()
    );

    write_exact_index_to_path(&mut first_minimizers, &header, output)?;

    let total_time = start_time.elapsed();
    info!("Completed diff in {:.2?}", total_time);

    Ok(())
}

/// Show index metadata.
pub fn info(index_path: &Path) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    info!("Deacon v{}; mode: info", version);

    if is_fuse_index_file(index_path) {
        let file = File::open(index_path)
            .context(format!("Failed to open index file {:?}", index_path))?;
        let mut reader = BufReader::new(file);
        let config = bincode::config::standard().with_fixed_int_encoding();
        let header: FuseIndexHeader = decode_from_std_read(&mut reader, config)
            .context("Failed to deserialise BFF header")?;
        header.validate()?;

        let file_size = std::fs::metadata(index_path).map(|m| m.len()).unwrap_or(0);
        let bits_per_key = if header.key_count > 0 {
            (file_size as f64 * 8.0) / header.key_count as f64
        } else {
            0.0
        };

        println!("Index information:");
        println!(
            "  Format: BFF (binary fuse filter, {}-bit fingerprints)",
            header.fingerprint_bits()
        );
        println!("  Format version: {}", header.format_version);
        println!("  K-mer length (k): {}", header.kmer_length);
        println!("  Window size (w): {}", header.window_size);
        println!("  Key count: {}", header.key_count);
        println!(
            "  File size: {} bytes (~{:.2} bits/key)",
            file_size, bits_per_key
        );
        let bits = header.fingerprint_bits();
        println!(
            "  False-positive rate: ~2^-{} (~{:.2e})",
            bits,
            2f64.powi(-(bits as i32))
        );
    } else {
        // Load exact index file
        let (minimizers, header) = load_minimizers_from_path(index_path)?;

        println!("Index information:");
        println!("  Format: exact (minimizer set)");
        println!("  Format version: {}", header.format_version);
        println!("  K-mer length (k): {}", header.kmer_length);
        println!("  Window size (w): {}", header.window_size);
        println!("  Distinct minimizer count: {}", minimizers.len());
    }

    let total_time = start_time.elapsed();
    info!("Loaded index info in {:.2?}", total_time);

    Ok(())
}

/// Drop low-complexity minimizers, or retain them if inverted.
pub fn filter(
    index_path: &Path,
    output: Option<&Path>,
    algorithm: crate::ComplexityAlgorithm,
    threshold: f32,
    invert: bool,
) -> Result<()> {
    crate::validate_unit_interval("complexity threshold", threshold)?;
    let start_time = Instant::now();
    if is_fuse_index_file(index_path) {
        anyhow::bail!("Complexity filtering is not supported on BFF indexes; use an exact index");
    }

    let (mut minimizers, header) = load_minimizers_from_path(index_path)?;
    let before = minimizers.len();
    info!(
        "Filtering index (k={}, w={}): {} minimizers, {} {} {}",
        header.kmer_length,
        header.window_size,
        before,
        algorithm,
        if invert { "<" } else { ">=" },
        threshold
    );

    minimizers.retain_complexity(header.kmer_length, algorithm, threshold, invert)?;
    let after = minimizers.len();

    write_exact_index_to_path(&mut minimizers, &header, output)?;

    info!(
        "Kept {} of {} minimizers ({} removed) in {:.2?}",
        after,
        before,
        before - after,
        start_time.elapsed()
    );

    Ok(())
}

/// Union minimizer indexes.
pub fn union(inputs: &[PathBuf], output: Option<&Path>) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    info!("Deacon v{}; mode: union", version);
    // Check input files
    if inputs.is_empty() {
        return Err(anyhow::anyhow!(
            "No input files provided for union operation"
        ));
    }

    // Read all headers first to determine total capacity needed
    let mut headers_and_counts = Vec::new();

    for path in inputs {
        reject_fuse_index(path, "union")?;
        let (header, count) = load_header_and_count(path)?;
        headers_and_counts.push((header, count));
    }

    // Get header from first file for output
    let header = &headers_and_counts[0].0;

    info!(
        "Combining indexes (k={}, w={})",
        header.kmer_length, header.window_size
    );

    // Verify all headers are compatible
    for (i, (file_header, _)) in headers_and_counts.iter().enumerate() {
        if file_header.kmer_length != header.kmer_length
            || file_header.window_size != header.window_size
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index {} has k={}, w={}, but first index has k={}, w={}",
                i,
                file_header.kmer_length,
                file_header.window_size,
                header.kmer_length,
                header.window_size
            ));
        }
    }

    // Load first index to determine type (u64 vs u128)
    let (mut all_minimizers, _) = load_minimizers_from_path(&inputs[0])?;
    info!("Index 1: loaded {} minimizers", all_minimizers.len());

    // Now load and merge remaining indexes
    for (i, path) in inputs.iter().enumerate().skip(1) {
        let (minimizers, _) = load_minimizers_from_path(path)?;
        let before_count = all_minimizers.len();

        // Merge minimizers (set union)
        all_minimizers.extend(minimizers);

        let expected_count = headers_and_counts[i].1;
        info!(
            "Index {}: {} minimizers, added {}, total: {}",
            i + 1,
            expected_count,
            all_minimizers.len() - before_count,
            all_minimizers.len()
        );
    }

    write_exact_index_to_path(&mut all_minimizers, header, output)?;

    let total_time = start_time.elapsed();
    info!(
        "United {} indexes with {} total minimizers in {:.2?}",
        inputs.len(),
        all_minimizers.len(),
        total_time
    );

    Ok(())
}

pub fn intersect(inputs: &[PathBuf], output: Option<&Path>) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    info!("Deacon v{}; mode: intersect", version);
    // Check inputs
    if inputs.len() < 2 {
        return Err(anyhow::anyhow!(
            "At least two input files are required for intersection operation"
        ));
    }

    // Read all headers first
    let mut headers_and_counts = Vec::new();

    for path in inputs {
        reject_fuse_index(path, "intersect")?;
        let (header, count) = load_header_and_count(path)?;
        headers_and_counts.push((header, count));
    }

    // Get header from first file for output
    let header = &headers_and_counts[0].0;

    info!(
        "Intersecting indexes (k={}, w={})",
        header.kmer_length, header.window_size
    );

    // Check header compat, allow w=1 passthrough
    for (i, (file_header, _)) in headers_and_counts.iter().enumerate() {
        if file_header.kmer_length != header.kmer_length
            || (file_header.window_size != header.window_size && file_header.window_size != 1)
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index {} has k={}, w={}, but first index has k={}, w={} (w must match or index must be w=1)",
                i,
                file_header.kmer_length,
                file_header.window_size,
                header.kmer_length,
                header.window_size
            ));
        }
    }

    // Load first index
    let (mut result_minimizers, _) = load_minimizers_from_path(&inputs[0])?;
    info!("Index 1: loaded {} minimizers", result_minimizers.len());

    // Intersect with remaining indexes
    for (i, path) in inputs.iter().enumerate().skip(1) {
        let (minimizers, _) = load_minimizers_from_path(path)?;

        // Intersect minimizers (set intersection)
        result_minimizers.intersect(&minimizers);

        let expected_count = headers_and_counts[i].1;
        info!(
            "Index {}: {} minimizers, retained {}, total: {}",
            i + 1,
            expected_count,
            result_minimizers.len(),
            result_minimizers.len()
        );
    }

    write_exact_index_to_path(&mut result_minimizers, header, output)?;

    let total_time = start_time.elapsed();
    info!(
        "Intersected {} indexes with {} common minimizers in {:.2?}",
        inputs.len(),
        result_minimizers.len(),
        total_time
    );

    Ok(())
}
