use crate::{FixedRapidHasher, RapidHashSet};
use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use serde::{Deserialize, Serialize};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

#[cfg(feature = "cli")]
use crate::IndexConfig;
#[cfg(feature = "cli")]
use crate::minimizers::{Buffers, KmerHasher};
#[cfg(feature = "cli")]
use paraseq::Record;
#[cfg(feature = "cli")]
use paraseq::prelude::{ParallelProcessor, ParallelReader};
#[cfg(feature = "cli")]
use parking_lot::Mutex;
#[cfg(feature = "cli")]
use rayon::prelude::*;
#[cfg(feature = "cli")]
use std::fs::File;
#[cfg(feature = "cli")]
use std::hash::BuildHasher;
#[cfg(feature = "cli")]
use std::path::PathBuf;
#[cfg(feature = "cli")]
use std::sync::{Arc, OnceLock};
#[cfg(feature = "cli")]
use std::time::Instant;

#[cfg(feature = "fetch")]
use indicatif::ProgressBar;

/// Index format version
pub const INDEX_FORMAT_VERSION: u8 = 3;

/// BFF format version. The on-disk magic is b"DBF" followed by the fingerprint
/// width in bits (16 or 32), so 16-/32-bit indexes are distinguishable.
pub const BFF_FORMAT_VERSION: u8 = 1;

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
        if k > 61 || k + w > 96 || !(k + w).is_multiple_of(2) {
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

/// Serialisable header for the BFF index file
#[derive(Serialize, Deserialize, Debug)]
pub struct BffHeader {
    /// b"DBF" prefix plus the fingerprint width in bits as the 4th byte
    pub magic: [u8; 4],
    pub format_version: u8,
    pub kmer_length: u8,
    pub window_size: u8,
    /// Distinct minimizers inserted (not BinaryFuse32::len())
    pub key_count: u64,
}

impl BffHeader {
    #[cfg(any(feature = "cli", test))]
    pub fn new(filter_bits: u8, kmer_length: u8, window_size: u8, key_count: u64) -> Self {
        BffHeader {
            magic: [b'D', b'B', b'F', filter_bits],
            format_version: BFF_FORMAT_VERSION,
            kmer_length,
            window_size,
            key_count,
        }
    }

    /// Fingerprint width in bits (16 or 32), encoded in the magic
    pub fn filter_bits(&self) -> u8 {
        self.magic[3]
    }

    /// Validate BFF header
    pub fn validate(&self) -> anyhow::Result<()> {
        if !self.magic.starts_with(b"DBF") {
            return Err(anyhow::anyhow!("Not a BFF index (bad magic bytes)"));
        }
        if self.format_version != BFF_FORMAT_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported BFF format version: {}",
                self.format_version
            ));
        }
        if !matches!(self.filter_bits(), 16 | 32) {
            return Err(anyhow::anyhow!(
                "Unsupported BFF fingerprint width: {} bits (expected 16 or 32)",
                self.filter_bits()
            ));
        }
        // BFF keys are raw u64 minimizers
        if self.kmer_length > 32 {
            return Err(anyhow::anyhow!(
                "Invalid BFF k-mer length: k={} (BFF supports k <= 32)",
                self.kmer_length
            ));
        }
        // Reuse exact k/w constraints
        IndexHeader::new(self.kmer_length, self.window_size).validate()
    }
}

/// Load just the header and count from an index file
#[cfg(feature = "cli")]
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

#[cfg(feature = "cli")]
static INDEX: OnceLock<(PathBuf, Arc<crate::MinimizerSet>, IndexHeader)> = OnceLock::new();

#[cfg(feature = "cli")]
pub fn current_index_path() -> Option<PathBuf> {
    INDEX.get().map(|(p, _m, _h)| p.clone())
}

#[cfg(feature = "cli")]
pub fn load_minimizers_cached(
    path: &Path,
) -> Result<(Arc<crate::MinimizerSet>, &'static IndexHeader)> {
    let (p, minimizers, header) = INDEX.get_or_init(|| {
        // Auto-detect exact vs BFF format
        let (m, h) = load_index_from_path_auto(path).unwrap();
        (path.to_owned(), Arc::new(m), h)
    });
    assert_eq!(
        p, path,
        "Currently, the server can only have one index loaded."
    );

    Ok((Arc::clone(minimizers), header))
}

/// Load minimizers from a reader (generic over any Read impl)
pub fn load_minimizers(reader: &mut impl Read) -> Result<(crate::MinimizerSet, IndexHeader)> {
    let mut reader = BufReader::with_capacity(1 << 20, reader);
    let config = bincode::config::standard().with_fixed_int_encoding();

    // Deserialise header
    let header: IndexHeader =
        decode_from_std_read(&mut reader, config).context("Failed to deserialise index header")?;
    header.validate()?;

    // Deserialise the count of minimizers (stored as u64 for cross-platform compatibility)
    let count: u64 = decode_from_std_read(&mut reader, config)
        .context("Failed to deserialise minimizer count")?;
    let count = count as usize;

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
            reader.read_exact(batch).context(format!(
                "Failed to load minimizer batch {}, index may be corrupt",
                i / B
            ))?;

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
            reader.read_exact(batch).context(format!(
                "Failed to load minimizer batch {}, index may be corrupt",
                i / B
            ))?;

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

    // Validate that we loaded the expected number of minimizers
    let loaded_count = minimizers.len();
    if loaded_count != count {
        return Err(anyhow::anyhow!(
            "Failed to load expected number of minimizers; expected {} and observed {}. Index may be corrupt",
            count,
            loaded_count
        ));
    }

    Ok((minimizers, header))
}

/// Load minimizers from an index file path
pub fn load_minimizers_from_path(path: &Path) -> Result<(crate::MinimizerSet, IndexHeader)> {
    let mut file =
        std::fs::File::open(path).context(format!("Failed to open index file {:?}", path))?;
    load_minimizers(&mut file)
}

/// Load a BFF index from a reader (synthesizes an IndexHeader from k/w)
pub fn load_bff(reader: &mut impl Read) -> Result<(crate::MinimizerSet, IndexHeader)> {
    let mut reader = BufReader::with_capacity(1 << 20, reader);
    let config = bincode::config::standard().with_fixed_int_encoding();

    // Header via serde, filter via native bincode
    let header: BffHeader =
        decode_from_std_read(&mut reader, config).context("Failed to deserialise BFF header")?;
    header.validate()?;

    // Decode the filter at the width recorded in the magic
    let filter = match header.filter_bits() {
        16 => crate::FuseFilterKind::Bits16(
            bincode::decode_from_std_read(&mut reader, config)
                .context("Failed to deserialise BFF filter, index may be corrupt")?,
        ),
        32 => crate::FuseFilterKind::Bits32(
            bincode::decode_from_std_read(&mut reader, config)
                .context("Failed to deserialise BFF filter, index may be corrupt")?,
        ),
        bits => {
            return Err(anyhow::anyhow!(
                "Unsupported BFF fingerprint width: {} bits (expected 16 or 32)",
                bits
            ));
        }
    };

    let index_header = IndexHeader::new(header.kmer_length, header.window_size);
    let set = crate::MinimizerSet::Fuse(crate::FuseFilter {
        filter,
        key_count: header.key_count as usize,
    });
    Ok((set, index_header))
}

/// Load an index from a reader, auto-detecting exact vs BFF format from the magic bytes
pub fn load_index_auto(reader: &mut impl Read) -> Result<(crate::MinimizerSet, IndexHeader)> {
    let mut magic = [0u8; 4];
    reader
        .read_exact(&mut magic)
        .context("Failed to read index header, file may be empty or truncated")?;
    // Re-chain the consumed magic bytes ahead of the stream
    let mut combined = (&magic[..]).chain(reader);
    if magic.starts_with(b"DBF") {
        load_bff(&mut combined)
    } else {
        load_minimizers(&mut combined)
    }
}

/// Load an index from a path, auto-detecting exact vs BFF format
pub fn load_index_from_path_auto(path: &Path) -> Result<(crate::MinimizerSet, IndexHeader)> {
    let mut file =
        std::fs::File::open(path).context(format!("Failed to open index file {:?}", path))?;
    load_index_auto(&mut file)
}

/// Minimizer bytes buffered before hitting the writer
const WRITE_CHUNK: usize = 4 << 20;

/// Write minimizers in byte-aligned packed format after the header and count
fn write_minimizers<K: Key>(
    minimizers: impl Iterator<Item = K>,
    count: usize,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    let mut writer: BufWriter<Box<dyn Write>> = match output_path {
        Some(path) if path.as_os_str() != "-" => BufWriter::new(Box::new(
            std::fs::File::create(path).context("Failed to create output file")?,
        )),
        _ => BufWriter::new(Box::new(io::stdout())),
    };

    // Use fixed-width little-endian encoding for all values.
    let config = bincode::config::standard().with_fixed_int_encoding();
    encode_into_std_write(header, &mut writer, config)
        .context("Failed to serialise index header")?;
    // Serialise the count of minimizers first (as u64 for cross-platform compatibility)
    encode_into_std_write(count as u64, &mut writer, config)
        .context("Failed to serialise minimizer count")?;

    let bytes = (header.kmer_length as usize).div_ceil(4);
    let mut buf: Vec<u8> = Vec::with_capacity(WRITE_CHUNK + 16);
    for minimizer in minimizers {
        minimizer.push_le(&mut buf, bytes);
        if buf.len() >= WRITE_CHUNK {
            writer
                .write_all(&buf)
                .context("Failed to write minimizers")?;
            buf.clear();
        }
    }
    writer
        .write_all(&buf)
        .context("Failed to write minimizers")?;
    writer.flush().context("Failed to flush index")
}

pub fn dump_minimizers(
    minimizers: &crate::MinimizerSet,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    match minimizers {
        crate::MinimizerSet::U64(set) => dump_set(set, header, output_path),
        crate::MinimizerSet::U128(set) => dump_set(set, header, output_path),
        crate::MinimizerSet::Fuse(_) => Err(anyhow::anyhow!(
            "Cannot serialise a BFF index in the exact index format"
        )),
    }
}

#[cfg(feature = "cli")]
fn dump_set<K: Key>(
    set: &RapidHashSet<K>,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    let order = LoadOrder::from_set(set);
    write_minimizers(order.iter(), set.len(), header, output_path)
}

/// Without threads there is nothing to gain from ordering the output
#[cfg(not(feature = "cli"))]
fn dump_set<K: Key>(
    set: &RapidHashSet<K>,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    write_minimizers(set.iter().copied(), set.len(), header, output_path)
}

/// Dump indexed minimizers to FASTA
/// Detect a BFF index by its magic bytes
#[cfg(feature = "cli")]
fn is_bff_file(path: &Path) -> bool {
    let mut magic = [0u8; 4];
    File::open(path)
        .ok()
        .map(|mut f| f.read_exact(&mut magic).is_ok() && magic.starts_with(b"DBF"))
        .unwrap_or(false)
}

/// Reject a binary fuse filter (.pidx) index for commands that require an exact (.idx) index
#[cfg(feature = "cli")]
fn reject_bff(path: &Path, operation: &str) -> Result<()> {
    if is_bff_file(path) {
        return Err(anyhow::anyhow!(
            "{:?} is a binary fuse filter (.pidx) index; {} requires an exact (.idx) index",
            path,
            operation
        ));
    }
    Ok(())
}

#[cfg(feature = "cli")]
pub fn dump(index_path: &Path, output_path: Option<&Path>) -> Result<()> {
    if is_bff_file(index_path) {
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
        crate::MinimizerSet::Fuse(_) => unreachable!("BFF dump is rejected by is_bff_file above"),
    }

    writer.flush()?;
    Ok(())
}

/// Freeze an exact index into a BFF (binary fuse filter) index (k<=32)
#[cfg(feature = "cli")]
pub fn freeze(index_path: &Path, output_path: Option<&Path>, bits: u8) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    eprintln!("Deacon v{}; mode: freeze", version);

    if !matches!(bits, 16 | 32) {
        return Err(anyhow::anyhow!(
            "Unsupported BFF fingerprint width: {} bits (expected 16 or 32)",
            bits
        ));
    }

    reject_bff(index_path, "freeze")?;

    let (minimizers, header) =
        load_minimizers_from_path(index_path).context("Failed to load index")?;

    // Collect unique u64 keys
    let keys: Vec<u64> = match minimizers {
        crate::MinimizerSet::U64(set) => set.into_iter().collect(),
        crate::MinimizerSet::U128(_) => {
            return Err(anyhow::anyhow!(
                "BFF supports k <= 32 (u64 minimizers); got k={} (u128). Build a k<=32 index first.",
                header.kmer_length
            ));
        }
        crate::MinimizerSet::Fuse(_) => {
            return Err(anyhow::anyhow!("Input is already a BFF index"));
        }
    };

    let key_count = keys.len();
    eprintln!(
        "Building {}-bit binary fuse filter from {} minimizers (k={}, w={})",
        bits, key_count, header.kmer_length, header.window_size
    );

    let bff_header = BffHeader::new(
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
    encode_into_std_write(&bff_header, &mut writer, config)
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
    eprintln!(
        "Wrote BFF: {} keys, {} filter bytes ({:.2} bits/key); false-positive rate ~2^-{}",
        key_count, filter_bytes, bits_per_key, bits
    );
    eprintln!("Completed freeze in {:.2?}", start_time.elapsed());
    Ok(())
}

#[cfg(feature = "cli")]
fn reader_with_inferred_batch_size(
    in_path: Option<&Path>,
) -> Result<paraseq::fastx::Reader<Box<dyn Read + Send>>> {
    let mut reader = paraseq::ReaderBuilder::optional_path(in_path)
        .build()
        .map_err(|e| anyhow::anyhow!("Failed to open input: {}", e))?;
    reader.update_batch_size_in_bp(256 * 1024)?;
    Ok(reader)
}

#[cfg(feature = "cli")]
use crate::filter::ProcessingStats;

/// Parts minimizers are partitioned into: by k-mer value for deduplication, then by hash
/// bucket range for writing in load order
#[cfg(feature = "cli")]
const SHARD_BITS: u32 = 10;
#[cfg(feature = "cli")]
const SHARDS: usize = 1 << SHARD_BITS;
/// Values a thread stages for one part before flushing them
#[cfg(feature = "cli")]
const STAGE_CAP: usize = 1024;

/// Bits a k-mer of length `kmer_length` occupies
#[cfg(feature = "cli")]
fn used_bits(kmer_length: u8) -> u32 {
    2 * kmer_length as u32
}
/// Values a shard may hold before it is sorted and deduplicated in place
#[cfg(feature = "cli")]
const COMPACT_MIN: usize = 1 << 22;

/// Minimizer width handled by the index writers
#[cfg_attr(not(feature = "cli"), allow(dead_code))]
pub(crate) trait Key: Copy + Ord + std::hash::Hash + Send + Sync + 'static {
    fn from_minimizers(minimizers: &crate::MinimizerVec) -> &Vec<Self>;
    fn shard(self, shift: u32) -> usize;
    fn push_le(self, out: &mut Vec<u8>, bytes: usize);
}

impl Key for u64 {
    #[inline]
    fn from_minimizers(minimizers: &crate::MinimizerVec) -> &Vec<Self> {
        match minimizers {
            crate::MinimizerVec::U64(vec) => vec,
            crate::MinimizerVec::U128(_) => unreachable!("u128 minimizers in a u64 build"),
        }
    }
    #[inline]
    fn shard(self, shift: u32) -> usize {
        (self >> shift) as usize
    }
    #[inline]
    fn push_le(self, out: &mut Vec<u8>, bytes: usize) {
        out.extend_from_slice(&self.to_le_bytes()[..bytes]);
    }
}

impl Key for u128 {
    #[inline]
    fn from_minimizers(minimizers: &crate::MinimizerVec) -> &Vec<Self> {
        match minimizers {
            crate::MinimizerVec::U128(vec) => vec,
            crate::MinimizerVec::U64(_) => unreachable!("u64 minimizers in a u128 build"),
        }
    }
    #[inline]
    fn shard(self, shift: u32) -> usize {
        (self >> shift) as usize
    }
    #[inline]
    fn push_le(self, out: &mut Vec<u8>, bytes: usize) {
        out.extend_from_slice(&self.to_le_bytes()[..bytes]);
    }
}

#[cfg(feature = "cli")]
struct Shard<K> {
    values: Vec<K>,
    /// Length at which the shard is compacted, doubling after each compaction
    limit: usize,
}

/// Minimizers partitioned by value range. Threads append to their own shard buffers
/// and only take a shard's lock to hand over a full buffer, so there is no global
/// set to contend on and no hashing: duplicates go once the shards are sorted.
#[cfg(feature = "cli")]
struct Shards<K> {
    shards: Vec<Mutex<Shard<K>>>,
    shift: u32,
}

#[cfg(feature = "cli")]
impl<K: Key> Shards<K> {
    fn new(kmer_length: u8) -> Self {
        Shards {
            shards: (0..SHARDS)
                .map(|_| {
                    Mutex::new(Shard {
                        values: Vec::new(),
                        limit: COMPACT_MIN,
                    })
                })
                .collect(),
            shift: used_bits(kmer_length).saturating_sub(SHARD_BITS),
        }
    }

    /// Hand a full staging buffer to its shard, compacting the shard if it has grown
    #[inline]
    fn flush(&self, shard: usize, staged: &mut Vec<K>) {
        let mut shard = self.shards[shard].lock();
        shard.values.extend_from_slice(staged);
        staged.clear();
        if shard.values.len() >= shard.limit {
            shard.values.sort_unstable();
            shard.values.dedup();
            shard.limit = (2 * shard.values.len()).max(COMPACT_MIN);
        }
    }

    /// Sort and deduplicate every shard in parallel
    fn finish(self) -> Vec<Vec<K>> {
        let mut values: Vec<Vec<K>> = self
            .shards
            .into_iter()
            .map(|shard| shard.into_inner().values)
            .collect();
        values.par_iter_mut().for_each(|shard| {
            shard.sort_unstable();
            shard.dedup();
        });
        values
    }
}

#[cfg(feature = "cli")]
fn staging<K>() -> Vec<Vec<K>> {
    (0..SHARDS).map(|_| Vec::with_capacity(STAGE_CAP)).collect()
}

/// Minimizers grouped by the slice of the hash table `load_minimizers` inserts them
/// into, so that its inserts stay local instead of jumping across the whole table. Each
/// group is sorted by value, which is uncorrelated with bucket order: that avoids
/// inserting in exact bucket order, which clusters, and makes the output deterministic.
/// Only affects how fast an index loads, never its contents.
#[cfg(feature = "cli")]
struct LoadOrder<K>(Vec<Vec<Vec<K>>>);

#[cfg(feature = "cli")]
impl<K: Key> LoadOrder<K> {
    /// Split value-sorted runs in parallel, keeping each group sorted, and freeing each
    /// run once it is split
    fn from_sorted(runs: Vec<Vec<K>>, count: usize) -> Self {
        LoadOrder(runs.into_par_iter().map(|run| split(&run, count)).collect())
    }

    fn from_set(set: &RapidHashSet<K>) -> Self {
        let mut groups = split(set, set.len());
        groups
            .par_iter_mut()
            .for_each(|group| group.sort_unstable());
        LoadOrder(vec![groups])
    }

    fn iter(&self) -> impl Iterator<Item = K> + '_ {
        (0..SHARDS).flat_map(|g| self.0.iter().flat_map(move |run| run[g].iter().copied()))
    }
}

#[cfg(feature = "cli")]
fn split<'a, K: Key>(minimizers: impl IntoIterator<Item = &'a K>, count: usize) -> Vec<Vec<K>> {
    // hashbrown allocates a power of two buckets at a 7/8 load factor
    let buckets = (count * 8 / 7).next_power_of_two();
    let shift = buckets.trailing_zeros().saturating_sub(SHARD_BITS);
    let mut groups = vec![Vec::new(); SHARDS];
    for &minimizer in minimizers {
        let bucket = FixedRapidHasher.hash_one(minimizer) as usize & (buckets - 1);
        groups[bucket >> shift].push(minimizer);
    }
    groups
}

#[cfg(feature = "cli")]
struct BuildIndexProcessor<'c, K: Key> {
    config: &'c IndexConfig,
    hasher: KmerHasher,
    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    /// Minimizers waiting to be handed over, one buffer per shard
    staged: Vec<Vec<K>>,
    // Global state
    shards: Arc<Shards<K>>,
    global_stats: Arc<Mutex<ProcessingStats>>,
}

#[cfg(feature = "cli")]
impl<K: Key> Clone for BuildIndexProcessor<'_, K> {
    fn clone(&self) -> Self {
        BuildIndexProcessor {
            config: self.config,
            hasher: self.hasher.clone(),
            buffers: self.buffers.clone(),
            local_stats: ProcessingStats::default(),
            staged: staging(),
            shards: Arc::clone(&self.shards),
            global_stats: Arc::clone(&self.global_stats),
        }
    }
}

#[cfg(feature = "cli")]
impl<Rf: Record, K: Key> ParallelProcessor<Rf> for BuildIndexProcessor<'_, K> {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        crate::minimizers::fill_minimizers(
            &seq,
            &self.hasher,
            self.config.kmer_length,
            self.config.window_size,
            &mut self.buffers,
        );

        let Self {
            buffers,
            staged,
            shards,
            ..
        } = self;
        let shift = shards.shift;
        for &minimizer in K::from_minimizers(&buffers.minimizers) {
            let shard = minimizer.shard(shift);
            let staged = &mut staged[shard];
            staged.push(minimizer);
            if staged.len() == STAGE_CAP {
                shards.flush(shard, staged);
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        // Tick to stderr once every Gbp
        let mut stats = self.global_stats.lock();
        stats.total_seqs += self.local_stats.total_seqs;
        stats.total_bp += self.local_stats.total_bp;

        if !self.config.quiet {
            let current_gb = stats.total_bp / 1_000_000_000;
            if current_gb > stats.last_reported {
                eprintln!(
                    "  Processed {} sequences ({}bp)",
                    stats.total_seqs, stats.total_bp
                );
                stats.last_reported = current_gb;
            }
        }

        self.local_stats = ProcessingStats::default();
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let Self { staged, shards, .. } = self;
        for (shard, staged) in staged.iter_mut().enumerate() {
            if !staged.is_empty() {
                shards.flush(shard, staged);
            }
        }
        Ok(())
    }
}

/// Build an index of minimizers from a fastx file
#[cfg(feature = "cli")]
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

    let header = IndexHeader::new(config.kmer_length, config.window_size);
    if config.kmer_length <= 32 {
        build_shards::<u64>(config, reader, &header)?;
    } else {
        build_shards::<u128>(config, reader, &header)?;
    }

    eprintln!("Completed build in {:.2?}", start_time.elapsed());
    Ok(())
}

#[cfg(feature = "cli")]
fn build_shards<K: Key>(
    config: &IndexConfig,
    reader: paraseq::fastx::Reader<Box<dyn Read + Send>>,
    header: &IndexHeader,
) -> Result<()> {
    let mut processor = BuildIndexProcessor::<K> {
        config,
        hasher: KmerHasher::new(config.kmer_length as usize),
        buffers: if config.kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        },
        local_stats: ProcessingStats::default(),
        staged: staging(),
        shards: Arc::new(Shards::new(config.kmer_length)),
        global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
    };
    reader.process_parallel(&mut processor, config.threads as usize)?;

    let stats = Arc::into_inner(processor.global_stats)
        .expect("stats outlived the workers")
        .into_inner();
    let shards = Arc::into_inner(processor.shards).expect("shards outlived the workers");

    // Own pool for the parallel phases, so the global one is left to the caller
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads as usize)
        .build()
        .context("Failed to initialise thread pool")?;
    let (count, order) = pool.install(|| {
        let shards = shards.finish();
        let count: usize = shards.iter().map(|shard| shard.len()).sum();
        (count, LoadOrder::from_sorted(shards, count))
    });
    eprintln!(
        "Indexed {} minimizers from {} record(s) ({}bp)",
        count, stats.total_seqs, stats.total_bp
    );

    write_minimizers(order.iter(), count, header, config.output_path.as_deref())
}

/// Minimizers found in the index being diffed
#[cfg(feature = "cli")]
#[derive(Clone)]
enum HitSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

#[cfg(feature = "cli")]
impl HitSet {
    /// An empty hit set matching the width of `set`
    fn empty_like(set: &crate::MinimizerSet) -> Self {
        match set {
            crate::MinimizerSet::U64(_) => HitSet::U64(RapidHashSet::default()),
            crate::MinimizerSet::U128(_) => HitSet::U128(RapidHashSet::default()),
            // load_minimizers rejects BFF files, so this is unreachable
            crate::MinimizerSet::Fuse(_) => unreachable!("diff does not operate on BFF indexes"),
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
    fn remove_from(&self, set: &mut crate::MinimizerSet) {
        match (self, set) {
            (HitSet::U64(hits), crate::MinimizerSet::U64(set)) => {
                for minimizer in hits {
                    set.remove(minimizer);
                }
            }
            (HitSet::U128(hits), crate::MinimizerSet::U128(set)) => {
                for minimizer in hits {
                    set.remove(minimizer);
                }
            }
            _ => unreachable!("minimizer width mismatch between hit set and index"),
        }
    }
}

#[cfg(feature = "cli")]
#[derive(Clone)]
struct DiffIndexProcessor<'a> {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    /// Index being subtracted from. Read-only while streaming, so every thread probes
    /// it without locking; the removal itself is a single pass once streaming is done.
    first: &'a crate::MinimizerSet,
    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    /// Hits seen by this thread since the last batch. Bounded by `first`, and usually a
    /// tiny fraction of it.
    local_hits: HitSet,
    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_hits: Arc<Mutex<HitSet>>,
}

#[cfg(feature = "cli")]
impl<Rf: Record> ParallelProcessor<Rf> for DiffIndexProcessor<'_> {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        if seq.len() < self.kmer_length as usize {
            return Ok(());
        }

        crate::minimizers::fill_minimizers_unchecked(
            &seq,
            &self.hasher,
            self.kmer_length,
            self.window_size,
            &mut self.buffers,
        );

        // Dispatch on width once, then probe lock-free, keeping only the hits
        match (self.first, &self.buffers.minimizers, &mut self.local_hits) {
            (crate::MinimizerSet::U64(first), crate::MinimizerVec::U64(vec), HitSet::U64(hits)) => {
                for &minimizer in vec.iter() {
                    if first.contains(&minimizer) {
                        hits.insert(minimizer);
                    }
                }
            }
            (
                crate::MinimizerSet::U128(first),
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
                eprintln!(
                    "  Processed {} sequences ({}bp), removed {} minimizers",
                    stats.total_seqs, stats.total_bp, hits
                );
                stats.last_reported = current_10gb;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

/// Stream minimizers from a FASTX file or stdin and remove those present in first_minimizers
#[cfg(feature = "cli")]
fn stream_diff_fastx(
    fastx_path: &Path,
    window_size: u8,
    first_header: &IndexHeader,
    threads: u16,
    first_minimizers: &mut crate::MinimizerSet,
) -> Result<(usize, usize)> {
    let path = fastx_path;
    let kmer_length = first_header.kmer_length();

    // Validate k-mer and window size constraints
    let temp_config = crate::IndexConfig {
        input_path: PathBuf::new(),
        kmer_length,
        window_size,
        output_path: None,
        threads: 0,
        quiet: false,
    };
    temp_config.validate()?;

    // w must match or be 1 (w=1 emits every k-mer, for exact masking)
    if window_size != first_header.window_size() && window_size != 1 {
        return Err(anyhow::anyhow!(
            "FASTX w={} must match first index w={} or be 1 (for exact k-mer subtraction)",
            window_size,
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

    // Read only index while streaming - remove minimizers in single final op
    let global_stats = Arc::new(Mutex::new(ProcessingStats::default()));
    let global_hits = Arc::new(Mutex::new(HitSet::empty_like(first_minimizers)));
    let start_time = Instant::now();

    let mut processor = DiffIndexProcessor {
        kmer_length,
        window_size,
        hasher: KmerHasher::new(kmer_length as usize),
        first: first_minimizers,
        buffers: if kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        },
        local_stats: ProcessingStats::default(),
        local_hits: HitSet::empty_like(first_minimizers),
        global_stats: global_stats.clone(),
        global_hits: global_hits.clone(),
    };

    reader.process_parallel(&mut processor, threads as usize)?;
    drop(processor);

    let stats = global_stats.lock().clone();
    global_hits.lock().remove_from(first_minimizers);

    let elapsed = start_time.elapsed();
    eprintln!(
        "Processed {} sequences ({}bp) from FASTX file in {:.2?} ({:.1} Mbp/s)",
        stats.total_seqs,
        stats.total_bp,
        elapsed,
        stats.total_bp as f64 / elapsed.as_secs_f64() / 1_000_000.0
    );

    Ok((stats.total_seqs as usize, stats.total_bp as usize))
}

/// Compute the set difference between two minimizer indexes (A - B)
#[cfg(feature = "cli")]
pub fn diff(
    first: &Path,
    second: &Path,
    window_size: Option<u8>,
    threads: u16,
    output: Option<&Path>,
) -> Result<()> {
    let start_time = Instant::now();

    reject_bff(first, "diff")?;
    reject_bff(second, "diff")?;

    // Load first file (always an index)
    let (mut first_minimizers, header) = load_minimizers_from_path(first)?;
    eprintln!("First index: loaded {} minimizers", first_minimizers.len());

    // Guess if second file is an index or FASTX file
    let second_minimizers = if let Some(w) = window_size {
        // An explicit window marks the second file as FASTX; k comes from the first index
        let before_count = first_minimizers.len();
        let (_seq_count, _total_bp) =
            stream_diff_fastx(second, w, &header, threads, &mut first_minimizers)?;

        // Report results
        eprintln!(
            "Removed {} minimizers, {} remaining",
            before_count - first_minimizers.len(),
            first_minimizers.len()
        );

        dump_minimizers(&first_minimizers, &header, output)?;

        let total_time = start_time.elapsed();
        eprintln!("Completed diff in {:.2?}", total_time);

        return Ok(());
    } else {
        // Try to load as index file first
        if let Ok((second_minimizers, second_header)) = load_minimizers_from_path(second) {
            // Second file is an index file
            eprintln!(
                "Second index: loaded {} minimizers",
                second_minimizers.len()
            );

            // k must match; w must match or be 1 (w=1 means every k-mer)
            if second_header.kmer_length() != header.kmer_length()
                || (second_header.window_size() != header.window_size()
                    && second_header.window_size() != 1)
            {
                return Err(anyhow::anyhow!(
                    "Incompatible headers: second index has k={}, w={}, but first index has k={}, w={} (w must match or second index must be w=1)",
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
            let w = header.window_size();

            // Count minimizers before diff
            let before_count = first_minimizers.len();

            let (_seq_count, _total_bp) =
                stream_diff_fastx(second, w, &header, threads, &mut first_minimizers)?;

            // Report results
            eprintln!(
                "Removed {} minimizers, {} remaining",
                before_count - first_minimizers.len(),
                first_minimizers.len()
            );

            dump_minimizers(&first_minimizers, &header, output)?;

            let total_time = start_time.elapsed();
            eprintln!("Completed diff in {:.2?}", total_time);

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
    eprintln!("Completed diff in {:.2?}", total_time);

    Ok(())
}

/// Show info about an index
#[cfg(feature = "cli")]
pub fn info(index_path: &Path) -> Result<()> {
    let start_time = Instant::now();

    if is_bff_file(index_path) {
        let file = File::open(index_path)
            .context(format!("Failed to open index file {:?}", index_path))?;
        let mut reader = BufReader::new(file);
        let config = bincode::config::standard().with_fixed_int_encoding();
        let header: BffHeader = decode_from_std_read(&mut reader, config)
            .context("Failed to deserialise BFF header")?;
        header.validate()?;

        let file_size = std::fs::metadata(index_path).map(|m| m.len()).unwrap_or(0);
        let bits_per_key = if header.key_count > 0 {
            (file_size as f64 * 8.0) / header.key_count as f64
        } else {
            0.0
        };

        eprintln!("Index information:");
        eprintln!(
            "  Format: BFF (binary fuse filter, {}-bit fingerprints)",
            header.filter_bits()
        );
        eprintln!("  Format version: {}", header.format_version);
        eprintln!("  K-mer length (k): {}", header.kmer_length);
        eprintln!("  Window size (w): {}", header.window_size);
        eprintln!("  Key count: {}", header.key_count);
        eprintln!(
            "  File size: {} bytes (~{:.2} bits/key)",
            file_size, bits_per_key
        );
        let bits = header.filter_bits();
        eprintln!(
            "  False-positive rate: ~2^-{} (~{:.2e})",
            bits,
            2f64.powi(-(bits as i32))
        );
    } else {
        // Load exact index file
        let (minimizers, header) = load_minimizers_from_path(index_path)?;

        eprintln!("Index information:");
        eprintln!("  Format: exact (minimizer set)");
        eprintln!("  Format version: {}", header.format_version);
        eprintln!("  K-mer length (k): {}", header.kmer_length());
        eprintln!("  Window size (w): {}", header.window_size());
        eprintln!("  Distinct minimizer count: {}", minimizers.len());
    }

    let total_time = start_time.elapsed();
    eprintln!("Loaded index info in {:.2?}", total_time);

    Ok(())
}

/// Discard minimizers below a complexity threshold (or keep only those below, if inverted)
#[cfg(feature = "cli")]
pub fn filter(
    index_path: &Path,
    output: Option<&Path>,
    algorithm: crate::ComplexityAlgorithm,
    threshold: f32,
    invert: bool,
) -> Result<()> {
    crate::validate_unit_interval("complexity threshold", threshold)?;
    let start_time = Instant::now();
    if is_bff_file(index_path) {
        anyhow::bail!("Complexity filtering is not supported on BFF indexes; use an exact index");
    }

    let (mut minimizers, header) = load_minimizers_from_path(index_path)?;
    let before = minimizers.len();
    eprintln!(
        "Filtering index (k={}, w={}): {} minimizers, {} {} {}",
        header.kmer_length(),
        header.window_size(),
        before,
        algorithm,
        if invert { "<" } else { ">=" },
        threshold
    );

    minimizers.retain_complexity(header.kmer_length(), algorithm, threshold, invert)?;
    let after = minimizers.len();

    dump_minimizers(&minimizers, &header, output)?;

    eprintln!(
        "Kept {} of {} minimizers ({} removed) in {:.2?}",
        after,
        before,
        before - after,
        start_time.elapsed()
    );

    Ok(())
}

/// Combine minimizer indexes (set union)
#[cfg(feature = "cli")]
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
        reject_bff(path, "union")?;
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
    let (mut all_minimizers, _) = load_minimizers_from_path(&inputs[0])?;
    eprintln!("Index 1: loaded {} minimizers", all_minimizers.len());

    // Now load and merge remaining indexes
    for (i, path) in inputs.iter().enumerate().skip(1) {
        let (minimizers, _) = load_minimizers_from_path(path)?;
        let before_count = all_minimizers.len();

        // Merge minimizers (set union)
        all_minimizers.extend(minimizers);

        let expected_count = headers_and_counts[i].1;
        eprintln!(
            "Index {}: {} minimizers, added {}, total: {}",
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

#[cfg(feature = "cli")]
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
        reject_bff(path, "intersect")?;
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

    // Check header compat, allow w=1 passthrough
    for (i, (file_header, _)) in headers_and_counts.iter().enumerate() {
        if file_header.kmer_length() != header.kmer_length()
            || (file_header.window_size() != header.window_size() && file_header.window_size() != 1)
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index {} has k={}, w={}, but first index has k={}, w={} (w must match or index must be w=1)",
                i,
                file_header.kmer_length(),
                file_header.window_size(),
                header.kmer_length(),
                header.window_size()
            ));
        }
    }

    // Load first index
    let (mut result_minimizers, _) = load_minimizers_from_path(&inputs[0])?;
    eprintln!("Index 1: loaded {} minimizers", result_minimizers.len());

    // Intersect with remaining indexes
    for (i, path) in inputs.iter().enumerate().skip(1) {
        let (minimizers, _) = load_minimizers_from_path(path)?;

        // Intersect minimizers (set intersection)
        result_minimizers.intersect(&minimizers);

        let expected_count = headers_and_counts[i].1;
        eprintln!(
            "Index {}: {} minimizers, retained {}, total: {}",
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

    let mut response = minreq::get(&url)
        .send_lazy()
        .context("Failed to download index")?;

    if response.status_code != 200 {
        anyhow::bail!("Failed to fetch index: HTTP {}", response.status_code);
    }

    let content_length = response
        .headers
        .iter()
        .find(|(name, _)| name == "content-length")
        .and_then(|(_, value)| value.parse::<u64>().ok())
        .unwrap_or(0);

    let pb = ProgressBar::new(content_length);

    let output_path = output
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| std::path::PathBuf::from(&filename));

    let mut temp_path = output_path.clone();
    temp_path.as_mut_os_string().push(".tmp");

    let mut file = std::fs::File::create(&temp_path).context("Failed to create temporary file")?;
    std::io::copy(&mut pb.wrap_read(&mut response), &mut file)
        .context("Failed to write index to file")?;

    std::fs::rename(&temp_path, &output_path).context("Failed to finalise index")?;

    pb.finish_and_clear();
    eprintln!("Index saved to: {}", output_path.display());

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

    #[test]
    fn test_bff_header_validation() {
        // Both 16- and 32-bit are valid widths
        assert!(BffHeader::new(16, 31, 15, 100).validate().is_ok());
        assert!(BffHeader::new(32, 31, 15, 100).validate().is_ok());

        // Bad magic
        let mut bad_magic = BffHeader::new(16, 31, 15, 100);
        bad_magic.magic = *b"XXXX";
        assert!(bad_magic.validate().is_err());

        // Bad version
        let mut bad_ver = BffHeader::new(16, 31, 15, 100);
        bad_ver.format_version = 99;
        assert!(bad_ver.validate().is_err());

        // Unsupported filter width
        let mut bad_bits = BffHeader::new(16, 31, 15, 100);
        bad_bits.magic[3] = 8;
        assert!(bad_bits.validate().is_err());

        // k > 32 is rejected (BFF keys on raw u64 minimizers)
        assert!(BffHeader::new(16, 41, 21, 100).validate().is_err());
    }

    #[test]
    fn test_bff_roundtrip_membership() {
        let keys: Vec<u64> = (0..10_000u64)
            .map(|i| i.wrapping_mul(0x9E3779B97F4A7C15))
            .collect();
        let filter = xorf::BinaryFuse32::try_from(&keys).unwrap();
        let header = BffHeader::new(32, 31, 15, keys.len() as u64);

        let config = bincode::config::standard().with_fixed_int_encoding();
        let mut buf = Vec::new();
        encode_into_std_write(&header, &mut buf, config).unwrap();
        bincode::encode_into_std_write(&filter, &mut buf, config).unwrap();

        let mut cursor = std::io::Cursor::new(&buf);
        let (set, idx_header) = load_bff(&mut cursor).unwrap();
        assert_eq!(idx_header.kmer_length(), 31);
        assert_eq!(idx_header.window_size(), 15);
        assert_eq!(set.len(), keys.len());
        assert!(set.is_u64());
        // No false negatives
        for &k in &keys {
            assert!(set.contains_u64(k));
        }

        // auto-detect dispatches to BFF
        let mut cursor2 = std::io::Cursor::new(&buf);
        let (set2, _) = load_index_auto(&mut cursor2).unwrap();
        assert!(matches!(set2, crate::MinimizerSet::Fuse(_)));
    }
}
