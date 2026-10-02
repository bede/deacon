//! Exact/BFF index headers and serialization.

use crate::index::{FixedRapidHasher, FuseFilter, FuseFilterKind, IndexStorage, RapidHashSet};
use crate::minimizers::{MinimizerShards, validate_k_w};
use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use serde::{Deserialize, Serialize};
use std::hash::BuildHasher;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use tracing::info;

/// Index format version
pub const INDEX_FORMAT_VERSION: u8 = 3;

/// BFF format version. The on-disk magic is b"DBF" followed by the fingerprint
/// width in bits (16 or 32), so 16-/32-bit indexes are distinguishable.
pub const FUSE_FORMAT_VERSION: u8 = 1;

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

    pub fn validate(&self) -> anyhow::Result<()> {
        if self.format_version != INDEX_FORMAT_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported index format version: {}",
                self.format_version
            ));
        }
        validate_k_w(self.kmer_length, self.window_size)
    }
}

/// Serialisable header for the BFF index file
#[derive(Serialize, Deserialize, Debug)]
pub struct FuseIndexHeader {
    /// b"DBF" prefix plus the fingerprint width in bits as the 4th byte
    pub magic: [u8; 4],
    pub format_version: u8,
    pub kmer_length: u8,
    pub window_size: u8,
    /// Distinct minimizers inserted (not BinaryFuse32::len())
    pub key_count: u64,
}

impl FuseIndexHeader {
    pub fn new(fingerprint_bits: u8, kmer_length: u8, window_size: u8, key_count: u64) -> Self {
        FuseIndexHeader {
            magic: [b'D', b'B', b'F', fingerprint_bits],
            format_version: FUSE_FORMAT_VERSION,
            kmer_length,
            window_size,
            key_count,
        }
    }

    /// Fingerprint width in bits (16 or 32), encoded in the magic
    pub fn fingerprint_bits(&self) -> u8 {
        self.magic[3]
    }

    pub fn validate(&self) -> anyhow::Result<()> {
        if !self.magic.starts_with(b"DBF") {
            return Err(anyhow::anyhow!("Not a BFF index (bad magic bytes)"));
        }
        if self.format_version != FUSE_FORMAT_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported BFF format version: {}",
                self.format_version
            ));
        }
        if !matches!(self.fingerprint_bits(), 16 | 32) {
            return Err(anyhow::anyhow!(
                "Unsupported BFF fingerprint width: {} bits (expected 16 or 32)",
                self.fingerprint_bits()
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

/// Load minimizers from a reader (generic over any Read impl)
pub(crate) fn load_minimizers(reader: &mut impl Read) -> Result<(IndexStorage, IndexHeader)> {
    let config = bincode::config::standard().with_fixed_int_encoding();

    // Deserialise header
    let header: IndexHeader =
        decode_from_std_read(reader, config).context("Failed to deserialise index header")?;
    header.validate()?;

    // Deserialise the count of minimizers (stored as u64 for cross-platform compatibility)
    let count: u64 =
        decode_from_std_read(reader, config).context("Failed to deserialise minimizer count")?;
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
        IndexStorage::ExactU64(set)
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
        IndexStorage::ExactU128(set)
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

/// Load exact minimizers from a file.
#[cfg(feature = "io")]
pub(crate) fn load_minimizers_from_path(path: &Path) -> Result<(IndexStorage, IndexHeader)> {
    let file =
        std::fs::File::open(path).context(format!("Failed to open index file {:?}", path))?;
    load_minimizers(&mut BufReader::with_capacity(1 << 20, file))
}

/// Load a BFF index from a reader (synthesizes an IndexHeader from k/w)
pub(crate) fn load_fuse(reader: &mut impl Read) -> Result<(IndexStorage, IndexHeader)> {
    let config = bincode::config::standard().with_fixed_int_encoding();

    // Header via serde, filter via native bincode
    let header: FuseIndexHeader =
        decode_from_std_read(reader, config).context("Failed to deserialise BFF header")?;
    header.validate()?;

    // Decode the filter at the width recorded in the magic
    let filter = match header.fingerprint_bits() {
        16 => FuseFilterKind::Bits16(
            bincode::decode_from_std_read(reader, config)
                .context("Failed to deserialise BFF filter, index may be corrupt")?,
        ),
        32 => FuseFilterKind::Bits32(
            bincode::decode_from_std_read(reader, config)
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
    let set = IndexStorage::Fuse(FuseFilter {
        filter,
        key_count: header.key_count as usize,
    });
    Ok((set, index_header))
}

/// Load one exact/BFF index, leaving trailing bytes unread.
/// Use `BufReader` for files and streams.
pub fn load_index(reader: &mut impl Read) -> Result<crate::Index> {
    let mut magic = [0u8; 4];
    reader
        .read_exact(&mut magic)
        .context("Failed to read index header, file may be empty or truncated")?;
    // Re-chain the consumed magic bytes ahead of the stream
    let mut combined = (&magic[..]).chain(reader);
    let (minimizers, header) = if magic.starts_with(b"DBF") {
        load_fuse(&mut combined)?
    } else {
        load_minimizers(&mut combined)?
    };
    Ok(crate::Index { minimizers, header })
}

/// Load an index from a path, auto-detecting exact vs BFF format
pub fn load_index_from_path(path: &Path) -> Result<crate::Index> {
    let file =
        std::fs::File::open(path).context(format!("Failed to open index file {:?}", path))?;
    load_index(&mut BufReader::with_capacity(1 << 20, file))
}

/// Write exact minimizers to a file or stdout.
#[cfg(feature = "io")]
pub(crate) fn write_exact_index_to_path(
    minimizers: &mut IndexStorage,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    let minimizer_list = match minimizers {
        IndexStorage::ExactU64(set) => MinimizerShards::U64(vec![sort_hashset(set)]),
        IndexStorage::ExactU128(set) => MinimizerShards::U128(vec![sort_hashset(set)]),
        IndexStorage::Fuse(_) => {
            return Err(anyhow::anyhow!(
                "Cannot serialise a BFF index in the exact index format"
            ));
        }
    };
    write_exact_index_shards_to_path(minimizer_list, header, output_path)
}

/// Write an exact index from minimizer shards.
#[cfg(feature = "io")]
pub(crate) fn write_exact_index_shards_to_path(
    minimizers: MinimizerShards,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    // Create writer based on output path
    let mut writer: BufWriter<Box<dyn Write>> = match output_path {
        Some(path) if path.as_os_str() != "-" => BufWriter::new(Box::new(
            std::fs::File::create(path).context("Failed to create output file")?,
        )),
        _ => BufWriter::new(Box::new(io::stdout())),
    };

    write_exact_index_shards(minimizers, header, &mut writer)?;
    writer.flush().context("Failed to flush index")
}

fn write_exact_index_shards(
    minimizers: MinimizerShards,
    header: &IndexHeader,
    mut writer: &mut impl Write,
) -> Result<()> {
    // Use fixed-width little-endian encoding for all values.
    let config = bincode::config::standard().with_fixed_int_encoding();
    encode_into_std_write(header, &mut writer, config)
        .context("Failed to serialise index header")?;

    // Serialise the count of minimizers first (as u64 for cross-platform compatibility)
    let count = minimizers.len();
    encode_into_std_write(count, &mut writer, config)
        .context("Failed to serialise minimizer count")?;

    // Serialise minimizers in byte-aligned packed format
    let bytes_per_minimizer = (header.kmer_length as usize).div_ceil(4);
    match minimizers {
        MinimizerShards::U64(minimizers) => {
            info!("Writing minimizers");
            for list in minimizers {
                for val in list {
                    // Write only the required bytes (little-endian)
                    let bytes = val.to_le_bytes();
                    writer
                        .write_all(&bytes[..bytes_per_minimizer])
                        .context("Failed to write minimizer")?;
                }
            }
        }
        MinimizerShards::U128(minimizers) => {
            info!("Writing minimizers");
            for list in minimizers {
                for val in list {
                    // Write only the required bytes (little-endian)
                    let bytes = val.to_le_bytes();
                    writer
                        .write_all(&bytes[..bytes_per_minimizer])
                        .context("Failed to write minimizer")?;
                }
            }
        }
    }
    Ok(())
}

/// Write an exact/BFF index without modifying it.
/// The caller must flush buffers and finish compression.
pub fn write_index(index: &crate::Index, writer: &mut impl Write) -> Result<()> {
    let header = &index.header;
    match &index.minimizers {
        IndexStorage::Fuse(filter) => {
            let config = bincode::config::standard().with_fixed_int_encoding();
            let header = FuseIndexHeader::new(
                filter.fingerprint_bits(),
                header.kmer_length,
                header.window_size,
                filter.key_count as u64,
            );
            encode_into_std_write(&header, writer, config)
                .context("Failed to serialise BFF header")?;
            match &filter.filter {
                FuseFilterKind::Bits16(filter) => {
                    bincode::encode_into_std_write(filter, writer, config)
                }
                FuseFilterKind::Bits32(filter) => {
                    bincode::encode_into_std_write(filter, writer, config)
                }
            }
            .context("Failed to serialise BFF filter")?;
            Ok(())
        }
        IndexStorage::ExactU64(set) => write_exact_index_shards(
            MinimizerShards::U64(vec![sorted_values(set)]),
            header,
            writer,
        ),
        IndexStorage::ExactU128(set) => write_exact_index_shards(
            MinimizerShards::U128(vec![sorted_values(set)]),
            header,
            writer,
        ),
    }
}

/// Write an index. `-` selects stdout.
pub fn write_index_to_path(index: &crate::Index, path: &Path) -> Result<()> {
    let mut writer: BufWriter<Box<dyn Write>> = if path.as_os_str() == "-" {
        BufWriter::new(Box::new(io::stdout()))
    } else {
        BufWriter::new(Box::new(std::fs::File::create(path).with_context(
            || format!("Failed to create index: {}", path.display()),
        )?))
    };
    write_index(index, &mut writer)?;
    writer.flush().context("Failed to flush index")
}

fn sorted_values<T: Copy + std::hash::Hash + Ord>(set: &RapidHashSet<T>) -> Vec<T> {
    let buckets = table_buckets(set.len());
    let mut values: Vec<_> = set.iter().copied().collect();
    values.sort_unstable_by_key(|v| (FixedRapidHasher.hash_one(v) as usize & (buckets - 1), *v));
    values
}

/// Hashbrown's bucket count for `len` items, guarded by `test_table_buckets_matches_hashbrown`
pub(crate) fn table_buckets(len: usize) -> usize {
    match len {
        0..4 => 4,
        4..8 => 8,
        8..15 => 16,
        _ => (len * 8 / 7).next_power_of_two(),
    }
}

/// Sort by (bucket, value) for deterministic output and fast loading.
/// Shrink first, then insertion-sort the mostly ordered iteration output.
#[cfg(any(feature = "io", test))]
fn sort_hashset<T>(set: &mut RapidHashSet<T>) -> Vec<T>
where
    T: Copy + std::hash::Hash + Ord + Send + Sync,
{
    info!("Sorting minimizers");
    // Shrink so iteration follows table_buckets(len). capacity() is unreliable after removals
    set.shrink_to_fit();

    let num_buckets = table_buckets(set.len());

    let bucket = |x: &T| -> usize { FixedRapidHasher.hash_one(x) as usize & (num_buckets - 1) };

    // Get a vec of values, and do a silly but fast insertion sort on it,
    // since the vec is already mostly sorted by bucket anyway.
    let mut vals: Vec<_> = set.iter().copied().collect();

    // insertion sort
    for i in 1..vals.len() {
        let x = (bucket(&vals[i]), vals[i]);

        // where does x go? iterate backwards
        let mut j = i;
        while j > 0 {
            let y = (bucket(&vals[j - 1]), vals[j - 1]);
            if x >= y {
                break;
            }
            j -= 1;
        }

        // rotate x to the right slot.
        if j != i {
            vals[j..=i].rotate_right(1);
        }
    }

    info!("Checking that values are sorted by hash");
    assert!(vals.is_sorted_by_key(|x| (bucket(x), x)));

    vals
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_creation() {
        let header = IndexHeader::new(31, 21);

        assert_eq!(header.format_version, 3);
        assert_eq!(header.kmer_length, 31);
        assert_eq!(header.window_size, 21);
    }

    /// Check hashbrown hasn't changed its table sizing which would break stuff
    #[test]
    fn test_table_buckets_matches_hashbrown() {
        for len in (1..2000).chain([1 << 16, 100_000]) {
            let buckets = table_buckets(len);
            let capacity = if buckets <= 8 {
                buckets - 1
            } else {
                buckets / 8 * 7
            };
            let set = RapidHashSet::<u64>::with_capacity_and_hasher(len, FixedRapidHasher);
            assert_eq!(set.capacity(), capacity, "u64 len={len}");
            let set = RapidHashSet::<u128>::with_capacity_and_hasher(len, FixedRapidHasher);
            assert_eq!(set.capacity(), capacity, "u128 len={len}");
        }
    }

    /// Tiny sets and sets with removal tombstones must not panic
    #[test]
    fn test_sort_hashset_small_and_after_removals() {
        for len in [0, 1, 2, 3, 4, 5, 7, 8, 14, 15, 1000, 100_000] {
            let mut set: RapidHashSet<u64> = (0..len + len / 10).collect();
            set.retain(|&x| x < len);
            let vals = sort_hashset(&mut set);
            assert_eq!(vals.len(), len as usize);
        }
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
        assert!(FuseIndexHeader::new(16, 31, 15, 100).validate().is_ok());
        assert!(FuseIndexHeader::new(32, 31, 15, 100).validate().is_ok());

        // Bad magic
        let mut bad_magic = FuseIndexHeader::new(16, 31, 15, 100);
        bad_magic.magic = *b"XXXX";
        assert!(bad_magic.validate().is_err());

        // Bad version
        let mut bad_ver = FuseIndexHeader::new(16, 31, 15, 100);
        bad_ver.format_version = 99;
        assert!(bad_ver.validate().is_err());

        // Unsupported filter width
        let mut bad_bits = FuseIndexHeader::new(16, 31, 15, 100);
        bad_bits.magic[3] = 8;
        assert!(bad_bits.validate().is_err());

        // k > 32 is rejected (BFF keys on raw u64 minimizers)
        assert!(FuseIndexHeader::new(16, 41, 21, 100).validate().is_err());
    }

    #[rstest::rstest]
    #[case(16)]
    #[case(32)]
    fn test_bff_roundtrip_membership(#[case] bits: u8) {
        let keys: Vec<u64> = (0..10_000u64)
            .map(|i| i.wrapping_mul(0x9E3779B97F4A7C15))
            .collect();
        let header = FuseIndexHeader::new(bits, 31, 15, keys.len() as u64);

        let config = bincode::config::standard().with_fixed_int_encoding();
        let mut buf = Vec::new();
        encode_into_std_write(&header, &mut buf, config).unwrap();
        match bits {
            16 => bincode::encode_into_std_write(
                xorf::BinaryFuse16::try_from(&keys).unwrap(),
                &mut buf,
                config,
            ),
            _ => bincode::encode_into_std_write(
                xorf::BinaryFuse32::try_from(&keys).unwrap(),
                &mut buf,
                config,
            ),
        }
        .unwrap();

        let mut cursor = std::io::Cursor::new(&buf);
        let (set, idx_header) = load_fuse(&mut cursor).unwrap();
        assert_eq!(idx_header.kmer_length, 31);
        assert_eq!(idx_header.window_size, 15);
        assert_eq!(set.len(), keys.len());
        let IndexStorage::Fuse(filter) = set else {
            panic!("Expected fuse index")
        };
        // No false negatives
        for &k in &keys {
            assert!(filter.contains(k));
        }

        // auto-detect dispatches to BFF
        let mut cursor2 = std::io::Cursor::new([buf.as_slice(), buf.as_slice()].concat());
        let index = load_index(&mut cursor2).unwrap();
        assert_eq!(cursor2.position(), buf.len() as u64);
        assert_eq!(load_index(&mut cursor2).unwrap().len(), keys.len());
        assert_eq!(
            index.kind(),
            crate::IndexKind::Fuse {
                fingerprint_bits: bits
            }
        );
        let mut encoded = Vec::new();
        write_index(&index, &mut encoded).unwrap();
        assert_eq!(encoded, buf);
        let mut kernel =
            crate::FilterKernel::new(std::sync::Arc::new(index), crate::FilterParams::default())
                .unwrap();
        let read = b"ACGTTGCAAGGCTTAACCGGTTACGATCGATCGGATCCTAGCTAGCTTAACCGGATCGTA";
        let counts = kernel.score_read(read);
        let diagnostics = kernel.classify_read_with_diagnostics(read);
        assert_eq!(counts, diagnostics.score);
    }
}
