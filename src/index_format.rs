use crate::{FixedRapidHasher, MinimizerSet, MinimizerVecVec, RapidHashSet};
use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use serde::{Deserialize, Serialize};
use std::hash::{BuildHasher, Hasher};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use tracing::info;

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
    #[cfg(any(feature = "io", test))]
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

/// Helper function to write minimizers to output file or stdout
pub fn dump_minimizers(
    minimizers: &mut crate::MinimizerSet,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    let minimizer_list = match minimizers {
        MinimizerSet::U64(set) => MinimizerVecVec::U64(vec![sort_hashset(set)]),
        MinimizerSet::U128(set) => MinimizerVecVec::U128(vec![sort_hashset(set)]),
        MinimizerSet::Fuse(_) => {
            return Err(anyhow::anyhow!(
                "Cannot serialise a BFF index in the exact index format"
            ));
        }
    };
    dump_minimizer_lists(minimizer_list, header, output_path)
}

/// Write an index for the given minimizers.
/// Takes a list of lists so it can be used with multiple shards.
pub(crate) fn dump_minimizer_lists(
    minimizers: MinimizerVecVec,
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
        MinimizerVecVec::U64(minimizers) => {
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
        MinimizerVecVec::U128(minimizers) => {
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

/// Hashbrown's bucket count for `len` items, guarded by `test_table_buckets_matches_hashbrown`
pub(crate) fn table_buckets(len: usize) -> usize {
    match len {
        0..4 => 4,
        4..8 => 8,
        8..15 => 16,
        _ => (len * 8 / 7).next_power_of_two(),
    }
}

/// Sort the values in the hashset by (bucket, value).
///
/// This way, the output is deterministic, and construction from an index file
/// is fast because the data is already in the right order.
///
/// This works by first shrinking the input hashset to fit so it has the correct final size.
/// Then, we collect all elements to a vector by iteration order,
/// which _mostly_ but not exactly returns them by order of target bucket.
/// We end with a naive insertion sort to precisely sort all values by their target bucket.
pub(crate) fn sort_hashset<T>(set: &mut RapidHashSet<T>) -> Vec<T>
where
    T: Copy + std::hash::Hash + Ord + Send + Sync,
{
    info!("Sorting minimizers");
    // Shrink so iteration follows table_buckets(len). capacity() is unreliable after removals
    set.shrink_to_fit();

    let num_buckets = table_buckets(set.len());

    let bucket = |x: &T| -> usize {
        let mut hasher = FixedRapidHasher::default().build_hasher();
        x.hash(&mut hasher);
        let hash = hasher.finish() as usize;
        hash & (num_buckets - 1)
    };

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
        assert_eq!(header.kmer_length(), 31);
        assert_eq!(header.window_size(), 21);
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
