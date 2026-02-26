use std::io::Cursor;
use std::sync::OnceLock;
use wasm_bindgen::prelude::*;

use deacon::minimizers::{Buffers, KmerHasher};
use deacon::{MinimizerSet, MinimizerVec, RapidHashSet};

struct LoadedIndex {
    minimizers: MinimizerSet,
    header: deacon::IndexHeader,
}

static INDEX: OnceLock<LoadedIndex> = OnceLock::new();

#[wasm_bindgen]
pub fn load_index(data: &[u8]) -> Result<(), JsValue> {
    let mut cursor = Cursor::new(data);
    let (minimizers, header) = deacon::load_minimizers(&mut cursor)
        .map_err(|e| JsValue::from_str(&format!("Failed to load index: {}", e)))?;

    INDEX
        .set(LoadedIndex { minimizers, header })
        .map_err(|_| JsValue::from_str("Index already loaded"))?;

    let idx = INDEX.get().unwrap();
    web_sys::console::log_1(
        &format!(
            "Loaded index: k={}, w={}, {} minimizers",
            idx.header.kmer_length(),
            idx.header.window_size(),
            idx.minimizers.len()
        )
        .into(),
    );
    Ok(())
}

#[wasm_bindgen]
pub fn get_index_info() -> Result<String, JsValue> {
    let idx = INDEX
        .get()
        .ok_or_else(|| JsValue::from_str("No index loaded"))?;
    Ok(format!(
        "k={}, w={}, {} minimizers",
        idx.header.kmer_length(),
        idx.header.window_size(),
        idx.minimizers.len()
    ))
}

/// Filter FASTQ data. Input may be gzip-compressed.
/// Returns filtered FASTQ as bytes.
#[wasm_bindgen]
pub fn filter_fastq(
    input: &[u8],
    is_gzipped: bool,
    deplete: bool,
    abs_threshold: usize,
    rel_threshold: f64,
) -> Result<Vec<u8>, JsValue> {
    let idx = INDEX
        .get()
        .ok_or_else(|| JsValue::from_str("No index loaded. Call load_index first."))?;

    let decompressed: Vec<u8>;
    let data: &[u8] = if is_gzipped {
        use flate2::read::GzDecoder;
        use std::io::Read;
        let mut decoder = GzDecoder::new(input);
        let mut buf = Vec::new();
        decoder
            .read_to_end(&mut buf)
            .map_err(|e| JsValue::from_str(&format!("Gzip decompression failed: {}", e)))?;
        decompressed = buf;
        &decompressed
    } else {
        input
    };

    let k = idx.header.kmer_length();
    let w = idx.header.window_size();
    let hasher = KmerHasher::new(k as usize);
    let mut buffers = if k <= 32 {
        Buffers::new_u64()
    } else {
        Buffers::new_u128()
    };

    let mut output = Vec::with_capacity(data.len());
    let mut pos = 0;
    let mut seqs_in: u64 = 0;
    let mut seqs_out: u64 = 0;

    // Simple 4-line FASTQ parser
    while pos < data.len() {
        // Line 1: header (starts with @)
        let header_start = pos;
        let header_end = memchr_newline(data, pos);
        if header_end >= data.len() {
            break;
        }

        // Line 2: sequence
        let seq_start = header_end + 1;
        let seq_end = memchr_newline(data, seq_start);
        if seq_end >= data.len() {
            break;
        }

        // Line 3: plus line
        let plus_start = seq_end + 1;
        let plus_end = memchr_newline(data, plus_start);
        if plus_end >= data.len() {
            break;
        }

        // Line 4: quality
        let qual_start = plus_end + 1;
        let qual_end = memchr_newline(data, qual_start);

        let seq = &data[seq_start..seq_end];
        seqs_in += 1;

        // Compute minimizers and count hits
        let keep = should_keep_sequence(
            seq,
            &hasher,
            k,
            w,
            &idx.minimizers,
            &mut buffers,
            abs_threshold,
            rel_threshold,
            deplete,
        );

        if keep {
            seqs_out += 1;
            // Write entire record (header + seq + plus + qual + newlines)
            let record_end = if qual_end < data.len() {
                qual_end + 1
            } else {
                qual_end
            };
            output.extend_from_slice(&data[header_start..record_end]);
        }

        pos = if qual_end < data.len() {
            qual_end + 1
        } else {
            data.len()
        };
    }

    web_sys::console::log_1(
        &format!(
            "Filtered: {}/{} sequences retained ({}%)",
            seqs_out,
            seqs_in,
            if seqs_in > 0 {
                (seqs_out as f64 / seqs_in as f64 * 100.0) as u64
            } else {
                0
            }
        )
        .into(),
    );

    Ok(output)
}

/// Find next newline position from start, returns data.len() if not found
fn memchr_newline(data: &[u8], start: usize) -> usize {
    for i in start..data.len() {
        if data[i] == b'\n' {
            return i;
        }
    }
    data.len()
}

/// Check if a sequence should be kept based on minimizer hits
fn should_keep_sequence(
    seq: &[u8],
    hasher: &KmerHasher,
    k: u8,
    w: u8,
    minimizers_set: &MinimizerSet,
    buffers: &mut Buffers,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
) -> bool {
    if seq.len() < k as usize {
        return deplete; // Too short: keep in deplete mode, skip in search mode
    }

    // Fill buffers with minimizers
    deacon::minimizers::fill_minimizers(seq, hasher, k, w, 0.0, buffers);

    let num_minimizers = buffers.minimizers.len();

    // Count distinct hits
    let hit_count = match (&buffers.minimizers, minimizers_set) {
        (MinimizerVec::U64(vec), MinimizerSet::U64(set)) => {
            let mut seen = RapidHashSet::default();
            for &m in vec {
                if set.contains(&m) {
                    seen.insert(m);
                }
            }
            seen.len()
        }
        (MinimizerVec::U128(vec), MinimizerSet::U128(set)) => {
            let mut seen = RapidHashSet::default();
            for &m in vec {
                if set.contains(&m) {
                    seen.insert(m);
                }
            }
            seen.len()
        }
        _ => 0,
    };

    // Calculate required hits (same logic as native deacon)
    let abs_required = abs_threshold;
    let rel_required = if num_minimizers == 0 {
        0
    } else {
        ((rel_threshold * num_minimizers as f64).round() as usize).max(1)
    };
    let required = abs_required.max(rel_required);

    if deplete {
        hit_count < required
    } else {
        hit_count >= required
    }
}
