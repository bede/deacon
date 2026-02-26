use std::io::{Cursor, Write};
use flate2::write::GzEncoder;
use flate2::Compression;
use wasm_bindgen::prelude::*;

use deacon::minimizers::{Buffers, KmerHasher};
use deacon::{MinimizerSet, MinimizerVec, RapidHashSet};
use needletail::parse_fastx_reader;

#[wasm_bindgen]
pub struct WasmIndex {
    minimizers: MinimizerSet,
    header: deacon::IndexHeader,
}

#[wasm_bindgen]
impl WasmIndex {
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<WasmIndex, JsValue> {
        console_error_panic_hook::set_once();
        let mut cursor = Cursor::new(data);
        let (minimizers, header) = deacon::load_minimizers(&mut cursor)
            .map_err(|e| JsValue::from_str(&format!("Failed to load index: {}", e)))?;
        Ok(WasmIndex { minimizers, header })
    }

    pub fn info(&self) -> String {
        format!(
            "k={}, w={}, {} minimizers",
            self.header.kmer_length(),
            self.header.window_size(),
            self.minimizers.len()
        )
    }
}

#[wasm_bindgen]
pub fn filter(
    index: &WasmIndex,
    input: &[u8],
    deplete: bool,
    abs_threshold: usize,
    rel_threshold: f64,
) -> Result<Vec<u8>, JsValue> {
    let k = index.header.kmer_length();
    let w = index.header.window_size();
    let hasher = KmerHasher::new(k as usize);
    let mut buffers = if k <= 32 {
        Buffers::new_u64()
    } else {
        Buffers::new_u128()
    };

    let mut encoder = GzEncoder::new(Vec::with_capacity(input.len()), Compression::new(2));
    let mut seqs_in: u64 = 0;
    let mut seqs_out: u64 = 0;

    let mut reader = parse_fastx_reader(Cursor::new(input))
        .map_err(|e| JsValue::from_str(&format!("Failed to parse input: {}", e)))?;

    while let Some(result) = reader.next() {
        let record = result.map_err(|e| JsValue::from_str(&format!("Parse error: {}", e)))?;
        seqs_in += 1;
        let seq = record.seq();

        let keep = if (seq.len()) < k as usize {
            deplete
        } else {
            deacon::minimizers::fill_minimizers(&seq, &hasher, k, w, 0.0, &mut buffers);
            let num_minimizers = buffers.minimizers.len();
            let hit_count = count_hits(&buffers.minimizers, &index.minimizers);
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
        };

        if keep {
            seqs_out += 1;
            write_record(&record, &seq, &mut encoder);
        }
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

    encoder.finish().map_err(|e| JsValue::from_str(&format!("Compression error: {}", e)))
}

fn count_hits(minimizers: &MinimizerVec, set: &MinimizerSet) -> usize {
    match (minimizers, set) {
        (MinimizerVec::U64(vec), MinimizerSet::U64(s)) => {
            let mut seen = RapidHashSet::default();
            for &m in vec {
                if s.contains(&m) {
                    seen.insert(m);
                }
            }
            seen.len()
        }
        (MinimizerVec::U128(vec), MinimizerSet::U128(s)) => {
            let mut seen = RapidHashSet::default();
            for &m in vec {
                if s.contains(&m) {
                    seen.insert(m);
                }
            }
            seen.len()
        }
        _ => 0,
    }
}

fn write_record(record: &needletail::parser::SequenceRecord, seq: &[u8], output: &mut impl Write) {
    if let Some(qual) = record.qual() {
        let _ = output.write_all(b"@");
        let _ = output.write_all(record.id());
        let _ = output.write_all(b"\n");
        let _ = output.write_all(seq);
        let _ = output.write_all(b"\n+\n");
        let _ = output.write_all(qual);
        let _ = output.write_all(b"\n");
    } else {
        let _ = output.write_all(b">");
        let _ = output.write_all(record.id());
        let _ = output.write_all(b"\n");
        let _ = output.write_all(seq);
        let _ = output.write_all(b"\n");
    }
}
