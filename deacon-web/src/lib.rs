use std::io::{Cursor, Write};
use flate2::write::GzEncoder;
use flate2::Compression;
use js_sys::{Object, Reflect, Uint8Array};
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
    let (output, _) = run_filter(index, input, deplete, abs_threshold, rel_threshold)?;
    Ok(output)
}

#[wasm_bindgen]
pub fn filter_with_stats(
    index: &WasmIndex,
    input: &[u8],
    deplete: bool,
    abs_threshold: usize,
    rel_threshold: f64,
) -> Result<JsValue, JsValue> {
    let (output, stats) = run_filter(index, input, deplete, abs_threshold, rel_threshold)?;
    let output_bytes = Uint8Array::from(output.as_slice());

    let stats_obj = Object::new();
    set_field(&stats_obj, "readsIn", stats.reads_in)?;
    set_field(&stats_obj, "readsOut", stats.reads_out)?;
    set_field(&stats_obj, "basesIn", stats.bases_in)?;
    set_field(&stats_obj, "basesOut", stats.bases_out)?;

    let result_obj = Object::new();
    Reflect::set(&result_obj, &"output".into(), &output_bytes.into())?;
    Reflect::set(&result_obj, &"stats".into(), &stats_obj.into())?;

    Ok(result_obj.into())
}

#[derive(Default)]
struct FilterStats {
    reads_in: u64,
    reads_out: u64,
    bases_in: u64,
    bases_out: u64,
}

fn run_filter(
    index: &WasmIndex,
    input: &[u8],
    deplete: bool,
    abs_threshold: usize,
    rel_threshold: f64,
) -> Result<(Vec<u8>, FilterStats), JsValue> {
    let k = index.header.kmer_length();
    let w = index.header.window_size();
    let hasher = KmerHasher::new(k as usize);
    let mut buffers = if k <= 32 {
        Buffers::new_u64()
    } else {
        Buffers::new_u128()
    };

    // Avoid preallocating output to full input size, which can panic/abort in WASM
    // for large inputs. Let the output buffer grow incrementally instead.
    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(2));
    let mut stats = FilterStats::default();

    let mut reader = parse_fastx_reader(Cursor::new(input))
        .map_err(|e| JsValue::from_str(&format!("Failed to parse input: {}", e)))?;

    while let Some(result) = reader.next() {
        let record = result.map_err(|e| JsValue::from_str(&format!("Parse error: {}", e)))?;
        stats.reads_in += 1;
        let seq = record.seq();
        stats.bases_in += seq.len() as u64;

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
            stats.reads_out += 1;
            stats.bases_out += seq.len() as u64;
            write_record(&record, &seq, &mut encoder);
        }
    }

    let output = encoder
        .finish()
        .map_err(|e| JsValue::from_str(&format!("Compression error: {}", e)))?;
    Ok((output, stats))
}

fn set_field(obj: &Object, key: &str, value: u64) -> Result<(), JsValue> {
    Reflect::set(obj, &key.into(), &JsValue::from_f64(value as f64))?;
    Ok(())
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
