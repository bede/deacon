use std::io::{Cursor, Write};
use std::sync::Arc;

use flate2::Compression;
use flate2::write::GzEncoder;
use js_sys::{Object, Reflect, Uint8Array};
use wasm_bindgen::prelude::*;

use deacon::minimizers::{Buffers, KmerHasher};
use deacon::{MinimizerSet, MinimizerVec, RapidHashSet};
use needletail::parse_fastx_reader;

struct WasmIndexInner {
    minimizers: MinimizerSet,
    header: deacon::IndexHeader,
}

#[wasm_bindgen]
pub struct WasmIndex {
    inner: Arc<WasmIndexInner>,
}

#[wasm_bindgen]
impl WasmIndex {
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<WasmIndex, JsValue> {
        console_error_panic_hook::set_once();
        let mut cursor = Cursor::new(data);
        let (minimizers, header) = deacon::load_minimizers(&mut cursor)
            .map_err(|e| JsValue::from_str(&format!("Failed to load index: {}", e)))?;
        Ok(WasmIndex {
            inner: Arc::new(WasmIndexInner { minimizers, header }),
        })
    }

    pub fn info(&self) -> String {
        format!(
            "k={}, w={}, {} minimizers",
            self.inner.header.kmer_length(),
            self.inner.header.window_size(),
            self.inner.minimizers.len()
        )
    }
}

#[wasm_bindgen]
pub struct WasmFilterSession {
    index: Arc<WasmIndexInner>,
    deplete: bool,
    abs_threshold: usize,
    rel_threshold: f64,
    hasher: KmerHasher,
    buffers: Buffers,
    parser: FastqChunkParser,
    stats: FilterStats,
}

#[wasm_bindgen]
impl WasmFilterSession {
    #[wasm_bindgen(constructor)]
    pub fn new(
        index: &WasmIndex,
        deplete: bool,
        abs_threshold: usize,
        rel_threshold: f64,
    ) -> WasmFilterSession {
        let k = index.inner.header.kmer_length();
        let hasher = KmerHasher::new(k as usize);
        let buffers = if k <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        WasmFilterSession {
            index: Arc::clone(&index.inner),
            deplete,
            abs_threshold,
            rel_threshold,
            hasher,
            buffers,
            parser: FastqChunkParser::new(),
            stats: FilterStats::default(),
        }
    }

    pub fn push_chunk(&mut self, chunk: &[u8]) -> Result<Vec<u8>, JsValue> {
        let mut output = Vec::new();
        let records = self
            .parser
            .push_chunk(chunk)
            .map_err(|e| JsValue::from_str(&e))?;
        for record in records {
            self.process_record(&record.header, &record.seq, &record.qual, &mut output)?;
        }
        Ok(output)
    }

    pub fn finish(&mut self) -> Result<Vec<u8>, JsValue> {
        let mut output = Vec::new();
        let records = self.parser.finish().map_err(|e| JsValue::from_str(&e))?;
        for record in records {
            self.process_record(&record.header, &record.seq, &record.qual, &mut output)?;
        }
        Ok(output)
    }

    pub fn stats(&self) -> Result<JsValue, JsValue> {
        let stats_obj = Object::new();
        set_field(&stats_obj, "readsIn", self.stats.reads_in)?;
        set_field(&stats_obj, "readsOut", self.stats.reads_out)?;
        set_field(&stats_obj, "basesIn", self.stats.bases_in)?;
        set_field(&stats_obj, "basesOut", self.stats.bases_out)?;
        Ok(stats_obj.into())
    }

    pub fn pending_bytes(&self) -> usize {
        self.parser.pending_bytes()
    }

    fn process_record(
        &mut self,
        header: &[u8],
        seq: &[u8],
        qual: &[u8],
        output: &mut Vec<u8>,
    ) -> Result<(), JsValue> {
        let k = self.index.header.kmer_length();
        let w = self.index.header.window_size();

        self.stats.reads_in += 1;
        self.stats.bases_in += seq.len() as u64;

        let keep = if seq.len() < k as usize {
            self.deplete
        } else {
            deacon::minimizers::fill_minimizers(seq, &self.hasher, k, w, 0.0, &mut self.buffers);
            let num_minimizers = self.buffers.minimizers.len();
            let hit_count = count_hits(&self.buffers.minimizers, &self.index.minimizers);
            let rel_required = if num_minimizers == 0 {
                0
            } else {
                ((self.rel_threshold * num_minimizers as f64).round() as usize).max(1)
            };
            let required = self.abs_threshold.max(rel_required);
            if self.deplete {
                hit_count < required
            } else {
                hit_count >= required
            }
        };

        if keep {
            self.stats.reads_out += 1;
            self.stats.bases_out += seq.len() as u64;
            write_fastq_record(output, header, seq, qual)
                .map_err(|e| JsValue::from_str(&format!("Output write error: {}", e)))?;
        }

        Ok(())
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

struct FastqChunkParser {
    pending: Vec<u8>,
    record_lines: Vec<Vec<u8>>, // always 0..=3 lines
}

struct FastqRecord {
    header: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl FastqChunkParser {
    fn new() -> Self {
        Self {
            pending: Vec::new(),
            record_lines: Vec::new(),
        }
    }

    fn pending_bytes(&self) -> usize {
        self.pending.len()
    }

    fn push_chunk(&mut self, chunk: &[u8]) -> Result<Vec<FastqRecord>, String> {
        if chunk.is_empty() {
            return Ok(Vec::new());
        }

        let mut records = Vec::new();
        self.pending.extend_from_slice(chunk);
        let mut cursor = 0usize;

        while let Some(rel_pos) = memchr::memchr(b'\n', &self.pending[cursor..]) {
            let line_end = cursor + rel_pos;
            let mut line = self.pending[cursor..line_end].to_vec();
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            if let Some(record) = self.push_line(line)? {
                records.push(record);
            }
            cursor = line_end + 1;
        }

        if cursor > 0 {
            self.pending.drain(..cursor);
        }

        Ok(records)
    }

    fn finish(&mut self) -> Result<Vec<FastqRecord>, String> {
        let mut records = Vec::new();
        if !self.pending.is_empty() {
            let mut line = std::mem::take(&mut self.pending);
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            if let Some(record) = self.push_line(line)? {
                records.push(record);
            }
        }

        if !self.record_lines.is_empty() {
            return Err(
                "Incomplete FASTQ record at end of stream (expected 4 lines per record)"
                    .to_string(),
            );
        }

        Ok(records)
    }

    fn push_line(&mut self, line: Vec<u8>) -> Result<Option<FastqRecord>, String> {
        if self.record_lines.is_empty() && line.is_empty() {
            // Skip blank lines between records.
            return Ok(None);
        }

        self.record_lines.push(line);

        if self.record_lines.len() < 4 {
            return Ok(None);
        }

        if self.record_lines[0].first() != Some(&b'@') {
            return Err(
                "Streaming mode currently supports FASTQ records only (header must start with '@')"
                    .to_string(),
            );
        }
        if self.record_lines[2].first() != Some(&b'+') {
            return Err("Invalid FASTQ record: third line must start with '+'".to_string());
        }

        let header = self.record_lines[0][1..].to_vec();
        let seq = self.record_lines[1].clone();
        let qual = self.record_lines[3].clone();

        if qual.len() != seq.len() {
            return Err("Invalid FASTQ record: sequence and quality lengths differ".to_string());
        }

        self.record_lines.clear();
        Ok(Some(FastqRecord { header, seq, qual }))
    }
}

fn run_filter(
    index: &WasmIndex,
    input: &[u8],
    deplete: bool,
    abs_threshold: usize,
    rel_threshold: f64,
) -> Result<(Vec<u8>, FilterStats), JsValue> {
    let k = index.inner.header.kmer_length();
    let w = index.inner.header.window_size();
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
            let hit_count = count_hits(&buffers.minimizers, &index.inner.minimizers);
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

fn write_fastq_record(
    output: &mut impl Write,
    header: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> std::io::Result<()> {
    output.write_all(b"@")?;
    output.write_all(header)?;
    output.write_all(b"\n")?;
    output.write_all(seq)?;
    output.write_all(b"\n+\n")?;
    output.write_all(qual)?;
    output.write_all(b"\n")?;
    Ok(())
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

#[cfg(test)]
mod tests {
    use super::FastqChunkParser;

    #[test]
    fn parser_handles_byte_sized_chunks() {
        let input = b"@r1\nACGT\n+\n!!!!\n@r2\nTGCA\n+\n####\n";
        let mut parser = FastqChunkParser::new();
        let mut records = Vec::new();

        for chunk in input.chunks(1) {
            let mut out = parser.push_chunk(chunk).unwrap();
            records.append(&mut out);
        }
        let mut tail = parser.finish().unwrap();
        records.append(&mut tail);

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, b"r1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[0].qual, b"!!!!");
        assert_eq!(records[1].header, b"r2");
        assert_eq!(records[1].seq, b"TGCA");
        assert_eq!(records[1].qual, b"####");
    }

    #[test]
    fn parser_accepts_final_line_without_trailing_newline() {
        let input = b"@r1\nACGT\n+\n!!!!";
        let mut parser = FastqChunkParser::new();
        parser.push_chunk(input).unwrap();
        let records = parser.finish().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, b"r1");
    }

    #[test]
    fn parser_rejects_incomplete_record() {
        let mut parser = FastqChunkParser::new();
        parser.push_chunk(b"@r1\nACGT\n+\n").unwrap();
        assert!(parser.finish().is_err());
    }

    #[test]
    fn parser_rejects_mismatched_quality_length() {
        let mut parser = FastqChunkParser::new();
        let err = match parser.push_chunk(b"@r1\nACGT\n+\n!!!\n") {
            Ok(_) => panic!("expected parser error for mismatched FASTQ quality length"),
            Err(err) => err,
        };
        assert!(err.contains("sequence and quality lengths differ"));
    }
}
