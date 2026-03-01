use std::io::Cursor;
use std::io::Write;
use std::sync::Arc;

use flate2::Compression;
use flate2::write::{GzEncoder, MultiGzDecoder};
use js_sys::{Object, Reflect};
use wasm_bindgen::prelude::*;

use deacon::minimizers::{Buffers, KmerHasher};
use deacon::{MinimizerSet, MinimizerVec, RapidHashSet};

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
            fmt_commas(self.inner.minimizers.len())
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
    parser: SeqChunkParser,
    stats: FilterStats,
    gz_decoder: Option<MultiGzDecoder<Vec<u8>>>,
    gz_encoder: Option<GzEncoder<Vec<u8>>>,
}

#[wasm_bindgen]
impl WasmFilterSession {
    #[wasm_bindgen(constructor)]
    pub fn new(
        index: &WasmIndex,
        deplete: bool,
        abs_threshold: usize,
        rel_threshold: f64,
        decompress_input: bool,
        compress_output: bool,
    ) -> WasmFilterSession {
        let k = index.inner.header.kmer_length();
        let hasher = KmerHasher::new(k as usize);
        let buffers = if k <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        let gz_decoder = if decompress_input {
            Some(MultiGzDecoder::new(Vec::new()))
        } else {
            None
        };
        let gz_encoder = if compress_output {
            Some(GzEncoder::new(Vec::new(), Compression::new(2)))
        } else {
            None
        };

        WasmFilterSession {
            index: Arc::clone(&index.inner),
            deplete,
            abs_threshold,
            rel_threshold,
            hasher,
            buffers,
            parser: SeqChunkParser::new(),
            stats: FilterStats::default(),
            gz_decoder,
            gz_encoder,
        }
    }

    pub fn push_chunk(&mut self, chunk: &[u8]) -> Result<Vec<u8>, JsValue> {
        // Step 1: Decompress if needed
        let input = if let Some(ref mut decoder) = self.gz_decoder {
            decoder
                .write_all(chunk)
                .map_err(|e| JsValue::from_str(&format!("Decompression error: {}", e)))?;
            std::mem::take(decoder.get_mut())
        } else {
            chunk.to_vec()
        };

        // Step 2: Parse FASTQ records and filter
        let mut raw_output = Vec::new();
        let records = self
            .parser
            .push_chunk(&input)
            .map_err(|e| JsValue::from_str(&e))?;
        for record in records {
            self.process_record(&record.header, &record.seq, record.qual.as_deref(), &mut raw_output)?;
        }

        // Step 3: Compress output if needed
        if let Some(ref mut encoder) = self.gz_encoder {
            if !raw_output.is_empty() {
                encoder
                    .write_all(&raw_output)
                    .map_err(|e| JsValue::from_str(&format!("Compression error: {}", e)))?;
            }
            Ok(std::mem::take(encoder.get_mut()))
        } else {
            Ok(raw_output)
        }
    }

    pub fn finish(&mut self) -> Result<Vec<u8>, JsValue> {
        let mut raw_output = Vec::new();

        // Finish decompression and process any remaining bytes
        if let Some(decoder) = self.gz_decoder.take() {
            let remaining = decoder
                .finish()
                .map_err(|e| JsValue::from_str(&format!("Decompression finish error: {}", e)))?;
            if !remaining.is_empty() {
                let records = self
                    .parser
                    .push_chunk(&remaining)
                    .map_err(|e| JsValue::from_str(&e))?;
                for record in records {
                    self.process_record(
                        &record.header,
                        &record.seq,
                        record.qual.as_deref(),
                        &mut raw_output,
                    )?;
                }
            }
        }

        // Finish parsing (flush any pending partial record)
        let records = self.parser.finish().map_err(|e| JsValue::from_str(&e))?;
        for record in records {
            self.process_record(&record.header, &record.seq, record.qual.as_deref(), &mut raw_output)?;
        }

        // Finish compression
        if let Some(mut encoder) = self.gz_encoder.take() {
            if !raw_output.is_empty() {
                encoder
                    .write_all(&raw_output)
                    .map_err(|e| JsValue::from_str(&format!("Compression error: {}", e)))?;
            }
            encoder
                .finish()
                .map_err(|e| JsValue::from_str(&format!("Compression finish error: {}", e)))
        } else {
            Ok(raw_output)
        }
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

    pub fn output_compressed(&self) -> bool {
        self.gz_encoder.is_some()
    }

    fn process_record(
        &mut self,
        header: &[u8],
        seq: &[u8],
        qual: Option<&[u8]>,
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
            write_record(output, header, seq, qual)
                .map_err(|e| JsValue::from_str(&format!("Output write error: {}", e)))?;
        }

        Ok(())
    }
}

#[derive(Default)]
struct FilterStats {
    reads_in: u64,
    reads_out: u64,
    bases_in: u64,
    bases_out: u64,
}

struct SeqRecord {
    header: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>, // None for FASTA
}

#[derive(Clone, Copy, PartialEq)]
enum SeqFormat {
    Fasta,
    Fastq,
}

struct SeqChunkParser {
    pending: Vec<u8>,
    format: Option<SeqFormat>,
    // FASTQ state
    record_lines: Vec<Vec<u8>>, // always 0..=3 lines
    // FASTA state
    fasta_header: Option<Vec<u8>>,
    fasta_seq: Vec<u8>,
}

impl SeqChunkParser {
    fn new() -> Self {
        Self {
            pending: Vec::new(),
            format: None,
            record_lines: Vec::new(),
            fasta_header: None,
            fasta_seq: Vec::new(),
        }
    }

    fn pending_bytes(&self) -> usize {
        self.pending.len()
    }

    fn push_chunk(&mut self, chunk: &[u8]) -> Result<Vec<SeqRecord>, String> {
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

    fn finish(&mut self) -> Result<Vec<SeqRecord>, String> {
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

        match self.format.unwrap_or(SeqFormat::Fastq) {
            SeqFormat::Fastq => {
                if !self.record_lines.is_empty() {
                    return Err(
                        "Incomplete FASTQ record at end of stream (expected 4 lines per record)"
                            .to_string(),
                    );
                }
            }
            SeqFormat::Fasta => {
                if let Some(header) = self.fasta_header.take() {
                    let seq = std::mem::take(&mut self.fasta_seq);
                    if !seq.is_empty() {
                        records.push(SeqRecord {
                            header,
                            seq,
                            qual: None,
                        });
                    }
                }
            }
        }

        Ok(records)
    }

    fn push_line(&mut self, line: Vec<u8>) -> Result<Option<SeqRecord>, String> {
        // Auto-detect format from the first non-empty line
        if self.format.is_none() {
            if line.is_empty() {
                return Ok(None);
            }
            match line.first() {
                Some(b'@') => self.format = Some(SeqFormat::Fastq),
                Some(b'>') => self.format = Some(SeqFormat::Fasta),
                _ => {
                    return Err(
                        "Cannot detect format: first record must start with '@' (FASTQ) or '>' (FASTA)"
                            .to_string(),
                    )
                }
            }
        }

        match self.format.unwrap() {
            SeqFormat::Fastq => self.push_line_fastq(line),
            SeqFormat::Fasta => self.push_line_fasta(line),
        }
    }

    fn push_line_fastq(&mut self, line: Vec<u8>) -> Result<Option<SeqRecord>, String> {
        if self.record_lines.is_empty() && line.is_empty() {
            return Ok(None);
        }

        self.record_lines.push(line);

        if self.record_lines.len() < 4 {
            return Ok(None);
        }

        if self.record_lines[0].first() != Some(&b'@') {
            return Err("Invalid FASTQ record: header must start with '@'".to_string());
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
        Ok(Some(SeqRecord {
            header,
            seq,
            qual: Some(qual),
        }))
    }

    fn push_line_fasta(&mut self, line: Vec<u8>) -> Result<Option<SeqRecord>, String> {
        if line.is_empty() {
            return Ok(None);
        }

        if line.first() == Some(&b'>') {
            // New header — emit previous record if any
            let prev = if let Some(prev_header) = self.fasta_header.take() {
                let seq = std::mem::take(&mut self.fasta_seq);
                if seq.is_empty() {
                    return Err("Empty FASTA sequence".to_string());
                }
                Some(SeqRecord {
                    header: prev_header,
                    seq,
                    qual: None,
                })
            } else {
                None
            };
            self.fasta_header = Some(line[1..].to_vec());
            self.fasta_seq.clear();
            Ok(prev)
        } else {
            // Sequence continuation line
            if self.fasta_header.is_none() {
                return Err("FASTA sequence data before header".to_string());
            }
            self.fasta_seq.extend_from_slice(&line);
            Ok(None)
        }
    }
}

fn fmt_commas(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result
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

fn write_record(
    output: &mut impl Write,
    header: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
) -> std::io::Result<()> {
    if let Some(qual) = qual {
        output.write_all(b"@")?;
        output.write_all(header)?;
        output.write_all(b"\n")?;
        output.write_all(seq)?;
        output.write_all(b"\n+\n")?;
        output.write_all(qual)?;
        output.write_all(b"\n")?;
    } else {
        output.write_all(b">")?;
        output.write_all(header)?;
        output.write_all(b"\n")?;
        output.write_all(seq)?;
        output.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::SeqChunkParser;

    #[test]
    fn fastq_parser_handles_byte_sized_chunks() {
        let input = b"@r1\nACGT\n+\n!!!!\n@r2\nTGCA\n+\n####\n";
        let mut parser = SeqChunkParser::new();
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
        assert_eq!(records[0].qual.as_deref(), Some(b"!!!!".as_slice()));
        assert_eq!(records[1].header, b"r2");
        assert_eq!(records[1].seq, b"TGCA");
        assert_eq!(records[1].qual.as_deref(), Some(b"####".as_slice()));
    }

    #[test]
    fn fastq_parser_accepts_final_line_without_trailing_newline() {
        let input = b"@r1\nACGT\n+\n!!!!";
        let mut parser = SeqChunkParser::new();
        parser.push_chunk(input).unwrap();
        let records = parser.finish().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, b"r1");
    }

    #[test]
    fn fastq_parser_rejects_incomplete_record() {
        let mut parser = SeqChunkParser::new();
        parser.push_chunk(b"@r1\nACGT\n+\n").unwrap();
        assert!(parser.finish().is_err());
    }

    #[test]
    fn fastq_parser_rejects_mismatched_quality_length() {
        let mut parser = SeqChunkParser::new();
        let err = match parser.push_chunk(b"@r1\nACGT\n+\n!!!\n") {
            Ok(_) => panic!("expected parser error for mismatched FASTQ quality length"),
            Err(err) => err,
        };
        assert!(err.contains("sequence and quality lengths differ"));
    }

    #[test]
    fn fasta_parser_single_record() {
        let input = b">seq1\nACGT\nTGCA\n";
        let mut parser = SeqChunkParser::new();
        parser.push_chunk(input).unwrap();
        let records = parser.finish().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, b"seq1");
        assert_eq!(records[0].seq, b"ACGTTGCA");
        assert!(records[0].qual.is_none());
    }

    #[test]
    fn fasta_parser_multiple_records() {
        let input = b">s1\nACGT\n>s2\nTGCA\nAAAA\n";
        let mut parser = SeqChunkParser::new();
        let mut records = parser.push_chunk(input).unwrap();
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, b"s1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[1].header, b"s2");
        assert_eq!(records[1].seq, b"TGCAAAAA");
    }

    #[test]
    fn fasta_parser_byte_sized_chunks() {
        let input = b">r1\nACGT\nTGCA\n>r2\nAAAA\n";
        let mut parser = SeqChunkParser::new();
        let mut records = Vec::new();
        for chunk in input.chunks(1) {
            records.append(&mut parser.push_chunk(chunk).unwrap());
        }
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].seq, b"ACGTTGCA");
        assert_eq!(records[1].seq, b"AAAA");
    }

    #[test]
    fn fasta_parser_no_trailing_newline() {
        let input = b">r1\nACGT";
        let mut parser = SeqChunkParser::new();
        parser.push_chunk(input).unwrap();
        let records = parser.finish().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq, b"ACGT");
    }
}
