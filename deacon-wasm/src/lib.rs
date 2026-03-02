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
            self.process_record(
                &record.header,
                &record.seq,
                record.qual.as_deref(),
                &mut raw_output,
            )?;
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
            self.process_record(
                &record.header,
                &record.seq,
                record.qual.as_deref(),
                &mut raw_output,
            )?;
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
    /// FASTQ: number of newlines seen in the current partial record (0..=3)
    newlines_in_partial: u8,
    /// Offset into `pending` where the current incomplete record starts
    record_start: usize,
    /// FASTA: accumulated sequence bytes for the current record
    fasta_seq: Vec<u8>,
    /// FASTA: header of the current record being accumulated
    fasta_header: Option<Vec<u8>>,
}

impl SeqChunkParser {
    fn new() -> Self {
        Self {
            pending: Vec::new(),
            format: None,
            newlines_in_partial: 0,
            record_start: 0,
            fasta_seq: Vec::new(),
            fasta_header: None,
        }
    }

    fn pending_bytes(&self) -> usize {
        self.pending.len()
    }

    fn push_chunk(&mut self, chunk: &[u8]) -> Result<Vec<SeqRecord>, String> {
        if chunk.is_empty() {
            return Ok(Vec::new());
        }

        self.pending.extend_from_slice(chunk);

        // Auto-detect format from the first non-whitespace byte
        if self.format.is_none() {
            let first = self.pending.iter().find(|&&b| b != b'\n' && b != b'\r');
            match first {
                Some(b'@') => self.format = Some(SeqFormat::Fastq),
                Some(b'>') => self.format = Some(SeqFormat::Fasta),
                Some(_) => return Err(
                    "Cannot detect format: first record must start with '@' (FASTQ) or '>' (FASTA)"
                        .to_string(),
                ),
                None => return Ok(Vec::new()), // only whitespace so far
            }
        }

        match self.format.unwrap() {
            SeqFormat::Fastq => self.scan_fastq(),
            SeqFormat::Fasta => self.scan_fasta(),
        }
    }

    fn finish(&mut self) -> Result<Vec<SeqRecord>, String> {
        // If there's trailing data without a final newline, append one to flush it
        if !self.pending.is_empty() {
            let last = *self.pending.last().unwrap();
            if last != b'\n' {
                self.pending.push(b'\n');
            }
        }

        let records = match self.format.unwrap_or(SeqFormat::Fastq) {
            SeqFormat::Fastq => {
                let records = self.scan_fastq()?;
                if self.newlines_in_partial > 0 || self.record_start < self.pending.len() {
                    // Check if remaining data is just whitespace
                    let remaining = &self.pending[self.record_start..];
                    let has_content = remaining.iter().any(|&b| b != b'\n' && b != b'\r');
                    if has_content {
                        return Err(
                            "Incomplete FASTQ record at end of stream (expected 4 lines per record)"
                                .to_string(),
                        );
                    }
                }
                records
            }
            SeqFormat::Fasta => {
                let mut records = self.scan_fasta()?;
                // Flush the last FASTA record
                if let Some(header) = self.fasta_header.take() {
                    let seq = std::mem::take(&mut self.fasta_seq);
                    records.push(SeqRecord {
                        header,
                        seq,
                        qual: None,
                    });
                }
                records
            }
        };

        self.pending.clear();
        self.record_start = 0;
        self.newlines_in_partial = 0;
        Ok(records)
    }

    /// Scan pending buffer for complete FASTQ records (4 newline-delimited lines each).
    fn scan_fastq(&mut self) -> Result<Vec<SeqRecord>, String> {
        let mut records = Vec::new();

        loop {
            // Find the next newline starting from where we left off scanning
            let search_start = self.record_start
                + if self.newlines_in_partial == 0 {
                    0
                } else {
                    // We need to skip past the newlines we've already counted
                    let mut skip = self.record_start;
                    let mut found = 0u8;
                    for &b in &self.pending[self.record_start..] {
                        skip += 1;
                        if b == b'\n' {
                            found += 1;
                            if found == self.newlines_in_partial {
                                break;
                            }
                        }
                    }
                    skip - self.record_start
                };

            // Count remaining newlines needed to complete a record
            let needed = 4 - self.newlines_in_partial;
            let mut pos = search_start;
            let mut found = 0u8;

            for &b in &self.pending[search_start..] {
                pos += 1;
                if b == b'\n' {
                    found += 1;
                    if found == needed {
                        break;
                    }
                }
            }

            if found < needed {
                // Not enough newlines for a complete record
                self.newlines_in_partial += found;
                break;
            }

            // We have a complete 4-line record from record_start..pos
            let record_bytes = &self.pending[self.record_start..pos];
            let record = Self::parse_fastq_record(record_bytes)?;
            records.push(record);

            // Skip any blank lines between records
            let mut next_start = pos;
            while next_start < self.pending.len()
                && (self.pending[next_start] == b'\n' || self.pending[next_start] == b'\r')
            {
                next_start += 1;
            }

            self.record_start = next_start;
            self.newlines_in_partial = 0;
        }

        // Compact: remove consumed data from pending
        if self.record_start > 0 {
            self.pending.drain(..self.record_start);
            self.record_start = 0;
        }

        Ok(records)
    }

    /// Parse a single FASTQ record from a byte slice containing exactly 4 lines (with newlines).
    fn parse_fastq_record(data: &[u8]) -> Result<SeqRecord, String> {
        let mut lines = [0usize; 5]; // start offsets of each line, plus end
        let mut li = 0;
        lines[0] = 0;

        for (i, &b) in data.iter().enumerate() {
            if b == b'\n' {
                li += 1;
                if li <= 4 {
                    lines[li] = i + 1;
                }
            }
        }

        let line = |n: usize| -> &[u8] {
            let start = lines[n];
            let mut end = if n + 1 <= li {
                lines[n + 1] - 1
            } else {
                data.len()
            };
            // Strip \r
            if end > start && data[end - 1] == b'\r' {
                end -= 1;
            }
            &data[start..end]
        };

        let header_line = line(0);
        let seq_line = line(1);
        let sep_line = line(2);
        let qual_line = line(3);

        if header_line.first() != Some(&b'@') {
            return Err("Invalid FASTQ record: header must start with '@'".to_string());
        }
        if sep_line.first() != Some(&b'+') {
            return Err("Invalid FASTQ record: third line must start with '+'".to_string());
        }
        if qual_line.len() != seq_line.len() {
            return Err("Invalid FASTQ record: sequence and quality lengths differ".to_string());
        }

        Ok(SeqRecord {
            header: header_line[1..].to_vec(),
            seq: seq_line.to_vec(),
            qual: Some(qual_line.to_vec()),
        })
    }

    /// Scan pending buffer for complete FASTA records (delimited by \n>).
    fn scan_fasta(&mut self) -> Result<Vec<SeqRecord>, String> {
        let mut records = Vec::new();

        loop {
            // Look for \n> which marks the start of a new record
            let search_from = if self.record_start == 0
                && self.fasta_header.is_none()
                && !self.pending.is_empty()
            {
                // First record: find the first > (should be at the start, possibly after whitespace)
                match self.pending.iter().position(|&b| b == b'>') {
                    Some(pos) => pos,
                    None => break,
                }
            } else if self.fasta_header.is_some() {
                // We have an active record; look for \n> to end it
                let start = self.record_start;
                let found = {
                    let mut pos = start;
                    loop {
                        match memchr::memchr(b'>', &self.pending[pos..]) {
                            Some(rel) => {
                                let abs = pos + rel;
                                // Valid record boundary if > is preceded by \n (or is at
                                // record_start after all preceding data was accumulated)
                                if abs > 0 && self.pending[abs - 1] == b'\n' {
                                    break Some(abs);
                                } else if abs == start && (self.fasta_seq.len() > 0 || abs == 0) {
                                    // > is right at record_start — previous data was already
                                    // accumulated (the \n was consumed/drained in a prior chunk)
                                    break Some(abs);
                                } else {
                                    // This > is inside sequence data; skip it
                                    pos = abs + 1;
                                }
                            }
                            None => break None,
                        }
                    }
                };
                match found {
                    Some(pos) => pos,
                    None => {
                        // No \n> found. Accumulate all complete lines into fasta_seq.
                        self.accumulate_fasta_seq();
                        break;
                    }
                }
            } else {
                break;
            };

            if self.fasta_header.is_some() {
                // Flush sequence data up to this boundary
                let end = if search_from > self.record_start {
                    search_from - 1 // exclude the \n before >
                } else {
                    self.record_start // nothing between record_start and >
                };
                for i in self.record_start..end {
                    let b = self.pending[i];
                    if b != b'\n' && b != b'\r' {
                        self.fasta_seq.push(b);
                    }
                }
                let header = self.fasta_header.take().unwrap();
                let seq = std::mem::take(&mut self.fasta_seq);
                records.push(SeqRecord {
                    header,
                    seq,
                    qual: None,
                });
            }

            // Parse the new header line starting at search_from
            let header_start = search_from + 1; // skip '>'
            match memchr::memchr(b'\n', &self.pending[header_start..]) {
                Some(rel) => {
                    let header_end = header_start + rel;
                    let mut header = self.pending[header_start..header_end].to_vec();
                    if header.last() == Some(&b'\r') {
                        header.pop();
                    }
                    self.fasta_header = Some(header);
                    self.record_start = header_end + 1;
                    self.fasta_seq.clear();
                }
                None => {
                    // Header line is incomplete; wait for more data
                    // Put record_start at the '>' so we re-process it next time
                    self.record_start = search_from;
                    break;
                }
            }
        }

        // Compact: remove consumed data from pending
        if self.record_start > 0 {
            self.pending.drain(..self.record_start);
            self.record_start = 0;
        }

        Ok(records)
    }

    /// Accumulate all complete lines from pending[record_start..] into fasta_seq,
    /// advancing record_start past consumed data.
    fn accumulate_fasta_seq(&mut self) {
        // Find the last newline — only accumulate complete lines
        let start = self.record_start;
        let last_nl = self.pending[start..].iter().rposition(|&b| b == b'\n');
        if let Some(last_nl) = last_nl {
            let end = start + last_nl;
            for i in start..end {
                let b = self.pending[i];
                if b != b'\n' && b != b'\r' {
                    self.fasta_seq.push(b);
                }
            }
            self.record_start = end + 1;
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

    #[test]
    fn fastq_quality_line_starting_with_at() {
        let input = b"@r1\nACGT\n+\n@!!!\n";
        let mut parser = SeqChunkParser::new();
        let records = parser.push_chunk(input).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, b"r1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[0].qual.as_deref(), Some(b"@!!!".as_slice()));
    }

    #[test]
    fn fastq_quality_line_starting_with_plus() {
        let input = b"@r1\nACGT\n+\n+!!!\n";
        let mut parser = SeqChunkParser::new();
        let records = parser.push_chunk(input).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].qual.as_deref(), Some(b"+!!!".as_slice()));
    }

    #[test]
    fn fastq_crlf_line_endings() {
        let input = b"@r1\r\nACGT\r\n+\r\n!!!!\r\n";
        let mut parser = SeqChunkParser::new();
        let mut records = parser.push_chunk(input).unwrap();
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, b"r1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[0].qual.as_deref(), Some(b"!!!!".as_slice()));
    }

    #[test]
    fn fasta_crlf_line_endings() {
        let input = b">s1\r\nACGT\r\nTGCA\r\n";
        let mut parser = SeqChunkParser::new();
        parser.push_chunk(input).unwrap();
        let records = parser.finish().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, b"s1");
        assert_eq!(records[0].seq, b"ACGTTGCA");
    }

    #[test]
    fn fastq_record_split_at_chunk_boundary() {
        let mut parser = SeqChunkParser::new();
        let mut records = Vec::new();
        // Header in first chunk, rest in second
        records.append(&mut parser.push_chunk(b"@r1\n").unwrap());
        assert_eq!(records.len(), 0);
        records.append(&mut parser.push_chunk(b"ACGT\n+\n!!!!\n").unwrap());
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].header, b"r1");
        assert_eq!(records[0].seq, b"ACGT");
    }

    #[test]
    fn fasta_record_split_at_chunk_boundary() {
        let mut parser = SeqChunkParser::new();
        let mut records = Vec::new();
        records.append(&mut parser.push_chunk(b">r1\nAC").unwrap());
        assert_eq!(records.len(), 0);
        records.append(&mut parser.push_chunk(b"GT\nTGCA\n>r2\nAAAA\n").unwrap());
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, b"r1");
        assert_eq!(records[0].seq, b"ACGTTGCA");
        assert_eq!(records[1].seq, b"AAAA");
    }

    #[test]
    fn fastq_long_sequence_spanning_chunks() {
        let seq = "A".repeat(10000);
        let qual = "!".repeat(10000);
        let input = format!("@long\n{}\n+\n{}\n", seq, qual);
        let bytes = input.as_bytes();
        let mut parser = SeqChunkParser::new();
        let mut records = Vec::new();
        // Feed in 137-byte chunks (odd size to stress boundaries)
        for chunk in bytes.chunks(137) {
            records.append(&mut parser.push_chunk(chunk).unwrap());
        }
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq.len(), 10000);
        assert_eq!(records[0].qual.as_ref().unwrap().len(), 10000);
    }

    #[test]
    fn fasta_long_sequence_spanning_chunks() {
        let seq_line = "ACGTACGT".repeat(125); // 1000 bp per line
        let input = format!(">long\n{}\n{}\n{}\n", seq_line, seq_line, seq_line);
        let bytes = input.as_bytes();
        let mut parser = SeqChunkParser::new();
        let mut records = Vec::new();
        for chunk in bytes.chunks(137) {
            records.append(&mut parser.push_chunk(chunk).unwrap());
        }
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq.len(), 3000);
    }

    #[test]
    fn fasta_blank_lines_between_records() {
        let input = b">s1\nACGT\n\n>s2\nTGCA\n\n";
        let mut parser = SeqChunkParser::new();
        let mut records = parser.push_chunk(input).unwrap();
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[1].seq, b"TGCA");
    }

    #[test]
    fn fasta_allows_empty_record_mid_stream() {
        let input = b">s1\nACGT\n>empty\n>s2\nTGCA\n";
        let mut parser = SeqChunkParser::new();
        let mut records = parser.push_chunk(input).unwrap();
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].header, b"s1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[1].header, b"empty");
        assert_eq!(records[1].seq, b"");
        assert_eq!(records[2].header, b"s2");
        assert_eq!(records[2].seq, b"TGCA");
    }

    #[test]
    fn fasta_allows_empty_record_at_eof() {
        let input = b">s1\nACGT\n>empty\n";
        let mut parser = SeqChunkParser::new();
        let mut records = parser.push_chunk(input).unwrap();
        records.append(&mut parser.finish().unwrap());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, b"s1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[1].header, b"empty");
        assert_eq!(records[1].seq, b"");
    }

    #[test]
    fn fastq_multiple_records_single_chunk() {
        let input = b"@r1\nAA\n+\n!!\n@r2\nCC\n+\n##\n@r3\nGG\n+\n$$\n";
        let mut parser = SeqChunkParser::new();
        let records = parser.push_chunk(input).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].header, b"r1");
        assert_eq!(records[1].header, b"r2");
        assert_eq!(records[2].header, b"r3");
    }
}
