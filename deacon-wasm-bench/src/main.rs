use std::fs;
use std::io::{Cursor, Write};
use std::time::Instant;

use flate2::write::MultiGzDecoder;

use deacon::minimizers::{Buffers, KmerHasher};
use deacon::{MinimizerSet, MinimizerVec, RapidHashSet};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <index.idx> <reads.fastq> [chunk_size_kb] [iterations]",
            args[0]
        );
        std::process::exit(1);
    }

    let index_path = &args[1];
    let reads_path = &args[2];
    let chunk_kb: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(256);
    let chunk_size = chunk_kb * 1024;
    let iterations: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(3);

    eprintln!("Loading index from {}...", index_path);
    let t0 = Instant::now();
    let index_data = fs::read(index_path).expect("failed to read index file");
    let mut cursor = Cursor::new(&index_data);
    let (minimizers, header) =
        deacon::load_minimizers(&mut cursor).expect("failed to parse index");
    eprintln!(
        "Index loaded in {:.2}s (k={}, w={}, {} minimizers)",
        t0.elapsed().as_secs_f64(),
        header.kmer_length(),
        header.window_size(),
        minimizers.len()
    );

    eprintln!("Loading reads from {}...", reads_path);
    let reads_data = fs::read(reads_path).expect("failed to read FASTQ file");
    let reads_mb = reads_data.len() as f64 / (1024.0 * 1024.0);
    let gz_input = reads_path.ends_with(".gz");
    eprintln!("Reads loaded: {:.1} MB (gz={})", reads_mb, gz_input);

    let k = header.kmer_length();
    let w = header.window_size();

    eprintln!(
        "\nBenchmark: chunk_size={}KB, iterations={}",
        chunk_kb, iterations
    );
    eprintln!("{:-<60}", "");

    let mut times = Vec::with_capacity(iterations);

    for i in 0..iterations {
        let t = Instant::now();
        let (reads_in, reads_out, _bases_in) =
            run_filter(&reads_data, chunk_size, &minimizers, k, w, true, 1, 0.0, gz_input);
        let secs = t.elapsed().as_secs_f64();
        let throughput = reads_mb / secs;
        times.push(secs);

        eprintln!(
            "  iter {}: {:.3}s ({:.1} MB/s) — {}/{} reads retained",
            i + 1,
            secs,
            throughput,
            reads_out,
            reads_in,
        );
    }

    let mean = times.iter().sum::<f64>() / times.len() as f64;
    let mean_throughput = reads_mb / mean;
    eprintln!("{:-<60}", "");
    eprintln!("Mean: {:.3}s ({:.1} MB/s)", mean, mean_throughput);

    println!("RESULT\t{:.3}\t{:.1}\t{}", mean, mean_throughput, chunk_kb);
}

fn run_filter(
    data: &[u8],
    chunk_size: usize,
    index: &MinimizerSet,
    k: u8,
    w: u8,
    deplete: bool,
    abs_threshold: usize,
    rel_threshold: f64,
    gz_input: bool,
) -> (u64, u64, u64) {
    let hasher = KmerHasher::new(k as usize);
    let mut buffers = if k <= 32 {
        Buffers::new_u64()
    } else {
        Buffers::new_u128()
    };
    let mut seen = if k <= 32 {
        MinimizerSet::U64(RapidHashSet::default())
    } else {
        MinimizerSet::U128(RapidHashSet::default())
    };

    let mut gz_decoder = if gz_input {
        Some(MultiGzDecoder::new(Vec::new()))
    } else {
        None
    };

    let mut pending: Vec<u8> = Vec::new();
    let mut record_start: usize = 0;
    let mut newlines_in_partial: u8 = 0;
    let mut reads_in: u64 = 0;
    let mut reads_out: u64 = 0;
    let mut bases_in: u64 = 0;
    let mut output_buf: Vec<u8> = Vec::new();

    for chunk in data.chunks(chunk_size) {
        let decompressed;
        let input = if let Some(ref mut decoder) = gz_decoder {
            decoder.write_all(chunk).expect("decompression error");
            decompressed = std::mem::take(decoder.get_mut());
            &decompressed[..]
        } else {
            chunk
        };
        pending.extend_from_slice(input);

        loop {
            let search_start =
                record_start + skip_counted_newlines(&pending, record_start, newlines_in_partial);
            let needed = 4 - newlines_in_partial;
            let mut pos = search_start;
            let mut found = 0u8;

            for &b in &pending[search_start..] {
                pos += 1;
                if b == b'\n' {
                    found += 1;
                    if found == needed {
                        break;
                    }
                }
            }

            if found < needed {
                newlines_in_partial += found;
                break;
            }

            let record_bytes = &pending[record_start..pos];

            if let Some((header, seq, qual)) = parse_fastq_fields(record_bytes) {
                reads_in += 1;
                bases_in += seq.len() as u64;

                let keep = if seq.len() < k as usize {
                    deplete
                } else {
                    deacon::minimizers::fill_minimizers(
                        seq, &hasher, k, w, 0.0, &mut buffers,
                    );
                    let num_minimizers = buffers.minimizers.len();
                    let hit_count = count_hits(&buffers.minimizers, index, &mut seen);
                    let rel_required = if num_minimizers == 0 {
                        0
                    } else {
                        ((rel_threshold * num_minimizers as f64).round() as usize).max(1)
                    };
                    let required = abs_threshold.max(rel_required);
                    if deplete {
                        hit_count < required
                    } else {
                        hit_count >= required
                    }
                };

                if keep {
                    reads_out += 1;
                    output_buf.clear();
                    output_buf.push(b'@');
                    output_buf.extend_from_slice(header);
                    output_buf.push(b'\n');
                    output_buf.extend_from_slice(seq);
                    output_buf.extend_from_slice(b"\n+\n");
                    output_buf.extend_from_slice(qual);
                    output_buf.push(b'\n');
                }
            }

            let mut next_start = pos;
            while next_start < pending.len()
                && (pending[next_start] == b'\n' || pending[next_start] == b'\r')
            {
                next_start += 1;
            }
            record_start = next_start;
            newlines_in_partial = 0;
        }

        if record_start > 0 {
            pending.drain(..record_start);
            record_start = 0;
        }
    }

    (reads_in, reads_out, bases_in)
}

fn skip_counted_newlines(pending: &[u8], record_start: usize, count: u8) -> usize {
    if count == 0 {
        return 0;
    }
    let mut skip = 0;
    let mut found = 0u8;
    for &b in &pending[record_start..] {
        skip += 1;
        if b == b'\n' {
            found += 1;
            if found == count {
                break;
            }
        }
    }
    skip
}

fn parse_fastq_fields(data: &[u8]) -> Option<(&[u8], &[u8], &[u8])> {
    let mut lines = [0usize; 5];
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
        if end > start && data[end - 1] == b'\r' {
            end -= 1;
        }
        &data[start..end]
    };

    let header_line = line(0);
    let seq_line = line(1);
    let qual_line = line(3);

    if header_line.first() != Some(&b'@') {
        return None;
    }

    Some((&header_line[1..], seq_line, qual_line))
}

fn count_hits(minimizers: &MinimizerVec, set: &MinimizerSet, seen: &mut MinimizerSet) -> usize {
    match (minimizers, set, seen) {
        (MinimizerVec::U64(vec), MinimizerSet::U64(s), MinimizerSet::U64(seen)) => {
            seen.clear();
            for &m in vec {
                if s.contains(&m) {
                    seen.insert(m);
                }
            }
            seen.len()
        }
        (MinimizerVec::U128(vec), MinimizerSet::U128(s), MinimizerSet::U128(seen)) => {
            seen.clear();
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
