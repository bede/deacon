use packed_seq::SeqVec;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[cfg(feature = "server")]
use reqwest::blocking::Client;

// JSON summary structure
#[derive(Serialize, Deserialize)]
pub struct FilterSummary {
    pub version: String,
    pub index: String,
    pub input: String,
    pub input2: Option<String>,
    pub output: String,
    pub output2: Option<String>,
    pub k: u8,
    pub w: u8,
    pub abs_threshold: usize,
    pub rel_threshold: f64,
    pub prefix_length: usize,
    pub deplete: bool,
    pub rename: bool,
    pub seqs_in: u64,
    pub seqs_out: u64,
    pub seqs_out_proportion: f64,
    pub seqs_removed: u64,
    pub seqs_removed_proportion: f64,
    pub bp_in: u64,
    pub bp_out: u64,
    pub bp_out_proportion: f64,
    pub bp_removed: u64,
    pub bp_removed_proportion: f64,
    pub time: f64,
    pub seqs_per_second: u64,
    pub bp_per_second: u64,
}

/// Get a summary string for the index used, either from the local path or by querying the server.
/// # Args:
/// * `minimizers_path`: Optional path to the local minimizer index.
/// * `server_address`: Optional server address to query for index version.
/// # Returns:
/// * A string summarizing the index used. If local, it's the path; if from server, it's "address:filename@hash".
pub fn get_summary_index(
    minimizers_path: &Option<&PathBuf>,
    server_address: &Option<String>,
) -> String {
    let index = match minimizers_path {
        Some(path) => path.to_string_lossy().to_string(),
        None => match &server_address {
            None => "No index or server specified".to_string(),
            Some(_addr) => {
                #[cfg(feature = "server")]
                {
                    let client = Client::new();
                    let response = client
                        .get(_addr.to_owned() + "/index_version")
                        .send()
                        .unwrap_or_else(|e| {
                            panic!("Failed to contact server at {}: {e}", _addr);
                        });
                    if response.status().is_success() {
                        _addr.to_owned()
                            + ":"
                            + &response.text().unwrap_or_else(|e| {
                                panic!("Failed to parse server response: {e}");
                            })
                    } else {
                        panic!("Server returned error: {}", response.status())
                    }
                }
                #[cfg(not(feature = "server"))]
                {
                    panic!("Server feature not enabled, cannot use server address");
                }
            }
        },
    };
    index
}

/// Calculate required hits based on absolute and relative thresholds
pub fn calculate_required_hits(
    abs_threshold: usize,
    rel_threshold: f64,
    total_minimizers: usize,
) -> usize {
    let abs_required = abs_threshold;
    let rel_required = if total_minimizers == 0 {
        0
    } else {
        ((rel_threshold * total_minimizers as f64).round() as usize).max(1)
    };
    abs_required.max(rel_required)
}

/// Check if sequence meets filtering criteria
pub fn meets_filtering_criteria(
    hit_count: usize,
    total_minimizers: usize,
    abs_threshold: usize,
    rel_threshold: f64,
    deplete: bool,
) -> bool {
    let required = calculate_required_hits(abs_threshold, rel_threshold, total_minimizers);
    if deplete {
        hit_count < required
    } else {
        hit_count >= required
    }
}

/// Check how many minimizers from the sequence are in the set of minimizer hashes
/// and optionally collect the matching k-mers for debugging. Unpaired equivalent of `pair_matches`.
/// Used in both local and remote filtering.
///
/// # Args:
/// * `minimizer_hashes`: Set of minimizer hashes to check against.
/// * `minimizer_values`: Minimizer hashes from the sequence.
/// * `positions`: Positions of the minimizers in the sequence.
/// * `effective_seq`: The effective sequence used for minimizer calculation.
/// * `kmer_length`: Length of the k-mers.
/// * `debug`: If true, collect matching k-mers.
/// # Returns:
/// * A tuple containing:
/// - The count of distinct minimizer hits.
/// - A vector of matching k-mers as strings (if debug is true).
pub fn sequence_matches(
    minimizer_hashes: &FxHashSet<u64>,
    minimizer_values: &[u64],
    positions: &[u32],
    effective_seq: &[u8],
    kmer_length: u8,
    debug: bool,
) -> (usize, Vec<String>) {
    // Count distinct minimizer hits and collect matching k-mers
    let mut seen_hits = FxHashSet::default();
    let mut hit_count = 0;
    let mut hit_kmers = Vec::new();

    // Should keep sequence if it meets filtering criteria
    for (i, &hash) in minimizer_values.iter().enumerate() {
        if minimizer_hashes.contains(&hash) && seen_hits.insert(hash) {
            hit_count += 1;
            // Extract the k-mer sequence at this position
            if debug && i < positions.len() {
                let pos = positions[i] as usize;
                let kmer = &effective_seq[pos..pos + kmer_length as usize];
                hit_kmers.push(String::from_utf8_lossy(kmer).to_string());
            }
        }
    }
    (hit_count, hit_kmers)
}

/// Check how many minimizers from both sequences in a pair are in the set of minimizer hashes
/// and optionally collect the matching k-mers for debugging. Paired equivalent of `sequence_matches`.
/// Used in both local and remote filtering.
///
/// # Args:
/// * `all_hashes`: Combined minimizer hashes from both sequences.
/// * `all_positions`: Combined positions of the minimizers in both sequences.
/// * `all_sequences`: Combined effective sequences from both sequences.
/// * `minimizer_hashes`: Set of minimizer hashes to check against.
/// * `kmer_length`: Length of the k-mers.
/// * `debug`: If true, collect matching k-mers.
/// # Returns:
/// * A tuple containing:
/// - The count of distinct minimizer hits across both sequences.
/// - A vector of matching k-mers as strings (if debug is true).
pub fn pair_matches(
    all_hashes: &Vec<u64>,
    all_positions: &Vec<u32>,
    all_sequences: &Vec<&[u8]>,
    minimizer_hashes: &FxHashSet<u64>,
    kmer_length: u8,
    debug: bool,
) -> (usize, Vec<String>) {
    let mut seen_hits_pair = FxHashSet::default();
    let mut pair_hit_count = 0;
    let mut hit_kmers = Vec::new();
    // Count hits and collect k-mers
    for (i, &hash) in all_hashes.iter().enumerate() {
        if minimizer_hashes.contains(&hash) && seen_hits_pair.insert(hash) {
            pair_hit_count += 1;
            if debug && i < all_positions.len() && i < all_sequences.len() {
                let pos = all_positions[i] as usize;
                let seq = all_sequences[i];
                if pos + kmer_length as usize <= seq.len() {
                    let kmer = &seq[pos..pos + kmer_length as usize];
                    hit_kmers.push(String::from_utf8_lossy(kmer).to_string());
                }
            }
        }
    }
    (pair_hit_count, hit_kmers)
}

/// Given a sequence, compute the minimizer hashes and positons
/// # Args:
/// * `seq`: The input sequence as a byte slice.
/// * `prefix_length`: If >0, only consider the first `prefix_length` bases of the sequence.
/// * `kmer_length`: The length of k-mers to consider for minimizers.
/// * `window_size`: The size of the sliding window to find minimizers.
/// # Returns:
/// * A tuple containing:
/// - A vector of minimizer hash values (u64).
/// - A vector of positions (u32) where each minimizer occurs in the sequence.
/// - A slice of the effective sequence used for minimizer calculation (after applying prefix length and trimming).
pub fn get_minimizer_hashes_and_positions(
    seq: &[u8],
    prefix_length: usize,
    kmer_length: u8,
    window_size: u8,
) -> (Vec<u64>, Vec<u32>, &[u8]) {
    if seq.len() < kmer_length as usize {
        return (Vec::new(), Vec::new(), &[]); // If too short, return nothing
    }

    // Apply prefix length limit if specified
    let effective_seq = if prefix_length > 0 && seq.len() > prefix_length {
        &seq[..prefix_length]
    } else {
        seq
    };

    // Trim the last newline character from `effective_seq` if it has one.
    let effective_seq = effective_seq.strip_suffix(b"\n").unwrap_or(effective_seq);

    let mut invalid_mask = Vec::new();
    let mut positions = Vec::new();
    let mut minimizer_values = Vec::new();

    // Pack the sequence into 2-bit representation.
    // Any non-ACGT characters are silently converted to 2-bit ACGT as well.
    // packed_seq.push_ascii(effective_seq);
    let packed_seq = packed_seq::PackedSeqVec::from_ascii(effective_seq);
    // let packed_seq = packed_seq::PackedSeqVec::from_ascii(effective_seq);

    // TODO: Extract this to some nicer helper function in packed_seq?
    // TODO: Use SIMD?
    // TODO: Should probably add some test for this.
    // +2: one to round up, and one buffer.
    invalid_mask.resize(packed_seq.len() / 64 + 2, 0);
    // let mut invalid_mask = vec![0u64; packed_seq.len() / 64 + 2];
    for i in (0..effective_seq.len()).step_by(64) {
        let mut mask = 0;
        for (j, b) in effective_seq[i..(i + 64).min(effective_seq.len())]
            .iter()
            .enumerate()
        {
            mask |=
                ((!matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) as u64) << j;
        }

        invalid_mask[i / 64] = mask;
    }

    // let mut positions = Vec::new();
    simd_minimizers::canonical_minimizer_positions(
        // packed_seq::AsciiSeq(&canonical_seq),
        packed_seq.as_slice(),
        kmer_length as usize,
        window_size as usize,
        &mut positions,
    );

    assert!(
        kmer_length <= 56,
        "Indexing the bitmask of invalid characters requires k<=56, but it is {}",
        kmer_length
    );

    // Filter positions to only include k-mers with ACGT bases
    positions.retain(|&pos| {
        // Extract bits pos .. pos+k from the bitmask.

        // mask of k ones in low positions.
        let mask = u64::MAX >> (64 - kmer_length);
        let byte = pos as usize / 8;
        let offset = pos as usize % 8;
        // The unaligned u64 read is OK, because we ensure that the underlying `Vec` always
        // has at least 8 bytes of padding at the end.
        let x = (unsafe { invalid_mask.as_ptr().byte_add(byte).read_unaligned() } >> offset) & mask;
        x == 0
    });

    // Hash valid positions
    if kmer_length > 32 {
        minimizer_values.extend(
            simd_minimizers::iter_canonical_minimizer_values_u128(
                packed_seq.as_slice(),
                kmer_length as usize,
                &positions,
            )
            .map(|kmer| xxhash_rust::xxh3::xxh3_64(&kmer.to_le_bytes())),
        );
    } else {
        minimizer_values.extend(
            simd_minimizers::iter_canonical_minimizer_values(
                packed_seq.as_slice(),
                kmer_length as usize,
                &positions,
            )
            .map(|kmer| xxhash_rust::xxh3::xxh3_64(&kmer.to_le_bytes())),
        );
    }

    (minimizer_values, positions, effective_seq)
}

pub fn get_paired_minimizer_hashes_and_positions<'a>(
    seq1: &'a [u8],
    seq2: &'a [u8],
    prefix_length: usize,
    kmer_length: u8,
    window_size: u8,
) -> (Vec<u64>, Vec<u32>, Vec<&'a [u8]>) {
    let mut all_hashes = Vec::new();
    let mut all_positions = Vec::new();
    let mut all_sequences = Vec::new();

    // Process read 1
    if seq1.len() >= kmer_length as usize {
        let (hashes, positions, effective_seq1) =
            get_minimizer_hashes_and_positions(seq1, prefix_length, kmer_length, window_size);
        all_hashes.extend(hashes);
        all_positions.extend(positions);
        all_sequences.extend(vec![effective_seq1; all_hashes.len() - all_positions.len()]);
    }

    // Process read 2
    if seq2.len() >= kmer_length as usize {
        let (hashes, positions, effective_seq2) =
            get_minimizer_hashes_and_positions(seq2, prefix_length, kmer_length, window_size);
        all_hashes.extend(hashes);
        all_positions.extend(positions);
        all_sequences.extend(vec![effective_seq2; all_hashes.len() - all_positions.len()]);
    }

    (all_hashes, all_positions, all_sequences)
}
