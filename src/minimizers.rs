use packed_seq::{SeqVec, unpack_base};

use crate::filter::Buffers;

pub const DEFAULT_KMER_LENGTH: u8 = 31;
pub const DEFAULT_WINDOW_SIZE: u8 = 15;

/// Supported complexity measures for k-mer filtering
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComplexityMeasure {
    /// Shannon entropy scaled to [0, 1] range (default)
    ScaledShannon,
    /// Linguistic complexity (vocabulary richness, W=4)
    Linguistic,
    /// SDust complexity (triplet redundancy scoring)
    SDust,
}

impl ComplexityMeasure {
    /// Parse complexity measure from string
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "scaled-shannon" => Ok(ComplexityMeasure::ScaledShannon),
            "linguistic" => Ok(ComplexityMeasure::Linguistic),
            "sdust" => Ok(ComplexityMeasure::SDust),
            _ => Err(format!(
                "Unknown complexity measure: '{}'. Supported: scaled-shannon, linguistic, sdust",
                s
            )),
        }
    }

    /// Get the string representation of the measure
    pub fn as_str(&self) -> &'static str {
        match self {
            ComplexityMeasure::ScaledShannon => "scaled-shannon",
            ComplexityMeasure::Linguistic => "linguistic",
            ComplexityMeasure::SDust => "sdust",
        }
    }

    /// Calculate complexity for a k-mer using the specified measure
    #[inline]
    pub fn calculate(&self, kmer: &[u8], kmer_length: u8) -> f32 {
        match self {
            ComplexityMeasure::ScaledShannon => calculate_scaled_shannon_entropy(kmer, kmer_length),
            ComplexityMeasure::Linguistic => calculate_linguistic_complexity(kmer, kmer_length),
            ComplexityMeasure::SDust => calculate_sdust_score(kmer, kmer_length),
        }
    }
}

impl Default for ComplexityMeasure {
    fn default() -> Self {
        ComplexityMeasure::ScaledShannon
    }
}

/// Decode u64 minimizer (2-bit canonical k-mer) to ASCII
pub fn decode_u64(minimizer: u64, k: u8) -> Vec<u8> {
    (0..k)
        .map(|i| {
            let base_bits = ((minimizer >> (2 * i)) & 0b11) as u8;
            unpack_base(base_bits)
        })
        .collect()
}

/// Decode u128 minimizer (2-bit canonical k-mer) to ASCII
pub fn decode_u128(minimizer: u128, k: u8) -> Vec<u8> {
    (0..k)
        .map(|i| {
            let base_bits = ((minimizer >> (2 * i)) & 0b11) as u8;
            unpack_base(base_bits)
        })
        .collect()
}

/// Canonical NtHash, with 1-bit rotations for backwards compatibility.
pub type KmerHasher = simd_minimizers::seq_hash::NtHasher<true, 1>;

/// Returns vector of all minimizers for a sequence
pub fn compute_minimizers(
    seq: &[u8],
    hasher: &KmerHasher,
    kmer_length: u8,
    window_size: u8,
    complexity_threshold: f32,
) -> crate::MinimizerVec {
    let mut buffers = if kmer_length <= 32 {
        Buffers::new_u64()
    } else {
        Buffers::new_u128()
    };
    fill_minimizers(
        seq,
        hasher,
        kmer_length,
        window_size,
        complexity_threshold,
        ComplexityMeasure::default(),
        &mut buffers,
    );
    buffers.minimizers
}

/// Calculate linguistic complexity using vocabulary richness (W=4)
/// Returns complexity between 0.0 and 1.0
#[inline]
fn calculate_linguistic_complexity(kmer: &[u8], kmer_length: u8) -> f32 {
    use rustc_hash::FxHashSet;

    // K-mers less than 10 bases long always pass filter
    if kmer_length < 10 {
        return 1.0;
    }

    let mut lc = 1.0f32;
    let k = kmer_length as usize;

    // For each word length w from 1 to 4
    for w in 1..=4 {
        let mut seen = FxHashSet::default();

        // Count distinct w-grams
        for i in 0..=(k - w) {
            seen.insert(&kmer[i..i + w]);
        }

        // Calculate Uw = (distinct w-grams) / (max possible w-grams)
        // Max possible is 4^w (all possible DNA w-grams)
        let max_possible = 4_usize.pow(w as u32);
        let u_w = seen.len() as f32 / max_possible as f32;

        lc *= u_w;
    }

    lc
}

/// Calculate sdust complexity using triplet redundancy
/// Returns normalized complexity score between 0.0 and 1.0 where:
/// - 1.0 = high complexity (all unique triplets)
/// - 0.0 = low complexity (homopolymer with all identical triplets)
/// Based on the symmetric DUST algorithm from Morgulis et al. 2006
#[inline]
fn calculate_sdust_score(kmer: &[u8], kmer_length: u8) -> f32 {
    const TRIPLET_SIZE: usize = 3;
    const NUM_TRIPLETS: usize = 64; // 4^3 possible DNA triplets

    // K-mers less than 3 bases long cannot form triplets, return high complexity
    if kmer_length < TRIPLET_SIZE as u8 {
        return 1.0;
    }

    let k = kmer_length as usize;

    // Count triplet frequencies using fixed array (A=0, C=1, G=2, T=3)
    let mut triplet_counts = [0u32; NUM_TRIPLETS];

    // Convert base to 2-bit representation
    #[inline]
    fn base_to_bits(base: u8) -> Option<u8> {
        match base {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    // Slide through k-mer and count triplets
    for i in 0..=(k - TRIPLET_SIZE) {
        // Convert triplet to index (0-63)
        if let (Some(b0), Some(b1), Some(b2)) = (
            base_to_bits(kmer[i]),
            base_to_bits(kmer[i + 1]),
            base_to_bits(kmer[i + 2]),
        ) {
            let triplet_idx = ((b0 as usize) << 4) | ((b1 as usize) << 2) | (b2 as usize);
            triplet_counts[triplet_idx] += 1;
        }
    }

    // Calculate raw sdust score: sum of (count - 1) for each triplet
    // This measures redundancy - higher raw score means more repetitive
    let mut raw_score = 0u32;
    for &count in &triplet_counts {
        if count > 1 {
            raw_score += count - 1;
        }
    }

    // Normalize to [0, 1] complexity scale
    // Max possible score for k-mer length k is (k - 3) when all triplets are identical
    // Transform: complexity = 1.0 - (raw_score / max_score)
    let max_score = (k - TRIPLET_SIZE) as f32;
    let complexity = 1.0 - (raw_score as f32 / max_score);

    complexity
}

/// Calculate scaled Shannon entropy using character frequency analysis
/// Returns scaled entropy between 0.0 and 1.0
#[inline]
fn calculate_scaled_shannon_entropy(kmer: &[u8], kmer_length: u8) -> f32 {
    // K-mers less than 10 bases long always pass filter
    if kmer_length < 10 {
        return 1.0;
    }

    // Count character frequencies using fixed array (faster than HashMap)
    let mut counts = [0u8; 4]; // A, C, G, T
    let mut total = 0u8;

    // Iterate only up to kmer_length to avoid bounds checks
    for &base in kmer.iter().take(kmer_length as usize) {
        match base {
            b'A' | b'a' => {
                counts[0] += 1;
                total += 1;
            }
            b'C' | b'c' => {
                counts[1] += 1;
                total += 1;
            }
            b'G' | b'g' => {
                counts[2] += 1;
                total += 1;
            }
            b'T' | b't' => {
                counts[3] += 1;
                total += 1;
            }
            _ => {} // Skip invalid characters
        }
    }

    if total == 0 {
        return 1.0; // All non-ACGT, don't filter
    }

    let total_f32 = total as f32;
    let mut entropy = 0.0;
    for &count in &counts {
        if count > 0 {
            let p = count as f32 / total_f32;
            entropy -= p * p.log2();
        }
    }

    // Scale entropy to [0, 1] range (max entropy for 4 bases is 2.0)
    entropy / 2.0
}

/// Fill a vector with minimizers, skipping k-mers with non-ACGT bases
/// and optionally filtering by complexity
pub(crate) fn fill_minimizers(
    seq: &[u8],
    hasher: &KmerHasher,
    kmer_length: u8,
    window_size: u8,
    complexity_threshold: f32,
    complexity_measure: ComplexityMeasure,
    buffers: &mut Buffers,
) {
    let Buffers {
        packed_nseq,
        positions,
        minimizers,
        cache,
    } = buffers;

    packed_nseq.seq.clear();
    packed_nseq.ambiguous.clear();
    minimizers.clear();
    positions.clear();

    // Skip if sequence is too short
    if seq.len() < kmer_length as usize {
        return;
    }

    // Pack the sequence into 2-bit representation.
    packed_nseq.seq.push_ascii(seq);
    packed_nseq.ambiguous.push_ascii(seq);

    // Get minimizer positions using simd-minimizers
    let out = simd_minimizers::canonical_minimizers(kmer_length as usize, window_size as usize)
        .hasher(hasher)
        .run_skip_ambiguous_windows_with_buf(packed_nseq.as_slice(), positions, cache);

    match minimizers {
        crate::MinimizerVec::U64(vec) => {
            vec.extend(
                out.pos_and_values_u64()
                    .filter(|&(pos, _val)| {
                        let pos_usize = pos as usize;
                        let kmer = &seq[pos_usize..pos_usize + kmer_length as usize];

                        // Check complexity constraint if threshold specified
                        if complexity_threshold == 0.0 {
                            true
                        } else {
                            let complexity = complexity_measure.calculate(kmer, kmer_length);
                            complexity >= complexity_threshold
                        }
                    })
                    .map(|(_pos, val)| val),
            );
        }
        crate::MinimizerVec::U128(vec) => {
            vec.extend(
                out.pos_and_values_u128()
                    .filter(|&(pos, _val)| {
                        let pos_usize = pos as usize;
                        let kmer = &seq[pos_usize..pos_usize + kmer_length as usize];

                        // Check complexity constraint if threshold specified
                        if complexity_threshold == 0.0 {
                            true
                        } else {
                            let complexity = complexity_measure.calculate(kmer, kmer_length);
                            complexity >= complexity_threshold
                        }
                    })
                    .map(|(_pos, val)| val),
            );
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_minimizers() {
        // Simple sequence test
        let seq = b"ACGTACGTACGT";
        let k = 5;
        let w = 3;
        let hasher = KmerHasher::new(k as usize);
        let minimizers = compute_minimizers(seq, &hasher, k, w, 0.0);

        // We should have at least one minimizer
        assert!(!minimizers.is_empty());

        // Test with a sequence shorter than k
        let short_seq = b"ACGT";
        let short_minimizers = compute_minimizers(short_seq, &hasher, k, w, 0.0);
        assert!(short_minimizers.is_empty());
    }

    #[test]
    fn test_calculate_scaled_entropy() {
        // Test short k-mers (should return 1.0 for k < 10)
        let short_kmer = b"ACGT";
        let entropy = calculate_scaled_shannon_entropy(short_kmer, 8);
        assert_eq!(entropy, 1.0, "Expected 1.0 for k-mer length < 10");

        // Test minimum entropy (homopolymer, 10bp)
        let min_entropy_kmer = b"AAAAAAAAAA";
        let entropy = calculate_scaled_shannon_entropy(min_entropy_kmer, 10);
        assert!(entropy < 0.1, "Expected very low entropy, got {}", entropy);

        // Test moderate entropy (alternating pattern, 10bp)
        let alt_entropy_kmer = b"ATATATATAT";
        let entropy = calculate_scaled_shannon_entropy(alt_entropy_kmer, 10);
        assert!(
            (0.5..1.0).contains(&entropy),
            "Expected moderate entropy, got {}",
            entropy
        );

        // Test maximum entropy (diverse 10bp)
        let max_entropy_kmer = b"ACGTACGTAC";
        let entropy = calculate_scaled_shannon_entropy(max_entropy_kmer, 10);
        assert!(
            entropy > 0.9,
            "Expected high entropy for diverse 10-mer, got {}",
            entropy
        );

        // Test realistic k-mer (31bp, default k)
        let realistic_kmer = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let entropy = calculate_scaled_shannon_entropy(realistic_kmer, 31);
        assert!(
            entropy > 0.9,
            "Expected high entropy for diverse 31-mer, got {}",
            entropy
        );
    }

    #[test]
    fn test_31mer_entropy_range() {
        // Test various 31-mers with different entropy values to demonstrate the range

        // Homopolymer - lowest entropy (31 A's)
        let homopolymer = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let entropy = calculate_scaled_shannon_entropy(homopolymer, 31);
        assert!(entropy < 0.01, "Homopolymer entropy = {}", entropy);

        // Mostly one base with minimal variation - low entropy
        let mostly_a = b"AAAAAAAAAAACAAAAAGAAAAATAAAAAAA";
        let entropy = calculate_scaled_shannon_entropy(mostly_a, 31);
        assert!(
            (0.25..=0.35).contains(&entropy),
            "Mostly A entropy = {}",
            entropy
        );

        // GC alternating - moderate entropy (2 bases, equal distribution)
        let gc_alternating = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
        let entropy = calculate_scaled_shannon_entropy(gc_alternating, 31);
        assert!(
            (0.45..=0.55).contains(&entropy),
            "GC alternating entropy = {}",
            entropy
        );

        // AT with G ending - moderate entropy (mostly 2 bases)
        let dinuc_repeat = b"ATATATATATATATATATATATATATATATG";
        let entropy = calculate_scaled_shannon_entropy(dinuc_repeat, 31);
        assert!(
            (0.55..=0.65).contains(&entropy),
            "AT+G repeat entropy = {}",
            entropy
        );

        // Trinucleotide repeat - high entropy (ACG repeated)
        let trinuc_repeat = b"ACGACGACGACGACGACGACGACGACGACGA";
        let entropy = calculate_scaled_shannon_entropy(trinuc_repeat, 31);
        assert!(
            (0.75..=0.85).contains(&entropy),
            "ACG repeat entropy = {}",
            entropy
        );

        // Four bases uneven distribution - high entropy
        let four_uneven = b"ACGTACGTACGTAAAACCCGGGTTTACGTAC";
        let entropy = calculate_scaled_shannon_entropy(four_uneven, 31);
        assert!(
            (0.8..=1.0).contains(&entropy),
            "Four bases uneven entropy = {}",
            entropy
        );

        // Complex pattern with all 4 bases - very high entropy
        let complex_repeat = b"AACCGGTTAACCGGTTAACCGGTTAACCGGT";
        let entropy = calculate_scaled_shannon_entropy(complex_repeat, 31);
        assert!(entropy >= 0.95, "Complex pattern entropy = {}", entropy);

        // Four bases perfectly balanced - maximum entropy
        let four_balanced = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let entropy = calculate_scaled_shannon_entropy(four_balanced, 31);
        assert!(entropy >= 0.95, "Four bases balanced entropy = {}", entropy);

        // Verify entropy ordering makes sense
        let one_base = calculate_scaled_shannon_entropy(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 31);
        let two_bases_even =
            calculate_scaled_shannon_entropy(b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG", 31);
        let three_bases = calculate_scaled_shannon_entropy(b"ACGACGACGACGACGACGACGACGACGACGA", 31);
        let four_bases = calculate_scaled_shannon_entropy(b"ACGTACGTACGTACGTACGTACGTACGTACG", 31);

        // Entropy should increase with base diversity
        assert!(
            one_base < two_bases_even,
            "1 base ({}) < 2 bases even ({})",
            one_base,
            two_bases_even
        );
        assert!(
            two_bases_even < three_bases,
            "2 bases even ({}) < 3 bases ({})",
            two_bases_even,
            three_bases
        );
        assert!(
            three_bases < four_bases,
            "3 bases ({}) < 4 bases ({})",
            three_bases,
            four_bases
        );

        // Verify threshold behavior: common thresholds like 0.01 should filter appropriately
        assert!(
            one_base < 0.01,
            "Homopolymer should be filtered at 0.01 threshold"
        );
        assert!(
            four_bases > 0.01,
            "High diversity should pass 0.01 threshold"
        );
        assert!(four_bases > 0.5, "High diversity should pass 0.5 threshold");
    }

    #[test]
    fn test_near_homopolymer_entropy() {
        // 30 A's + 1 T, entropy ~0.1028
        let near_homopolymer = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT";
        let entropy = calculate_scaled_shannon_entropy(near_homopolymer, 31);
        assert!(entropy < 0.5, "Entropy {:.4} should be < 0.5", entropy);
        assert!(entropy < 0.15, "Entropy {:.4} should be < 0.15", entropy);
    }

    #[test]
    fn test_decode_minimizer_not_complement() {
        // Test decode_u64 returns original k-mer or its rc
        let test_kmer = b"GCTGAGAGCGGCTGTGGCCTCTGTCTGCTGC";
        let k = 31;
        let w = 15;

        // Pad for length
        let mut test_seq = Vec::new();
        test_seq.extend_from_slice(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); // 31 As before
        test_seq.extend_from_slice(test_kmer);
        test_seq.extend_from_slice(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); // 31 As after

        let hasher = KmerHasher::new(k as usize);
        let minimizers = compute_minimizers(&test_seq, &hasher, k as u8, w as u8, 0.0);

        assert!(!minimizers.is_empty(), "Should have at least one minimizer");

        let reverse_complement = |s: &str| -> String {
            s.chars()
                .rev()
                .map(|c| match c {
                    'A' => 'T',
                    'T' => 'A',
                    'G' => 'C',
                    'C' => 'G',
                    _ => c,
                })
                .collect()
        };

        let test_seq_str = String::from_utf8_lossy(&test_seq).to_string();

        // Check all decoded minimizers appear in test sequence as fwd or rc
        let mut found_valid = false;
        match &minimizers {
            crate::MinimizerVec::U64(vec) => {
                for &value in vec {
                    let decoded = String::from_utf8_lossy(&decode_u64(value, k as u8)).to_string();
                    let decoded_revcomp = reverse_complement(&decoded);

                    if test_seq_str.contains(&decoded) || test_seq_str.contains(&decoded_revcomp) {
                        found_valid = true;
                    } else {
                        panic!(
                            "Minimizer '{}' not found as forward or revcomp in test sequence",
                            decoded
                        );
                    }
                }
            }
            crate::MinimizerVec::U128(_) => panic!("Expected U64 for k=31"),
        }

        assert!(found_valid, "No valid decoded minimizers found");
    }

    #[test]
    fn test_decode_canonical_revcomp_smaller() {
        // Test k-mer where the rc is lexicographically smaller becomes canconical
        let test_kmer = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let k = 31;
        let w = 15;

        let mut test_seq = Vec::new();
        test_seq.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
        test_seq.extend_from_slice(test_kmer);
        test_seq.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");

        let hasher = KmerHasher::new(k as usize);
        let minimizers = compute_minimizers(&test_seq, &hasher, k as u8, w as u8, 0.0);

        assert!(!minimizers.is_empty(), "Should have at least one minimizer");

        let reverse_complement = |s: &str| -> String {
            s.chars()
                .rev()
                .map(|c| match c {
                    'A' => 'T',
                    'T' => 'A',
                    'G' => 'C',
                    'C' => 'G',
                    _ => c,
                })
                .collect()
        };

        let test_seq_str = String::from_utf8_lossy(&test_seq).to_string();

        match &minimizers {
            crate::MinimizerVec::U64(vec) => {
                for &value in vec {
                    let decoded = String::from_utf8_lossy(&decode_u64(value, k as u8)).to_string();
                    let decoded_revcomp = reverse_complement(&decoded);

                    assert!(
                        test_seq_str.contains(&decoded) || test_seq_str.contains(&decoded_revcomp),
                        "Decoded minimizer '{}' not found in test sequence as forward or reverse complement",
                        decoded
                    );
                }
            }
            crate::MinimizerVec::U128(_) => panic!("Expected U64 for k=31"),
        }
    }

    #[test]
    fn test_decode_edge_cases() {
        let k = 31;
        let w = 15;
        let hasher = KmerHasher::new(k as usize);

        let reverse_complement = |s: &str| -> String {
            s.chars()
                .rev()
                .map(|c| match c {
                    'A' => 'T',
                    'T' => 'A',
                    'G' => 'C',
                    'C' => 'G',
                    _ => c,
                })
                .collect()
        };

        // Test cases: (description, k-mer)
        let test_cases = vec![
            ("All A's", b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
            ("All C's", b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
            ("All G's", b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
            ("All T's", b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
            ("Near palindrome", b"ACGTACGTACGTACGTACGTACGTACGTACG"),
            ("AT repeat", b"ATATATATATATATATATATATATATATATA"),
            ("GC repeat", b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG"),
        ];

        for (desc, test_kmer) in test_cases {
            let mut test_seq = Vec::new();
            test_seq.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"); // Pad with N's
            test_seq.extend_from_slice(test_kmer);
            test_seq.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");

            let minimizers = compute_minimizers(&test_seq, &hasher, k as u8, w as u8, 0.0);

            if minimizers.is_empty() {
                continue; // Skip if no minimizers (e.g., all N's filtered)
            }

            let test_seq_str = String::from_utf8_lossy(&test_seq).to_string();

            match &minimizers {
                crate::MinimizerVec::U64(vec) => {
                    for &value in vec {
                        let decoded =
                            String::from_utf8_lossy(&decode_u64(value, k as u8)).to_string();
                        let decoded_revcomp = reverse_complement(&decoded);

                        assert!(
                            test_seq_str.contains(&decoded)
                                || test_seq_str.contains(&decoded_revcomp),
                            "{}: Decoded minimizer '{}' not found in test sequence",
                            desc,
                            decoded
                        );
                    }
                }
                crate::MinimizerVec::U128(_) => panic!("Expected U64 for k=31"),
            }
        }
    }

    #[test]
    fn test_decode_u128_long_kmer() {
        // Test long k-mers with u128
        let test_kmer = b"ACGTACGTACGTACGTACGTACGTACGTACGTA"; // 33bp
        let k = 33;
        let w = 17; // Odd window, l = 33 + 17 - 1 = 49 (odd)

        let mut test_seq = Vec::new();
        test_seq.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"); // Pad with T's instead of N's
        test_seq.extend_from_slice(test_kmer);
        test_seq.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

        let hasher = KmerHasher::new(k as usize);
        let minimizers = compute_minimizers(&test_seq, &hasher, k as u8, w as u8, 0.0);

        assert!(!minimizers.is_empty(), "Should have at least one minimizer");

        let reverse_complement = |s: &str| -> String {
            s.chars()
                .rev()
                .map(|c| match c {
                    'A' => 'T',
                    'T' => 'A',
                    'G' => 'C',
                    'C' => 'G',
                    _ => c,
                })
                .collect()
        };

        let test_seq_str = String::from_utf8_lossy(&test_seq).to_string();

        match &minimizers {
            crate::MinimizerVec::U128(vec) => {
                for &value in vec {
                    let decoded = String::from_utf8_lossy(&decode_u128(value, k as u8)).to_string();
                    let decoded_revcomp = reverse_complement(&decoded);

                    assert!(
                        test_seq_str.contains(&decoded) || test_seq_str.contains(&decoded_revcomp),
                        "Decoded u128 minimizer '{}' not found in test sequence as forward or reverse complement",
                        decoded
                    );
                }
            }
            crate::MinimizerVec::U64(_) => panic!("Expected U128 for k=33"),
        }
    }

    #[test]
    fn test_linguistic_complexity() {
        // Test homopolymer - should have very low complexity
        let homopolymer = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 31 A's
        let lc_homo = calculate_linguistic_complexity(homopolymer, 31);
        assert!(
            lc_homo < 0.0001,
            "Homopolymer LC should be near 0, got {}",
            lc_homo
        );

        // Test repeating pattern - should have low complexity
        let repeat = b"ACGTACGTACGTACGTACGTACGTACGTACG"; // 31bp repeating pattern
        let lc_repeat = calculate_linguistic_complexity(repeat, 31);
        assert!(
            lc_repeat < 0.001,
            "Repeating pattern LC should be low, got {}",
            lc_repeat
        );

        // Test more diverse sequence - should have higher complexity than homopolymer
        let diverse = b"ACGTAGCTAGCTAGTCGATCGTAGCTAGCTA"; // 31bp more random
        let lc_diverse = calculate_linguistic_complexity(diverse, 31);
        assert!(
            lc_diverse > 0.001,
            "Diverse sequence LC should be > 0.001, got {}",
            lc_diverse
        );
        assert!(lc_diverse > lc_homo, "Diverse should be > homopolymer");
        assert!(lc_diverse > lc_repeat, "Diverse should be > repeat");

        // Test k<10 always returns 1.0
        let short = b"ACGTACGT"; // 8bp
        let lc = calculate_linguistic_complexity(short, 8);
        assert_eq!(lc, 1.0, "Short k-mer (k<10) should return 1.0, got {}", lc);
    }

    #[test]
    fn test_sdust_complexity() {
        // Test bounds: min (homopolymer) and max (perfect) complexity
        let homopolymer = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 31 A's, all AAA triplets
        assert_eq!(calculate_sdust_score(homopolymer, 31), 0.0);

        // Perfect complexity 31mer discovered through simulation
        let perfect = b"TTATAACTTCGGCATCTCCCTGTCACAGAGT"; // All 29 triplets unique
        assert_eq!(calculate_sdust_score(perfect, 31), 1.0);

        // Test relative ordering of complexity for different repeat types
        let dinuc = b"ATATATATATATATATATATATATATATATAT";
        let trinuc = b"ACGACGACGACGACGACGACGACGACGACGA";
        let diverse = b"ACGTAGCTGATCGATCGTAGCTAGCTAGCTA";

        let c_homo = calculate_sdust_score(homopolymer, 31);
        let c_dinuc = calculate_sdust_score(&dinuc[..31], 31);
        let c_trinuc = calculate_sdust_score(trinuc, 31);
        let c_diverse = calculate_sdust_score(diverse, 31);

        assert!(c_homo < c_dinuc && c_dinuc < c_trinuc && c_trinuc < c_diverse);

        // Test edge cases
        assert_eq!(calculate_sdust_score(b"AC", 2), 1.0); // Too short for triplets
        assert_eq!(calculate_sdust_score(b"AAACCC", 6), 1.0); // All unique triplets
        assert_eq!(calculate_sdust_score(b"AAAAAA", 6), 0.0); // All same triplet
    }

    #[test]
    fn test_complexity_measure_enum() {
        // Test from_str parsing
        assert_eq!(
            ComplexityMeasure::from_str("scaled-shannon").unwrap(),
            ComplexityMeasure::ScaledShannon
        );
        assert_eq!(
            ComplexityMeasure::from_str("linguistic").unwrap(),
            ComplexityMeasure::Linguistic
        );
        assert_eq!(
            ComplexityMeasure::from_str("sdust").unwrap(),
            ComplexityMeasure::SDust
        );
        assert!(ComplexityMeasure::from_str("invalid").is_err());

        // Test as_str
        assert_eq!(ComplexityMeasure::ScaledShannon.as_str(), "scaled-shannon");
        assert_eq!(ComplexityMeasure::Linguistic.as_str(), "linguistic");
        assert_eq!(ComplexityMeasure::SDust.as_str(), "sdust");

        // Test calculate method dispatch
        let test_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 31 A's
        let shannon = ComplexityMeasure::ScaledShannon.calculate(test_seq, 31);
        let linguistic = ComplexityMeasure::Linguistic.calculate(test_seq, 31);
        let sdust = ComplexityMeasure::SDust.calculate(test_seq, 31);

        // Verify all measures agree that homopolymer has low complexity
        assert!(
            shannon < 0.1,
            "Shannon entropy of homopolymer should be low"
        );
        assert!(
            linguistic < 0.0001,
            "Linguistic complexity of homopolymer should be near 0"
        );
        assert!(
            sdust < 0.01,
            "SDust complexity of homopolymer should be near 0"
        );
    }
}
