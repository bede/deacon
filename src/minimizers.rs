use packed_seq::{PackedNSeqVec, SeqVec, u32x8, unpack_base};

pub const DEFAULT_KMER_LENGTH: u8 = 31;
pub const DEFAULT_WINDOW_SIZE: u8 = 15;

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
type KmerHasher = simd_minimizers::seq_hash::NtHasher<true, 1>;

/// Require positive odd k/w, k <= 61 and k+w <= 96.
pub(crate) fn validate_k_w(kmer_length: u8, window_size: u8) -> anyhow::Result<()> {
    let (k, w) = (kmer_length as usize, window_size as usize);
    if k == 0 || w == 0 || k > 61 || k + w > 96 || k.is_multiple_of(2) || !(k + w).is_multiple_of(2)
    {
        anyhow::bail!(
            "Invalid k-w combination: k={}, w={}, k+w={} (constraints: k and w positive and odd, k<=61, k+w<=96)",
            k,
            w,
            k + w
        );
    }
    Ok(())
}

/// Canonical minimizers for fixed k/w with reusable buffers.
///
/// Values use two bits per base: u64 for k <= 32, otherwise u128.
/// Accepts lowercase and skips windows containing non-ACGT bases.
#[derive(Clone)]
pub struct Minimizers {
    params: Parameters,
    buffers: Buffers,
}

#[derive(Clone)]
struct Parameters {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
}

#[derive(Clone)]
struct Buffers {
    packed_nseq: PackedNSeqVec,
    positions: Vec<u32>,
    values: MinimizerVec,
    cache: (simd_minimizers::Cache, Vec<u32x8>, Vec<u32x8>),
}

impl Minimizers {
    pub fn new(kmer_length: u8, window_size: u8) -> anyhow::Result<Self> {
        validate_k_w(kmer_length, window_size)?;
        Ok(Self {
            params: Parameters {
                kmer_length,
                window_size,
                hasher: KmerHasher::new(kmer_length as usize),
            },
            buffers: Buffers {
                packed_nseq: PackedNSeqVec {
                    seq: Default::default(),
                    ambiguous: Default::default(),
                },
                positions: Vec::new(),
                values: if kmer_length <= 32 {
                    MinimizerVec::U64(Vec::new())
                } else {
                    MinimizerVec::U128(Vec::new())
                },
                cache: Default::default(),
            },
        })
    }

    pub fn kmer_length(&self) -> u8 {
        self.params.kmer_length
    }

    pub fn window_size(&self) -> u8 {
        self.params.window_size
    }

    /// Minimizers in sequence order. Empty for sequences shorter than k.
    pub fn compute(&mut self, seq: &[u8]) -> &MinimizerVec {
        self.clear();
        self.extend(seq);
        &self.buffers.values
    }

    /// Append minimizers, pooling mates.
    #[inline]
    pub(crate) fn extend(&mut self, seq: &[u8]) {
        if seq.len() >= self.params.kmer_length as usize {
            self.params.extend(seq, &mut self.buffers);
        }
    }

    #[inline]
    pub(crate) fn clear(&mut self) {
        self.buffers.values.clear();
    }

    #[inline]
    pub(crate) fn values(&self) -> &MinimizerVec {
        &self.buffers.values
    }
}

impl Parameters {
    #[inline]
    fn extend(&self, seq: &[u8], buffers: &mut Buffers) {
        let Buffers {
            packed_nseq,
            positions,
            values,
            cache,
        } = buffers;
        packed_nseq.seq.clear();
        packed_nseq.ambiguous.clear();
        positions.clear();
        packed_nseq.seq.push_ascii(seq);
        packed_nseq.ambiguous.push_ascii(seq);

        let out = simd_minimizers::canonical_minimizers(
            self.kmer_length as usize,
            self.window_size as usize,
        )
        .hasher(&self.hasher)
        .run_skip_ambiguous_windows_with_buf(packed_nseq.as_slice(), positions, cache);

        match values {
            MinimizerVec::U64(vec) => vec.extend(out.pos_and_values_u64().map(|(_pos, val)| val)),
            MinimizerVec::U128(vec) => vec.extend(out.pos_and_values_u128().map(|(_pos, val)| val)),
        }
    }
}

/// Returns scaled entropy between 0.0 and 1.0
#[inline]
pub(crate) fn calculate_scaled_entropy(kmer: &[u8], kmer_length: u8) -> f32 {
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

/// kdust: max-normalised DUST triplet score in [0,1] of a packed minimizer
#[inline]
pub(crate) fn calculate_kdust(code: u128, kmer_length: u8) -> f32 {
    // k=3 divides by zero, k<3 has no triplets
    if kmer_length < 4 {
        return 1.0;
    }
    let k = kmer_length as usize;
    let mut counts = [0u8; 64];
    let mut score = 0u32;
    let mut tri = 0usize;
    for i in 0..k {
        tri = ((tri << 2) | ((code >> (2 * i)) & 0b11) as usize) & 0b11_1111;
        if i >= 2 {
            score += counts[tri] as u32;
            counts[tri] += 1;
        }
    }
    let l = (k - 2) as f32;
    1.0 - score as f32 / (l * (l - 1.0) / 2.0)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[serde(rename_all = "lowercase")]
pub enum ComplexityAlgorithm {
    #[default]
    Kdust,
    Shannon,
}

impl std::fmt::Display for ComplexityAlgorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            ComplexityAlgorithm::Kdust => "kdust",
            ComplexityAlgorithm::Shannon => "shannon",
        })
    }
}

/// Packed minimizers: u64 for k <= 32, otherwise u128.
#[derive(Clone)]
pub enum MinimizerVec {
    U64(Vec<u64>),
    U128(Vec<u128>),
}

impl MinimizerVec {
    pub fn clear(&mut self) {
        match self {
            MinimizerVec::U64(v) => v.clear(),
            MinimizerVec::U128(v) => v.clear(),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            MinimizerVec::U64(v) => v.len(),
            MinimizerVec::U128(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            MinimizerVec::U64(v) => v.is_empty(),
            MinimizerVec::U128(v) => v.is_empty(),
        }
    }
}

/// Packed minimizer shards for building and writing indexes.
#[derive(Clone)]
pub(crate) enum MinimizerShards {
    U64(Vec<Vec<u64>>),
    U128(Vec<Vec<u128>>),
}

impl MinimizerShards {
    pub fn len(&self) -> usize {
        match self {
            MinimizerShards::U64(v) => v.iter().map(|inner| inner.len()).sum(),
            MinimizerShards::U128(v) => v.iter().map(|inner| inner.len()).sum(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compute_minimizers(seq: &[u8], k: u8, w: u8) -> MinimizerVec {
        Minimizers::new(k, w).unwrap().compute(seq).clone()
    }

    #[test]
    fn test_compute_minimizers() {
        let mut minimizers = Minimizers::new(5, 3).unwrap();
        assert!(!minimizers.compute(b"ACGTACGTACGT").is_empty());
        // Reuse buffers for a short sequence.
        assert!(minimizers.compute(b"ACGT").is_empty());
        assert!(Minimizers::new(4, 3).is_err());
    }

    #[test]
    fn non_acgt_mapped_to_n() {
        let minimizers = |base| {
            let mut seq = b"ACGT".repeat(40);
            seq[64] = base;
            match compute_minimizers(&seq, 31, 15) {
                MinimizerVec::U64(vec) => vec,
                MinimizerVec::U128(_) => unreachable!(),
            }
        };
        let expected = minimizers(b'N');

        for &base in b"nRrYySsWwKkMmBbDdHhVvUu-" {
            assert_eq!(minimizers(base), expected, "base {}", base as char);
        }
    }

    #[test]
    fn test_kdust_short_kmers_are_not_nan() {
        // k=3 divides by zero without the guard
        for k in [1u8, 2, 3] {
            let score = calculate_kdust(0, k);
            assert!(score.is_finite(), "k={k} scored {score}");
            assert_eq!(score, 1.0);
        }
        // k=4 is smallest valid k
        assert_eq!(calculate_kdust(0, 4), 0.0); // AAAA
    }

    #[test]
    fn test_calculate_scaled_entropy() {
        // Test short k-mers (should return 1.0 for k < 10)
        let short_kmer = b"ACGT";
        let entropy = calculate_scaled_entropy(short_kmer, 8);
        assert_eq!(entropy, 1.0, "Expected 1.0 for k-mer length < 10");

        // Test minimum entropy (homopolymer, 10bp)
        let min_entropy_kmer = b"AAAAAAAAAA";
        let entropy = calculate_scaled_entropy(min_entropy_kmer, 10);
        assert!(entropy < 0.1, "Expected very low entropy, got {}", entropy);

        // Test moderate entropy (alternating pattern, 10bp)
        let alt_entropy_kmer = b"ATATATATAT";
        let entropy = calculate_scaled_entropy(alt_entropy_kmer, 10);
        assert!(
            (0.5..1.0).contains(&entropy),
            "Expected moderate entropy, got {}",
            entropy
        );

        // Test maximum entropy (diverse 10bp)
        let max_entropy_kmer = b"ACGTACGTAC";
        let entropy = calculate_scaled_entropy(max_entropy_kmer, 10);
        assert!(
            entropy > 0.9,
            "Expected high entropy for diverse 10-mer, got {}",
            entropy
        );

        // Test realistic k-mer (31bp, default k)
        let realistic_kmer = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let entropy = calculate_scaled_entropy(realistic_kmer, 31);
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
        let entropy = calculate_scaled_entropy(homopolymer, 31);
        assert!(entropy < 0.01, "Homopolymer entropy = {}", entropy);

        // Mostly one base with minimal variation - low entropy
        let mostly_a = b"AAAAAAAAAAACAAAAAGAAAAATAAAAAAA";
        let entropy = calculate_scaled_entropy(mostly_a, 31);
        assert!(
            (0.25..=0.35).contains(&entropy),
            "Mostly A entropy = {}",
            entropy
        );

        // GC alternating - moderate entropy (2 bases, equal distribution)
        let gc_alternating = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
        let entropy = calculate_scaled_entropy(gc_alternating, 31);
        assert!(
            (0.45..=0.55).contains(&entropy),
            "GC alternating entropy = {}",
            entropy
        );

        // AT with G ending - moderate entropy (mostly 2 bases)
        let dinuc_repeat = b"ATATATATATATATATATATATATATATATG";
        let entropy = calculate_scaled_entropy(dinuc_repeat, 31);
        assert!(
            (0.55..=0.65).contains(&entropy),
            "AT+G repeat entropy = {}",
            entropy
        );

        // Trinucleotide repeat - high entropy (ACG repeated)
        let trinuc_repeat = b"ACGACGACGACGACGACGACGACGACGACGA";
        let entropy = calculate_scaled_entropy(trinuc_repeat, 31);
        assert!(
            (0.75..=0.85).contains(&entropy),
            "ACG repeat entropy = {}",
            entropy
        );

        // Four bases uneven distribution - high entropy
        let four_uneven = b"ACGTACGTACGTAAAACCCGGGTTTACGTAC";
        let entropy = calculate_scaled_entropy(four_uneven, 31);
        assert!(
            (0.8..=1.0).contains(&entropy),
            "Four bases uneven entropy = {}",
            entropy
        );

        // Complex pattern with all 4 bases - very high entropy
        let complex_repeat = b"AACCGGTTAACCGGTTAACCGGTTAACCGGT";
        let entropy = calculate_scaled_entropy(complex_repeat, 31);
        assert!(entropy >= 0.95, "Complex pattern entropy = {}", entropy);

        // Four bases perfectly balanced - maximum entropy
        let four_balanced = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let entropy = calculate_scaled_entropy(four_balanced, 31);
        assert!(entropy >= 0.95, "Four bases balanced entropy = {}", entropy);

        // Verify entropy ordering makes sense
        let one_base = calculate_scaled_entropy(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 31);
        let two_bases_even = calculate_scaled_entropy(b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG", 31);
        let three_bases = calculate_scaled_entropy(b"ACGACGACGACGACGACGACGACGACGACGA", 31);
        let four_bases = calculate_scaled_entropy(b"ACGTACGTACGTACGTACGTACGTACGTACG", 31);

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
        let entropy = calculate_scaled_entropy(near_homopolymer, 31);
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

        let minimizers = compute_minimizers(&test_seq, k as u8, w as u8);

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
            MinimizerVec::U64(vec) => {
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
            MinimizerVec::U128(_) => panic!("Expected U64 for k=31"),
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

        let minimizers = compute_minimizers(&test_seq, k as u8, w as u8);

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
            MinimizerVec::U64(vec) => {
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
            MinimizerVec::U128(_) => panic!("Expected U64 for k=31"),
        }
    }

    #[test]
    fn test_decode_edge_cases() {
        let k = 31;
        let w = 15;

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

            let minimizers = compute_minimizers(&test_seq, k as u8, w as u8);

            if minimizers.is_empty() {
                continue; // Skip if no minimizers (e.g., all N's filtered)
            }

            let test_seq_str = String::from_utf8_lossy(&test_seq).to_string();

            match &minimizers {
                MinimizerVec::U64(vec) => {
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
                MinimizerVec::U128(_) => panic!("Expected U64 for k=31"),
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

        let minimizers = compute_minimizers(&test_seq, k as u8, w as u8);

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
            MinimizerVec::U128(vec) => {
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
            MinimizerVec::U64(_) => panic!("Expected U128 for k=33"),
        }
    }
}
