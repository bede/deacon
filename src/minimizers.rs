use packed_seq::AsciiSeq;
use xxhash_rust::xxh3;

pub const DEFAULT_KMER_LENGTH: u8 = 31;
pub const DEFAULT_WINDOW_SIZE: u8 = 15;

/// Check if nucleotide is valid ACGT (case insensitive)
#[inline]
fn is_valid_acgt(nucleotide: u8) -> bool {
    matches!(
        nucleotide,
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
    )
}

/// Check if k-mer contains only ACGT nucleotides
#[inline]
fn kmer_contains_only_acgt(kmer: &[u8]) -> bool {
    kmer.iter().all(|&b| is_valid_acgt(b))
}

/// Canonicalise IUPAC ambiguous nucleotides to ACGT
#[inline]
fn canonicalise_nucleotide(nucleotide: u8) -> u8 {
    match nucleotide {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' => b'T',
        b'R' | b'r' => b'G',
        b'Y' | b'y' => b'C',
        b'S' | b's' => b'G',
        b'W' | b'w' => b'A',
        b'K' | b'k' => b'G',
        b'M' | b'm' => b'C',
        b'B' | b'b' => b'C',
        b'D' | b'd' => b'G',
        b'H' | b'h' => b'C',
        b'V' | b'v' => b'G',
        b'N' | b'n' => b'C',
        _ => b'C',
    }
}

/// Canonicalise a sequence
fn canonicalise_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&nucleotide| canonicalise_nucleotide(nucleotide))
        .collect()
}

/// Returns vector of all minimizer hashes for a sequence
pub fn compute_minimizer_hashes(
    seq: &[u8],
    kmer_length: u8,
    window_size: u8,
    entropy_threshold: f32,
) -> Vec<u64> {
    let mut hashes = Vec::new();
    fill_minimizer_hashes(
        seq,
        kmer_length,
        window_size,
        &mut hashes,
        entropy_threshold,
    );
    hashes
}

/// Calculate scaled entropy using character frequency analysis
/// Returns scaled entropy between 0.0 and 1.0
#[inline]
fn calculate_scaled_entropy(kmer: &[u8], kmer_length: u8) -> f32 {
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

/// Fill a vector with minimizer hashes, skipping k-mers with non-ACGT bases
/// and optionally filtering by scaled entropy
pub fn fill_minimizer_hashes(
    seq: &[u8],
    kmer_length: u8,
    window_size: u8,
    hashes: &mut Vec<u64>,
    entropy_threshold: f32,
) {
    hashes.clear();

    // Skip if sequence is too short
    if seq.len() < kmer_length as usize {
        return;
    }

    let canonical_seq = canonicalise_sequence(seq);

    // Get minimizer positions using simd-minimizers
    let mut positions = Vec::new();
    simd_minimizers::canonical_minimizer_positions(
        AsciiSeq(&canonical_seq),
        kmer_length as usize,
        window_size as usize,
        &mut positions,
    );

    // Filter positions to only include k-mers with ACGT bases and sufficient entropy
    let valid_positions: Vec<u32> = positions
        .into_iter()
        .filter(|&pos| {
            let pos_usize = pos as usize;
            let kmer = &seq[pos_usize..pos_usize + kmer_length as usize];

            // Check ACGT constraint
            if !kmer_contains_only_acgt(kmer) {
                return false;
            }

            // Check scaled entropy constraint if threshold specified
            if entropy_threshold == 0.0 {
                true
            } else {
                let entropy = calculate_scaled_entropy(kmer, kmer_length);
                entropy >= entropy_threshold
            }
        })
        .collect();

    if kmer_length > 32 {
        hashes.extend(
            simd_minimizers::iter_canonical_minimizer_values_u128(
                AsciiSeq(&canonical_seq),
                kmer_length as usize,
                &valid_positions,
            )
            .map(|kmer| xxh3::xxh3_64(&kmer.to_le_bytes())),
        );
    } else {
        hashes.extend(
            simd_minimizers::iter_canonical_minimizer_values(
                AsciiSeq(&canonical_seq),
                kmer_length as usize,
                &valid_positions,
            )
            .map(|kmer| xxh3::xxh3_64(&kmer.to_le_bytes())),
        );
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonicalise_nucleotide() {
        // Test basic nucleotides
        assert_eq!(canonicalise_nucleotide(b'A'), b'A');
        assert_eq!(canonicalise_nucleotide(b'C'), b'C');
        assert_eq!(canonicalise_nucleotide(b'G'), b'G');
        assert_eq!(canonicalise_nucleotide(b'T'), b'T');

        // Test lowercase
        assert_eq!(canonicalise_nucleotide(b'a'), b'A');
        assert_eq!(canonicalise_nucleotide(b'c'), b'C');

        // Test ambiguous nucleotides
        assert_eq!(canonicalise_nucleotide(b'R'), b'G');
        assert_eq!(canonicalise_nucleotide(b'Y'), b'C');
        assert_eq!(canonicalise_nucleotide(b'S'), b'G');
        assert_eq!(canonicalise_nucleotide(b'W'), b'A');
        assert_eq!(canonicalise_nucleotide(b'K'), b'G');
        assert_eq!(canonicalise_nucleotide(b'M'), b'C');
        assert_eq!(canonicalise_nucleotide(b'B'), b'C');
        assert_eq!(canonicalise_nucleotide(b'D'), b'G');
        assert_eq!(canonicalise_nucleotide(b'H'), b'C');
        assert_eq!(canonicalise_nucleotide(b'V'), b'G');
        assert_eq!(canonicalise_nucleotide(b'N'), b'C');
    }

    #[test]
    fn test_canonicalise_sequence() {
        let seq = b"ACGTN";
        let canonical = canonicalise_sequence(seq);
        assert_eq!(canonical, b"ACGTC");

        let seq = b"RYMKWS";
        let canonical = canonicalise_sequence(seq);
        assert_eq!(canonical, b"GCCGAG");
    }

    #[test]
    fn test_compute_minimizer_hashes() {
        // Simple sequence test
        let seq = b"ACGTACGTACGT";
        let k = 5;
        let w = 3;

        let hashes = compute_minimizer_hashes(seq, k, w, 0.0);

        // We should have at least one minimizer hash
        assert!(!hashes.is_empty());

        // Test with a sequence shorter than k
        let short_seq = b"ACGT";
        let short_hashes = compute_minimizer_hashes(short_seq, k, w, 0.0);
        assert!(short_hashes.is_empty());
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
        assert!(entropy < 0.1, "Expected very low entropy, got {entropy}");

        // Test moderate entropy (alternating pattern, 10bp)
        let alt_entropy_kmer = b"ATATATATAT";
        let entropy = calculate_scaled_entropy(alt_entropy_kmer, 10);
        assert!(
            (0.5..1.0).contains(&entropy),
            "Expected moderate entropy, got {entropy}"
        );

        // Test maximum entropy (diverse 10bp)
        let max_entropy_kmer = b"ACGTACGTAC";
        let entropy = calculate_scaled_entropy(max_entropy_kmer, 10);
        assert!(
            entropy > 0.9,
            "Expected high entropy for diverse 10-mer, got {entropy}"
        );

        // Test realistic k-mer (31bp, default k)
        let realistic_kmer = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let entropy = calculate_scaled_entropy(realistic_kmer, 31);
        assert!(
            entropy > 0.9,
            "Expected high entropy for diverse 31-mer, got {entropy}"
        );
    }

    #[test]
    fn test_31mer_entropy_range() {
        // Test various 31-mers with different entropy values to demonstrate the range

        // Homopolymer - lowest entropy (31 A's)
        let homopolymer = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let entropy = calculate_scaled_entropy(homopolymer, 31);
        assert!(entropy < 0.01, "Homopolymer entropy = {entropy}");

        // Mostly one base with minimal variation - low entropy
        let mostly_a = b"AAAAAAAAAAACAAAAAGAAAAATAAAAAAA";
        let entropy = calculate_scaled_entropy(mostly_a, 31);
        assert!(
            (0.25..=0.35).contains(&entropy),
            "Mostly A entropy = {entropy}"
        );

        // GC alternating - moderate entropy (2 bases, equal distribution)
        let gc_alternating = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
        let entropy = calculate_scaled_entropy(gc_alternating, 31);
        assert!(
            (0.45..=0.55).contains(&entropy),
            "GC alternating entropy = {entropy}"
        );

        // AT with G ending - moderate entropy (mostly 2 bases)
        let dinuc_repeat = b"ATATATATATATATATATATATATATATATG";
        let entropy = calculate_scaled_entropy(dinuc_repeat, 31);
        assert!(
            (0.55..=0.65).contains(&entropy),
            "AT+G repeat entropy = {entropy}"
        );

        // Trinucleotide repeat - high entropy (ACG repeated)
        let trinuc_repeat = b"ACGACGACGACGACGACGACGACGACGACGA";
        let entropy = calculate_scaled_entropy(trinuc_repeat, 31);
        assert!(
            (0.75..=0.85).contains(&entropy),
            "ACG repeat entropy = {entropy}"
        );

        // Four bases uneven distribution - high entropy
        let four_uneven = b"ACGTACGTACGTAAAACCCGGGTTTACGTAC";
        let entropy = calculate_scaled_entropy(four_uneven, 31);
        assert!(
            (0.8..=1.0).contains(&entropy),
            "Four bases uneven entropy = {entropy}"
        );

        // Complex pattern with all 4 bases - very high entropy
        let complex_repeat = b"AACCGGTTAACCGGTTAACCGGTTAACCGGT";
        let entropy = calculate_scaled_entropy(complex_repeat, 31);
        assert!(entropy >= 0.95, "Complex pattern entropy = {entropy}");

        // Four bases perfectly balanced - maximum entropy
        let four_balanced = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        let entropy = calculate_scaled_entropy(four_balanced, 31);
        assert!(entropy >= 0.95, "Four bases balanced entropy = {entropy}");

        // Verify entropy ordering makes sense
        let one_base = calculate_scaled_entropy(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 31);
        let two_bases_even = calculate_scaled_entropy(b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG", 31);
        let three_bases = calculate_scaled_entropy(b"ACGACGACGACGACGACGACGACGACGACGA", 31);
        let four_bases = calculate_scaled_entropy(b"ACGTACGTACGTACGTACGTACGTACGTACG", 31);

        // Entropy should increase with base diversity
        assert!(
            one_base < two_bases_even,
            "1 base ({one_base}) < 2 bases even ({two_bases_even})"
        );
        assert!(
            two_bases_even < three_bases,
            "2 bases even ({two_bases_even}) < 3 bases ({three_bases})"
        );
        assert!(
            three_bases < four_bases,
            "3 bases ({three_bases}) < 4 bases ({four_bases})"
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
        assert!(entropy < 0.5, "Entropy {entropy:.4} should be < 0.5");
        assert!(entropy < 0.15, "Entropy {entropy:.4} should be < 0.15");
    }
}
