use packed_seq::AsciiSeq;
use simd_minimizers::scalar::canonical_minimizer_positions_scalar;
use xxhash_rust::xxh3;

pub const DEFAULT_KMER_LENGTH: usize = 31;
pub const DEFAULT_WINDOW_SIZE: usize = 15;
pub const SIMD_LEN_THRESHOLD: usize = 1000;

/// Canonicalize IUPAC ambiguous nucleotides to ACGT based on position, avoiding compositional bias?
#[inline]
fn canonicalize_nucleotide(nucleotide: u8, position: usize) -> u8 {
    match nucleotide {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' => b'T',
        b'R' | b'r' => {
            if position % 2 == 0 {
                b'A'
            } else {
                b'G'
            }
        }
        b'Y' | b'y' => {
            if position % 2 == 0 {
                b'C'
            } else {
                b'T'
            }
        }
        b'S' | b's' => {
            if position % 2 == 0 {
                b'G'
            } else {
                b'C'
            }
        }
        b'W' | b'w' => {
            if position % 2 == 0 {
                b'A'
            } else {
                b'T'
            }
        }
        b'K' | b'k' => {
            if position % 2 == 0 {
                b'G'
            } else {
                b'T'
            }
        }
        b'M' | b'm' => {
            if position % 2 == 0 {
                b'A'
            } else {
                b'C'
            }
        }
        b'B' | b'b' => match position % 3 {
            0 => b'C',
            1 => b'G',
            _ => b'T',
        },
        b'D' | b'd' => match position % 3 {
            0 => b'A',
            1 => b'G',
            _ => b'T',
        },
        b'H' | b'h' => match position % 3 {
            0 => b'A',
            1 => b'C',
            _ => b'T',
        },
        b'V' | b'v' => match position % 3 {
            0 => b'A',
            1 => b'C',
            _ => b'G',
        },
        b'N' | b'n' => match position % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        },
        _ => b'A', // Default to A for anything else
    }
}

/// Canonicalize a sequence
fn canonicalize_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .enumerate()
        .map(|(i, &nucleotide)| canonicalize_nucleotide(nucleotide, i))
        .collect()
}

/// Returns vector of all minimizer hashes for a sequence
pub fn compute_minimizer_hashes(seq: &[u8], kmer_length: usize, window_size: usize) -> Vec<u64> {
    // Skip if sequence is short
    if seq.len() < kmer_length {
        return Vec::new();
    }
    let canonical_seq = canonicalize_sequence(seq);
    // Get minimizer positions using simd-minimizers or scalar fallback
    let mut positions = Vec::with_capacity(seq.len() / window_size + 1);
    if canonical_seq.len() >= kmer_length + SIMD_LEN_THRESHOLD {
        simd_minimizers::canonical_minimizer_positions(
            AsciiSeq(&canonical_seq),
            kmer_length,
            window_size,
            &mut positions,
        );
    } else {
        canonical_minimizer_positions_scalar(
            AsciiSeq(&canonical_seq),
            kmer_length,
            window_size,
            &mut positions,
        );
    }
    // Convert positions to hash values with xxh3_64
    positions
        .iter()
        .map(|&pos| {
            let pos = pos as usize;
            if pos + kmer_length <= canonical_seq.len() {
                hash_canonical_kmer(&canonical_seq[pos..(pos + kmer_length)])
            } else {
                0 // Should never happen with valid minimizer positions
            }
        })
        .collect()
}

/// Fill a vector with minimizer hashes
pub fn fill_minimizer_hashes(
    seq: &[u8],
    kmer_length: usize,
    window_size: usize,
    hashes: &mut Vec<u64>,
) {
    hashes.clear();

    // Skip  if sequence is too short
    if seq.len() < kmer_length {
        return;
    }

    let canonical_seq = canonicalize_sequence(seq);

    // Get minimizer positions using simd-minimizers or scalar fallback
    let mut positions = Vec::new();
    if canonical_seq.len() >= kmer_length + SIMD_LEN_THRESHOLD {
        simd_minimizers::canonical_minimizer_positions(
            AsciiSeq(&canonical_seq),
            kmer_length,
            window_size,
            &mut positions,
        );
    } else {
        canonical_minimizer_positions_scalar(
            AsciiSeq(&canonical_seq),
            kmer_length,
            window_size,
            &mut positions,
        );
    }

    // Convert positions to hash values using xxh3_64
    for &pos in &positions {
        let pos = pos as usize;
        if pos + kmer_length <= canonical_seq.len() {
            hashes.push(hash_canonical_kmer(
                &canonical_seq[pos..(pos + kmer_length)],
            ));
        }
    }
}

/// Hash forward k-mers
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    xxh3::xxh3_64(kmer)
}

/// Hash function for canonical k-mers (min of forward and reverse complement)
#[inline]
fn hash_canonical_kmer(kmer: &[u8]) -> u64 {
    let forward_hash = hash_kmer(kmer);
    let reverse_hash = hash_reverse_complement(kmer);

    // Return min of both hashes
    std::cmp::min(forward_hash, reverse_hash)
}

/// Hash revcomp k-mers
fn hash_reverse_complement(kmer: &[u8]) -> u64 {
    let mut state = xxh3::Xxh3::new();

    for &nucleotide in kmer.iter().rev() {
        let complement = match nucleotide {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        };
        state.update(&[complement]);
    }

    state.digest()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonicalize_nucleotide() {
        // Test basic nucleotides
        assert_eq!(canonicalize_nucleotide(b'A', 0), b'A');
        assert_eq!(canonicalize_nucleotide(b'C', 0), b'C');
        assert_eq!(canonicalize_nucleotide(b'G', 0), b'G');
        assert_eq!(canonicalize_nucleotide(b'T', 0), b'T');

        // Test lowercase
        assert_eq!(canonicalize_nucleotide(b'a', 0), b'A');
        assert_eq!(canonicalize_nucleotide(b'c', 0), b'C');

        // Test ambiguous nucleotides (position-dependent)
        assert_eq!(canonicalize_nucleotide(b'R', 0), b'A'); // position 0 -> A
        assert_eq!(canonicalize_nucleotide(b'R', 1), b'G'); // position 1 -> G

        // Test N at different positions
        assert_eq!(canonicalize_nucleotide(b'N', 0), b'A');
        assert_eq!(canonicalize_nucleotide(b'N', 1), b'C');
        assert_eq!(canonicalize_nucleotide(b'N', 2), b'G');
        assert_eq!(canonicalize_nucleotide(b'N', 3), b'T');
    }

    #[test]
    fn test_canonicalize_sequence() {
        let seq = b"ACGTN";
        let canonical = canonicalize_sequence(seq);
        assert_eq!(canonical, b"ACGTA"); // N at position 4 becomes A

        let seq = b"RYMKWS";
        let canonical = canonicalize_sequence(seq);
        assert_eq!(canonical[0], b'A'); // R at position 0 becomes A
        assert_eq!(canonical[1], b'T'); // Y at position 1 becomes T (not C) since position 1 % 2 == 1
    }

    #[test]
    fn test_compute_minimizer_hashes() {
        // Simple sequence test
        let seq = b"ACGTACGTACGT";
        let k = 5;
        let w = 3;

        let hashes = compute_minimizer_hashes(seq, k, w);

        // We should have at least one minimizer hash
        assert!(!hashes.is_empty());

        // Test with a sequence shorter than k
        let short_seq = b"ACGT";
        let short_hashes = compute_minimizer_hashes(short_seq, k, w);
        assert!(short_hashes.is_empty());
    }

    #[test]
    fn test_hash_canonical_kmer() {
        // Test canonical hashing (min of forward/reverse)
        let kmer = b"ACGTA";
        let rc_kmer = b"TACGT";

        let forward_hash = hash_kmer(kmer);
        let reverse_hash = hash_kmer(rc_kmer);
        let canonical_hash = hash_canonical_kmer(kmer);

        // Canonical hash should be minimum of the two
        assert_eq!(canonical_hash, std::cmp::min(forward_hash, reverse_hash));
    }
}
