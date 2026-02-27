# packed-seq

[![crates.io](https://img.shields.io/crates/v/packed-seq)](https://crates.io/crates/packed-seq)
[![docs](https://img.shields.io/docsrs/packed-seq)](https://docs.rs/packed-seq)

A library for constructing and iterating packed `PackedSeq` DNA sequences that
handles the encoding, decoding, and complements of packed bases.

ASCII `ACTG` representations of DNA and general ASCII text are also supported via the `Seq` trait.

Fast SIMD-based iteration over sequences is supported by splitting the sequence
into 8 (slightly overlapping) chunks and iterating those in parallel in a
memory-efficient way.

**Paper:**
Please cite the
[`simd-minimizers`](https://github.com/rust-seq/simd-minimizers) paper, for which this
crate was developed:

- SimdMinimizers: Computing Random Minimizers, fast.  
  Ragnar Groot Koerkamp, Igor Martayan.
  SEA 2025 [https://doi.org/10.4230/LIPIcs.SEA.2025.20](doi.org/10.4230/LIPIcs.SEA.2025.20)

## Requirements

This library requires AVX2 or NEON instruction sets, which, on x64, requires
either `target-cpu=native` or `target-cpu=x86-64-v3`.
See [this README](https://github.com/ragnargrootkoerkamp/ensure_simd) for details and [this
blog](https://curiouscoding.nl/posts/distributing-rust-simd-binaries/) for background.
The same restrictions apply when using packed-seq in a larger project.

``` sh
RUSTFLAGS="-C target-cpu=native" cargo run --release
```

## Usage example

Full documentation can be found on [docs.rs](https://docs.rs/packed-seq).

```rust
use packed_seq::{SeqVec, Seq, AsciiSeqVec, PackedSeqVec, pack_char};
// Plain ASCII sequence.
let seq = b"ACTGCAGCGCATATGTAGT";
// ASCII DNA sequence.
let ascii_seq = AsciiSeqVec::from_ascii(seq);
// Packed DNA sequence.
let packed_seq = PackedSeqVec::from_ascii(seq);
assert_eq!(ascii_seq.len(), packed_seq.len());
// Iterate the ASCII characters.
let characters: Vec<u8> = seq.iter_bp().collect();
assert_eq!(characters, seq);

// Iterate the bases with 0..4 values.
let bases: Vec<u8> = seq.iter().copied().map(pack_char).collect();
assert_eq!(bases, vec![0,1,2,3,1,0,3,1,3,1,0,2,0,2,3,2,0,3,2]);
let ascii_bases: Vec<u8> = ascii_seq.as_slice().iter_bp().collect();
assert_eq!(ascii_bases, bases);
let packed_bases: Vec<u8> = ascii_seq.as_slice().iter_bp().collect();
assert_eq!(packed_bases, bases);

// Iterate over 8 chunks at the same time.
let seq = b"AAAACCTTGGTTACTG"; // plain ASCII sequence
// chunks:  ^ ^ ^ ^ ^ ^ ^ ^
let (par_iter, padding) = seq.as_slice().par_iter_bp(1);
let mut par_iter_u8 = par_iter.map(|x| x.as_array_ref().map(|c| c as u8));
assert_eq!(par_iter_u8.next(), Some(*b"AACTGTAT"));
assert_eq!(par_iter_u8.next(), Some(*b"AACTGTCG"));
assert_eq!(par_iter_u8.next(), None);
```
