# simd-minimizers

[![crates.io](https://img.shields.io/crates/v/simd-minimizers)](https://crates.io/crates/simd-minimizers)
[![docs](https://img.shields.io/docsrs/simd-minimizers)](https://docs.rs/simd-minimizers)

A SIMD-accelerated library to compute random minimizers.

It can compute all the minimizers of a human genome in 4 seconds using a single thread.
It also provides a *canonical* version that ensures that a sequence and its reverse-complement always select the same positions, which takes 6 seconds on a human genome.

This crate builds on [`packed_seq`](https://github.com/rust-seq/packed-seq) and
[`seq-hash`](https://github.com/rust-seq/seq-hash).
 
The underlying algorithm is described in the following
[**paper**](https://doi.org/10.4230/LIPIcs.SEA.2025.20): 

- SimdMinimizers: Computing random minimizers, fast.
  Ragnar Groot Koerkamp, Igor Martayan
  SEA 2025 [doi.org/10.4230/LIPIcs.SEA.2025.20](https://doi.org/10.4230/LIPIcs.SEA.2025.20)

## Requirements

This library requires AVX2 or NEON instruction sets, which, on x64, requires
either `target-cpu=native` or `target-cpu=x86-64-v3`.
See [this README](https://github.com/ragnargrootkoerkamp/ensure_simd) for details and [this
blog](https://curiouscoding.nl/posts/distributing-rust-simd-binaries/) for background.
The same restrictions apply when using simd-minimizers in a larger project.

``` sh
RUSTFLAGS="-C target-cpu=native" cargo run --release
```

## Usage example

Full documentation can be found on [docs.rs](https://docs.rs/simd-minimizers).

```rust
use packed_seq::{PackedSeqVec, SeqVec};

let seq = b"ACGTGCTCAGAGACTCAGAGGA";
let packed_seq = PackedSeqVec::from_ascii(seq);

let k = 5;
let w = 7;
let hasher = <seq_hash::NtHasher>::new(k);

// Simple usage with default hasher, returning only positions.
let minimizer_positions = canonical_minimizer_positions(packed_seq.as_slice(), k, w);
assert_eq!(minimizer_positions, vec![0, 7, 9, 15]);

// Advanced usage with custom hasher, super-kmer positions, and minimizer values as well.
let mut minimizer_positions = Vec::new();
let mut super_kmers = Vec::new();
let minimizer_vals: Vec<u64> = canonical_minimizers(k, w)
    .hasher(&hasher)
    .super_kmers(&mut super_kmers)
    .run(packed_seq.as_slice(), &mut minimizer_positions)
    .values_u64()
    .collect();

// Compute _syncmers_ positions and values instead:
let mut syncmer_positions = Vec::new();
// List of (k+w-1)-mer values.
let syncmer_vals: Vec<u64> = canonical_syncmers(k, w)
    .run(packed_seq.as_slice(), &mut syncmer_positions)
    .values_u64()
    .collect();
```

## Benchmarks

Benchmarks can be found in the `bench` directory in the GitHub repository.

`bench/benches/bench.rs` contains benchmarks used in [this blogpost](https://curiouscoding.nl/posts/fast-minimizers/).

`bench/src/bin/paper.rs` contains benchmarks used in the paper.

Note that the benchmarks require some nightly features, you can install the latest nightly version with

```sh
rustup install nightly
```

To replicate results from the paper, go into `bench` and run
```sh
RUSTFLAGS="-C target-cpu=native" cargo +nightly run --release
python eval.py
```

The human genome we use is from the T2T consortium, and available by following
the first link [here](https://github.com/marbl/CHM13?tab=readme-ov-file#t2t-chm13v20-t2t-chm13y).
