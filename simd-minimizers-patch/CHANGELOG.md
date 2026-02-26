# Changelog

<!-- next-header -->

## git

## 2.3.0
- feat: add support for open and closed syncmers.
- code: simplify/cleanup dedup intrinsics
- ci: skip some tests in debug mode
- misc: `cargo update` to get `ensure_simd`

## 2.2.0
- Feat: add `run_with_buf` that reuses internal buffers.
- Perf: more `inline(always)`.
- Perf: avoid `u32x8::splat` on NEON.
- Doc: clarify `target-cpu=native` in readme.

## 2.1.0
- **Breaking:**: bump `packed-seq` to `4.1` with `PackedNSeq`.
- Add `scalar` feature that forwards to `packed-seq`, to disable warnings about
  lack of SIMD features when native compilation is not enabled.

## 2.0.1 (yanked)
- Drop accidental `#[feature(type-changing-struct-update)]`

## 2.0.0 (yanked)
- **Breaking**: Migrate to `packed-seq` `4.0` with `PaddedIt`.
- **Breaking**: Improve simd versions for short sequences by reusing allocated buffers.
- **Feature**: Use `seq-hash` crate to cleanly support multiple kmer hashers; add `*_with_hasher` function variants.
- **Feature**: Add support for `_skip_ambiguous_windows` variants that omit minimizers for windows containing non-ACTG (such as NYR).
- Move `simd-minimizers` crate from subdir into the repo root.
- Cleanups (use Rust 2024 edition; drop some public-but-unused functions).
- Improve simd versions for short sequences by preventing over-allocating the output vector.

## 1.4
- Make `NtHash` and `MulHash` seedable.

## 1.3
- Update to `packed-seq` `3.2.1` for `u128` kmer value support.
- Add `iter_{canonical}_minimizer_values_u128` to iterate over `u128` kmer values.

## 1.2
- Fix #10: Add `iter_{canonical}_minimizer_values` to convert positions into `u64` kmer values.
- Update to `packed-seq` `3.0`.
- Fix to properly initialize arrays when collecting super-k-mers.
- Update `packed-seq` to support non-byte offsets.

## 1.1
- Update `packed-seq` to `2.0`, which uses tuples of (simd iterator, padding),
  instead of a separate scalar iterator over the tail.

## 1.0
- Initial release.
