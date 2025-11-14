# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.13.1] - 2025-11-14

### Changed

- `deacon index fetch` now uses `minreq` rather than `ureq` to download indexes, removing dependency on `rustls` & `ring`, which caused a curious build error on Bioconda arm64 MacOS runners, meaning that 0.13.0 was released via Cargo only.

## [0.13.0] - 2025-11-11

### Added

- Command `deacon index fetch` for downloading prebuilt indexes by name. If no index name is specified, the `panhuman-1` index is downloaded.
- Parallel gzip compression of output files with automatic ~1:1 thread allocation between filtering and compression tasks if .gz extensions for `--output` (`-o`) and `--output2` (`-O`) arguments are detected.
  - Automatic thread allocation can be overriden using new `--compression-threads` argument.
  - ~3x faster filtering when reading and writing gzip-compressed Illumina FASTQs.

## [0.12.0] - 2025-10-16

### Added

- Command `deacon index intersect` for finding the intersection of two or more minimizer indexes.
-  Command `deacon index dump` for extracting minimizers from an index as plain text (FASTA).
- Command `deacon cite` showing citation info.

### Changed
- Graceful handling of empty compressed input files.
- Fixes bug where `--debug` incorrectly showed the complement of the hitting minimizer.
- Uses paraseq 0.4.3, addressing a bug identified in paraseq 0.4.2 causing FASTQ records without a trailing newline byte to be ignored.
- 2x increase in filtering throughput on arm64 / MacOS systems enabled by a series of optimisations in latest versions of libraries packed-seq, seq-hash and simd-minimizers.


## [0.11.0] - 2025-10-07

### Added

- Local (socket) server mode, enabling successive filter commands to be handled by a persistent server process for low latency filtering.
- Support for longer k-mers of up to length 61, where k+w ≤ 96 (packed-seq 4.1.1). 

### Changed

- Much faster paired read filtering, particularly from separate input files which are now decompressed in parallel.
- Faster indexing and index loading.
- While minimizers containing non-ACGT nucleotides were already discarded, minimizer selection could still be influenced by non-ACGT nucleotides present in the window, occasionally impacting results. Enabled by changes in simd-minimizers ≥ 2.0, entire windows containing non-ACGT nucleotides are discarded. Records containing non-ACGT nucleotides may therefore be classified differently in this release.
- Redesigned index format (v3).
  - Index now stores 'concrete' k-mers using using 2*k rather than 64bit `xxh3` k-mer hashes.
  - Eliminates [tiny] risk of false positive matches caused by xxh3 collisions.
  - Serialised k-mers are byte-aligned, balancing efficient storage and deserialisation speed.
  - Index disk footprint reduced by 10%.
  - Paves way for painless future adoption of faster HashSet implementations.
  - Paves way for future index introspection functionality.
- `RapidHashSet` (`rapidhash::fast`) replaces combined use of `xxHash` (`xxh3`) and `FxHashSet` 
- Fails gracefully given empty input files.
- Bugfix for paired read I/O.
- Feature gating for reduced compile times.

### Removed

- Removed `--capacity` argument, which was easily misused for little performance benefit.

## [0.10.0] - 2025-09-01

### Added

- Support for k-mer length up to 57 (previously 32).

## [0.9.0] - 2025-08-15

### Changed

- Performance optimisations deliver up to 80% faster filtering with unchanged accuracy.

  - \>2Gbp/s with uncompressed long read input.

  - \>500Mbp/s with gzip-compressed long read input.

## [0.8.1] - 2025-08-14

### Changed

- Fixes bug handling multiline FASTA input introduced in 0.8.0.
- Fixes bug handling paired reads introduced in 0.8.0 which could lead to mispaired read output.

## [0.8.0] - 2025-08-11

### Added

- Added new independent absolute (`-a`) and relative (`-r`) match thresholds with respective default values of 2 and 0.01 (1%). The new default relative threshold improves specificity for long sequences over the previous absolute-only default threshold without affecting short read accuracy. These replace the previous dual purpose `-m` parameter which could accept _either_ an absolute (integer) threshold _or_ a relative (float) threshold. 
- `deacon index` now offers the ability to discard minimizers with information content below a specified scaled Shannon `--entropy` (`-e`) threshold. This is disabled by default.
- `deacon filter` now has a `--debug` mode which prints all records with minimizer matches to stderr including the matched minimizer sequence(s).
- The default worst-case hash table capacity preallocation used in `deacon index union` operations can now be overriden with the new `--capacity` (`-c`) argument, in similar fashion to `deacon index build`. 

### Changed

- Filtering performance has improved dramatically on multicore systems due to improved work allocation using the Paraseq library. Filtering at >1Gbp/s is possible with uncompressed long sequences, and >500Mbp/s is achievable on many systems with Gzip-compressed long reads.
- Minimizers containing ambiguous nucleotides are now ignored.

### Removed

- The filtering argument `--matches` (`-m`) has been removed and replaced with `--abs-threshold` (`-a`) and `--rel-threshold` (`-r`).




## [0.7.0] - 2025-07-08

### Added
- `deacon index diff` optionally accepts a fastx file or stream in place of a second index. This enables index masking using massive sequence collections without the need to first index them.

### Changed
- Deacon uses the recently added `simd-minimizers::iter_canonical_minimizer_values()`, increasing filtering speed by up to 50% on Linux/x86_64 systems. Speeds of 1Gbp/s are now possible with uncompressed FASTA input.
  - Index format is now version 2. Existing indexes must be rebuilt for use with this version. A new version of the panhuman-1 index is available from Zenodo and object storage. Attempting to load an incompatible index throws an error.
- Position-dependent IUPAC ambiguous base canonicalisation was replaced with a simpler and faster fixed mapping, meaning that records containing ambiguous IUPAC bases may be classified differently to before.
- `deacon index union` now automatically preallocates the required hash table capacity, eliminating slowdowns when combining indexes.
- Compatible minimizer _k_ and _w_ is now validated (k+w-1 must be odd) prior to indexing.
- Default index capacity is now 400M (Was 500).


## [0.6.0] - 2025-06-25

### Added

- Support for the .xz compression format via liblzma.
- Adjustable filter output compression level with `--compression-level`.
- Report fields `seqs_out_proportion` and `bp_out_proportion`.

### Changed

- Use zlib-rs for much faster gzip decompression.
- Displays number and proportion of _retained_ reads and base pairs during filtering.

## [0.5.0] - 2025-06-11

### Added

- `--deplete` (`-d`) flag to remove index matches.
- Support for relative thresholds (floats between 0.0 and 1.0) for required minimizer hits to `--matches`.
- `-O` short argument name for `--output2`.
- Tests.

### Changed

- Default filtering behaviour now _passes_ index matches. Use `--deplete` (`-d`) to remove matches.
- Renamed `--nucleotides` (`-n`) to `--prefix-length` (`-p`).
- Renamed `--report` to `--summary`.

### Removed

- `--invert` argument has been removed.

## [0.4.0] - 2025-05-23

### Added

- Non-interleaved paired output file support (`--output2`).

### Changed

- Faster indexing.
- Renamed `--log` argument to `--report`.
- Filter stats are now always sent to stderr (a json report can be written wherever one chooses).
- For paired input sequences, identical minimizer hits in both mates of a read pair are now counted only once.

## [0.3.0] - 2025-05-09

### Added

- Parallel filtering.
  - Up to 10x faster from initial testing.
  - Configurable with new `--threads` parameter.
  - Uses available CPU cores by default (`--threads 0`).
- Tests.

### Changed

- Default minimizer parameters changed to k=31 and w=15.

## [0.2.0] - 2025-04-23

### Added

- Paired read support using either third positional argument or interleaved stdin.

### Changed

- Faster indexing.
- More accurate default parameters.
- Optional argument changes.
- Refactored lib.rs.
- Dependency updates.
  - Bincode2.

## [0.1.0] - 2025-03-14

### Added

- Initial experimental release.
