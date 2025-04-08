[![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon) [![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon)

# Deacon

A minimizer-based filter for nucleotide sequences in FASTA or FASTQ format, built for efficient host depletion. Default behaviour removes query sequences with two or more minimizers present in the index. Filters at ~50Mbp/s using a single Apple M1 core and indexes the human genome in under 60s. Peak memory usage is ~2.5GB for a human genome with default parameters. Accuracy benchmarks will be published soon.

The sensitivity/specificity/memory tradeoff can be tuned using *k*-mer length (`-k`), minimizer window length (`-w`), and match threshold (`-m`). Filtering speed may be increased by considering only the first `-n` bases per query sequence. Uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers) for accelerated minimizer computation. This project is currently unstable and under active development.

## Install

### conda/mamba/pixi  [![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon)

```bash
conda install -c bioconda deacon
```

### cargo [![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)

```bash
cargo install deacon  # Requires Rust toolchain
```

## Usage

### Indexing

Supports FASTA[.gz] input files and outputs to stdout or file (`-o`).

```bash
deacon index build chm13v2.fa > human.k31w21.idx
deacon index build -k 41 -m 27 chm13v2.fa > human.k41w27.idx  # Custom minimizer k and w
```

### Filtering

Supports FASTA or FASTQ input from stdin or file and outputs to stdout or file. Paired sequences are supported as either separate files or interleaved stdin, and are  written in interleaved format to either stdout or file. Gzip (.gz) and Zstandard (.zst) compression formats are detected automatically. Piping uncompressed FASTA/Q to pigz is advisable in order to avoid compression bottlenecks when writing gzip output directly.

```bash
deacon filter human.idx reads.fq.gz -o filt.fastq  # File input & output
zcat reads.fq.gz | deacon filter human.idx | pigz > filt.fq.gz  # Fast gzip
zcat reads.fq.gz | deacon -n 1000 filter human.idx | pigz > filt.fq.gz  # Faster
deacon filter -m 3 human.idx reads.fq.gz | pigz > filt.fq.gz  # More precise
deacon filter -m 1 human.idx reads.fq.gz | pigz > filt.fq.gz  # More sensitive
deacon filter human.idx r1.fq.gz r2.fq.gz > filt12.fastq  # Paired file input
zcat r12.fq.gz | deacon filter human.idx - - > filt12.fastq  # Interleaved stdin
deacon filter human.idx reads.fq.gz --log log.json > filt.fq  # Log results JSON
```

## Reports

Use `--log results.json` to save a filtering summary to a JSON file:
```json
{
  "version": "0.2.0",
  "index": "data/chm13v2.k31w21.idx",
  "input1": "data/HG02334.1m.fastq.gz",
  "input2": null,
  "output": "-",
  "k": 31,
  "w": 21,
  "m": 1,
  "n": 0,
  "invert": false,
  "rename": false,
  "seqs_in": 1000000,
  "seqs_out": 13452,
  "seqs_removed": 986548,
  "seqs_removed_proportion": 0.986548,
  "bp_in": 5477122928,
  "bp_out": 5710050,
  "bp_removed": 5471412878,
  "bp_removed_proportion": 0.9989574727324798,
  "time": 125.755103875,
  "seqs_per_second": 7951,
  "bp_per_second": 43553881
}
```

### Set operations on indexes

- Use `deacon index union 1.idx 2.idx > 1+2.idx` to nonredundantly combine two (or more) deacon minimizer indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 2.idx from 1.idx. Useful for masking.
