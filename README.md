[![Crates.io version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)
[![Conda version](https://img.shields.io/conda/v/bioconda/deacon?style=flat-square&label=bioconda&color=blue)](https://anaconda.org/bioconda/deacon)
[![Crates.io downloads](https://img.shields.io/crates/d/deacon?color=orange&label=crates.io%20downloads&style=flat-square)](https://crates.io/crates/deacon)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/deacon.svg?style=flat-square&label=conda%20downloads&color=blue)](https://anaconda.org/bioconda/deacon)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/tools/list?search=deacon)
[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

<div align="center"><img src="deacon.png" width="180" alt="Logo"></div>

# Deacon

Deacon filters DNA sequences in FASTA/Q files and streams using SIMD-accelerated minimizer comparison with query sequence(s), emitting either matching sequences (**search mode**), or sequences without matches (**deplete mode**). Sequences match when they share enough distinct minimizers with the indexed query to exceed chosen absolute and relative thresholds. Query size has little impact on filtering speed, enabling ultrafast search and depletion with gene-, genome- and pangenome-scale queries using a laptop. Deacon filters uncompressed FASTA/Q at **gigabases per second** on recent AMD, Intel (`x86_64`), and Apple `arm64` systems. Built with panhuman host depletion in mind—yet broadly useful for searching large sequence collections—Deacon delivers [leading classification accuracy](https://doi.org/10.1101/2025.06.09.658732) for host depletion and unrivalled speed using 5GB of RAM.

Default parameters are carefully chosen but easily changed. Classification sensitivity, specificity and memory requirements may be tuned by varying *k*-mer length (`-k`), window size (`-w`), absolute match threshold (`-a`) and relative match threshold (`-r`) . Minimizer `k` and `w` are chosen at query index time, while the match thresholds can be chosen at filter time. Matching sequences are those that share enough distinct minimizers with the indexed query to exceed *both* the absolute threshold (`-a`, default 2 shared minimizers) and the relative threshold (`-r`, default 0.01 [1%] shared minimizers). For paired sequences, hits in either mate counts towards a single match threshold for the pair. Deacon reports filtering performance during execution and optionally writes a JSON `--summary` upon completion. Sequences can optionally be renamed using `--rename` for privacy and smaller file sizes. Deacon fully supports stdin/out (uncompressed fastx) and natively handles .gz, .zst and .xz fastx file IO, as well as BINSEQ CBQ (.cbq & .cba).

Benchmarks for panhuman host depletion of complex microbial metagenomes are described in a [preprint](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1). Deacon with the `panhuman-1` (*k*=31, w=15) index exhibited the highest balanced accuracy for both long and short simulated reads. Deacon was less specific only than Hostile for short reads.

## Use cases

- Depletion of human or other host genome sequences in FASTQ reads or streams.
- Ultrafast binary classification of genes, genomes or pangenomes in terabase genome catalogues like [AllTheBacteria](https://allthebacteria.org/) without tedious pre-indexing.

## Install

### Conda/mamba/pixi  [![Conda version](https://img.shields.io/conda/v/bioconda/deacon?style=flat-square&label=bioconda&color=blue)](https://anaconda.org/bioconda/deacon)

```bash
conda install -c bioconda deacon
```

### Cargo [![Crates.io version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)

```bash
RUSTFLAGS="-C target-cpu=native" cargo install deacon
```

> [!IMPORTANT]
> Cargo installation requires [Rust 1.88 or newer](https://rust-lang.org/tools/install/). Update using `rustup update`.

### Docker [![Crates.io version](https://img.shields.io/badge/install%20with-docker-important.svg?style=flat-square&logo=docker)](https://biocontainers.pro/tools/deacon)

Containers are available from the [BioContainers registry](https://biocontainers.pro/tools/deacon).

```bash
docker pull quay.io/biocontainers/deacon:0.15.0--hdd79491_0
```

## Quickstart

### Ultrafast panhuman host depletion

```bash
# Download validated 3GB human pangenome index (version 0.13.0 or later)
deacon index fetch panhuman-1

# Deplete long reads
deacon filter -d panhuman-1.k31w15.idx reads.fq -o filt.fq

# Deplete short paired reads
deacon filter -d panhuman-1.k31w15.idx reads.r1.fq.gz reads.r2.fq.gz -o filt.r1.fq.gz -O filt.r2.fq.gz
```

### Ultrafast gene/genome/pangenome search

```bash
deacon index build amr-genes.fa > amr-genes.idx
deacon filter amr-genes.idx AllTheBacteria.fa.zst > hits.fa
```

*N.B. Indexing a 3Gbp human genome takes ~30s using 18GB of RAM with default parameters. Filtering uses 5GB.*

## Prebuilt indexes

Prebuilt pangenome indexes are provided. These can be downloaded using the links below, or with `deacon index fetch <name>`.  These may alternatively be reproduced using workflows in [deacon-indexes](https://github.com/bede/deacon-indexes). 

|                          Name & URL                          |                         Composition                          | Minimizers  | Subtracted minimizers | Size  | Date    |
| :----------------------------------------------------------: | :----------------------------------------------------------: | ----------- | --------------------- | ----- | ------- |
| **`panhuman-1` (*k*=31, *w*=15)** [Cloud](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w15.idx), [Zenodo](https://zenodo.org/records/17288185) | [HPRC Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index) ∪ [`CHM13v2.0`](https://www.ncbi.nlm.nih.gov/assembly/11828891) ∪ [`GRCh38.p14`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40) - bacteria (FDA-ARGOS) - viruses (RefSeq) | 409,907,949 | 20,671 (**0.0050%**)  | 3.3GB | 2025-04 |
| **`panmouse-1` (*k*=31, *w*=15)** [Cloud](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panmouse-1.k31w15.idx), [Zenodo](https://zenodo.org/records/17699167) | [`GRCm39`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27) ∪ [`PRJEB47108`](https://www.ebi.ac.uk/ena/browser/view/PRJEB47108?show=sequences) - bacteria (FDA-ARGOS) - viruses (RefSeq) | 551,041,865 | 9,866 (**0.0018%**)   | 4.4GB | 2025-11 |

> [!NOTE]
>
> **Index compatibility.** Deacon `0.11.0` and above uses index format version 3. Using version 3 indexes with older Deacon versions and vice versa triggers an error. Prebuilt indexes in legacy formats are archived in object storage and Zenodo to ensure  reproducibility. To download indexes in legacy formats, replace the `/3/` in any prebuilt index download URL with either `/2/` or `/1/`  accordingly.
>
> - Deacon **`0.11.0`** and above uses index format version **`3`**
> - Deacon **`0.7.0`** through to **`0.10.0`** used index format version **`2`**
> - Deacon **`0.1.0`** through to **`0.6.0`** used index format version **`1`**

## Usage

### Filtering

The main command `deacon filter` accepts an index path followed by up to two sequence file paths, depending on whether input sequences originate from stdin, a single file, or paired input files. Indexes are built with `deacon index build`.  Paired inputs are supported as either two separate or one interleaved file/stream when using `--interleaved`, and may be written either to separate paired output files or one interleaved file. For paired sequences, distinct minimizer hits originating from either mate are counted. Paired read headers can be validated using `--check-pairs`. By default, input sequences must meet both an absolute threshold of 2 minimizer hits (`-a 2`) and a relative threshold of 1% of minimizers (`-r 0.01`) to pass the filter. Filtering can be inverted for e.g. host depletion using the `--deplete` (`-d`) flag. Use `--inverse-output` to write discarded records at the same time; primary and inverse outputs must both be FASTX or both be CBQ/CBA, though their compression and paired-file layouts may differ. Gzip, Zstandard, and xz compressed FASTX formats are detected automatically by file extension. Single and paired [BINSEQ](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v2) CBQ files are natively supported for both input and output. CBQ input is detected by content rather than file extension, while a `.cbq` output path writes CBQ, and a `.cba` output path writes CBQ without quality.

#### Examples

```bash
# Keep only sequences matching a collection of genes
deacon index build genes.fa > genes.idx
deacon filter genes.idx sequences.fa.gz -o matches.fa.gz

# Host depletion using the panhuman-1 index
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz

# Write retained and discarded reads
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz -i host.fq.gz

# High sensitivity host depletion with absolute threshold of 1 and no relative threshold
deacon filter -d -a 1 -r 0 panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz

# High specificity 10% relative match threshold
deacon filter -d -r 0.1 panhuman-1.k31w15.idx reads.fq.gz > filt.fq.gz

# Stdin and stdout
zcat reads.fq.gz | deacon filter -d panhuman-1.k31w15.idx > filt.fq

# True multithreaded gzip decompression with rapidgzip
rapidgzip -dc reads.fq.gz | deacon filter -d panhuman-1.k31w15.idx > filt.fq

# Zstandard compression
deacon filter -d panhuman-1.k31w15.idx reads.fq.zst -o filt.fq.zst

# Paired reads
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz > filt12.fq
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz -o filt.r1.fq.gz -O filt.r2.fq.gz
deacon filter panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz \
  -o matches.r1.fq.gz -O matches.r2.fq.gz -i rest.r1.fq.gz -I rest.r2.fq.gz

# Validate matching Illumina record names while filtering paired sequences
deacon filter -d --check-pairs panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz -o filt.r1.fq.gz -O filt.r2.fq.gz

# Interleaved paired reads (file or stdin)
deacon filter -d --interleaved panhuman-1.k31w15.idx r12.fq.gz > filt12.fq
zcat r12.fq.gz | deacon filter -d panhuman-1.k31w15.idx - - > filt12.fq

# Save summary JSON
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz -s summary.json

# BINSEQ CBQ input and output; paired records are stored in one file
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.cbq
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz -o filt12.cbq
deacon filter -d panhuman-1.k31w15.idx filt12.cbq -o filt12.fq.gz

# A .cba extension writes cbq without quality values, if present
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.cba

# Replace read headers with incrementing integers
deacon filter -d -R panhuman-1.k31w15.idx reads.fq.gz > filt.fq

# Only look for minimizer hits inside the first 1000bp per record
deacon filter -d -p 1000 panhuman-1.k31w15.idx reads.fq.gz > filt.fq

# Output FASTA regardless of input format (discards quality scores)
deacon filter -d -f panhuman-1.k31w15.idx reads.fq.gz > filt.fa

# Debug mode: see sequences with minimizer hits in stderr
deacon filter -d --debug panhuman-1.k31w15.idx reads.fq.gz > filt.fq
```

> [!NOTE]
>
> `deacon filter` uses 8 threads by default. Using more threads (e.g.  `--threads 16`) can accelerate filtering given sufficient resources, especially with uncompressed sequences whose processing is not rate limited by decompression. Since version `0.13.0`, Deacon writes gzipped output files (e.g `-o out.fastq.gz`) in parallel, providing particular practical benefit for gzipped paired reads. If output file(s) ending in `.gz` are detected, total `--threads` are allocated 1:1 to compression and filtering tasks respectively. Gzip compression thread allocation can be overriden with `--compression-threads`.

### Indexing

```bash
# Index one FASTA/FASTQ file
deacon index build genome.fa.gz > genome.idx

# Index many FASTA/FASTQ files using stdin
zcat *.fa.gz | deacon index build - > genomes.idx
```

`deacon index build` accepts either a FASTA/FASTQ file or a stdin stream (`-`), enabling convenient indexing of compressed sequences in one or many files with a single step. Indexing a human genome takes a few seconds. Indexing uses 2-4x as much RAM as filtering. For indexing large collections approaching terabase scale—such as mammalian pangenomes—it may be practical to index genomes individually in parallel and later combine them using the `deacon index union` set operation, described below.

#### Set operations

A differentiating feature of Deacon is the ease of combining, subtracting and intersecting minimizer indexes. For example, `deacon index diff`can be used to subtract shared minimizers between target and host genomes when building custom indexes for host depletion.

- Use `deacon index union 1.idx 2.idx 3.idx… > 1+2+3.idx` to succinctly combine two or more indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 2.idx from 1.idx. Useful for masking out shared minimizer content between e.g. target and host genomes.
  - `deacon index diff` also supports subtracting minimizers from an index using a fastx file or stream directly, e.g. `deacon index diff 1.idx 2.fa.gz > 1-2.idx` or `zcat *.fa.gz | deacon index diff 1.idx - > 1-2.idx`. This enables diffing with larger-than-memory sequence collections if desired.
  - Specifying `w=1` subtracts every single k-mer in the second index/FASTX from the first.

- Use `deacon index intersect 1.idx 2.idx… > 1∩2.idx` to find the intersection of minimizers in two or more indexes.
  - Indexes with `w=1` are accepted in any position after the first, whatever the first index's `w`, reporting which minimizers of the first occur as exact k-mers in them.

#### Inspecting indexes

- Use `deacon index info 1.idx` to display index information including minimizer *k* and *w* parameters, number of minimizers, and index format version.
- Use `deacon index dump 1.idx > 1.fa` to dump a minimizer index to FASTA.

## Command line reference

### Filtering

```bash
$ deacon filter -h
Retain or deplete sequence records with sufficient minimizer hits to the index

Usage: deacon filter [OPTIONS] <INDEX> [INPUT] [INPUT2]

Arguments:
  <INDEX>   Path to minimizer index file
  [INPUT]   Optional path to fastx or binseq cbq file (or - for stdin) [default: -]
  [INPUT2]  Optional path to second paired fastx file

Options:
  -a, --abs-threshold <ABS_THRESHOLD>
          Minimum absolute number of minimizer hits for a match [default: 2]
  -r, --rel-threshold <REL_THRESHOLD>
          Minimum relative proportion (0.0-1.0) of minimizer hits for a match [default: 0.01]
  -p, --prefix-length <PREFIX_LENGTH>
          Search only the first N nucleotides per sequence (0 = entire sequence) [default: 0]
  -c, --complexity-threshold <COMPLEXITY_THRESHOLD>
          Ignore minimizers below this kdust complexity threshold (0.0-1.0)
  -d, --deplete
          Discard matching sequences (invert filtering behaviour)
  -R, --rename
          Replace sequence headers with incrementing numbers (deterministic with --ordered)
  -o, --output <OUTPUT>
          Path to output file (fastx to stdout by default; detects .gz, .zst, .xz, .cbq, .cba)
  -O, --output2 <OUTPUT2>
          Optional path to second paired output fastx file (detects .gz, .zst, .xz)
  -i, --inverse-output <INVERSE_OUTPUT>
          Path to inverse output file; container format must match --output
  -I, --inverse-output2 <INVERSE_OUTPUT2>
          Optional path to second paired inverse output fastx file (detects .gz, .zst, .xz)
  -s, --summary <SUMMARY>
          Path to JSON summary output file
  -t, --threads <THREADS>
          Number of threads (0 = auto) [default: 8]
      --compression-threads <COMPRESSION_THREADS>
          Number of threads used for output compression (0 = auto) [default: 0]
      --compression-level <COMPRESSION_LEVEL>
          Output compression level (1-9 for gz & xz; 1-22 for zstd including cbq) [default: 2]
      --cbq-block-size <CBQ_BLOCK_SIZE>
          cbq output block size in MiB (or cbq input block size if higher) [default: 16]
      --discard-quality
          Emit fasta or quality-free cbq regardless of input format
      --interleaved
          Treat INPUT as interleaved paired records from single file or stdin
      --ordered
          Preserve input record ordering (deterministic, slightly slower)
      --check-pairs
          Validate paired record names (Illumina CASAVA or /1 /2 suffixes)
      --debug
          Emit sequences with minimizer hits to stderr
  -q, --quiet
          Suppress progress reporting
  -h, --help
          Print help
```

### Indexing

```bash
$ deacon index -h
Build, inspect, compose and fetch minimizer indexes

Usage: deacon index <COMMAND>

Commands:
  build      Index minimizers contained within a fastx file
  union      Combine multiple minimizer indexes (A ∪ B…)
  intersect  Intersect multiple minimizer indexes (A ∩ B…)
  diff       Subtract minimizers in one index from another (A - B)
  dump       Dump minimizer index to fasta
  filter     Discard minimizers below a complexity threshold
  info       Show index information
  freeze     Freeze an index into a binary fuse filter (BFF) index (k<=32)
  fetch      Fetch a pre-built index from remote storage
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help

```

```bash
$ deacon index build -h
Index minimizers contained within a fastx file

Usage: deacon index build [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Path to input fastx file (or - for stdin; supports gz, zst and xz compression)

Options:
  -k <KMER_LENGTH>         K-mer length used for indexing (k+w-1 must be <= 96 and odd) [default: 31]
  -w <WINDOW_SIZE>         Minimizer window size used for indexing [default: 15]
  -o, --output <OUTPUT>    Path to output file (stdout if not specified)
  -t, --threads <THREADS>  Number of execution threads (0 = auto) [default: 8]
  -q, --quiet              Suppress progress reporting
  -h, --help               Print help
```


## Filtering summary statistics

Use `-s summary.json` to save detailed filtering statistics:
```json
{
  "version": "deacon 0.17.0",
  "index": "panhuman-1.k31w15.idx",
  "input": "HG02334.100MB.fastq.gz",
  "input2": null,
  "output": "-",
  "output2": null,
  "inverse_output": null,
  "inverse_output2": null,
  "k": 31,
  "w": 15,
  "abs_threshold": 2,
  "rel_threshold": 0.01,
  "prefix_length": 0,
  "deplete": true,
  "rename": false,
  "ordered": false,
  "check_pairs": false,
  "seqs_in": 37500,
  "seqs_out": 454,
  "seqs_out_proportion": 0.012106666666666667,
  "seqs_removed": 37046,
  "seqs_removed_proportion": 0.9878933333333333,
  "bp_in": 141474280,
  "bp_out": 227079,
  "bp_out_proportion": 0.001605090338682056,
  "bp_removed": 141247201,
  "bp_removed_proportion": 0.9983949096613179,
  "time": 0.446945667,
  "seqs_per_second": 84129,
  "bp_per_second": 317392822,
  "seqs_per_second_total": 83902,
  "bp_per_second_total": 316535745
}
```

## Server mode

From version `0.11.0`, it is possible to eliminate index loading overhead at the start of each filter operation by preloading the index in the memory of a local server process. Subsequent filtering commands with `--use-server` are executed by the server process using a UNIX socket. Having started a server process, the index of the first filtering command it receives persists in memory for the life of that server process, enabling subsequent filter commands to be served rapidly without HashSet construction overhead.

```bash
# Start the server
deacon server start

# The first filter command loads the index as usual
deacon --use-server filter ref.idx reads.fq > /dev/null

# Subsequent filter commands use the existing index stored in memory
deacon --use-server filter ref.idx reads.fq -o filt.fq -s summary.json

# Stop the server
deacon --use-server server stop
```

## Python bindings

Since `0.16.0`, Deacon can also be used as a [Python package distributed via PyPI](https://pypi.org/project/deacon/). Like server mode, this enables multihreaded filtering with index reuse. Refer to the [Python readme](https://github.com/bede/deacon/blob/main/deacon-py/python/README.md) for more information.

## Workflow manager integration

### Nextflow (nf-core)

- Modules [`deacon_index`](https://nf-co.re/modules/deacon_index) and [`deacon_filter`](https://nf-co.re/modules/deacon_filter)
- Subworkflow [`fastq_index_filter_deacon`](https://nf-co.re/subworkflows/fastq_index_filter_deacon/)

### Galaxy

- Tool Shed suite [`suite_deacon`](https://toolshed.g2.bx.psu.edu/view/iuc/suite_deacon), installable into any Galaxy instance

## Limitations

Deacon is benchmarked using both long reads and 2x150bp short reads. Classification accuracy deteriorates with very short reads however – if classifying old 50bp Illumina reads for example, consider using `--abs-threshold 1` (`-a 1`). 



## Citation

[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

>  Bede Constantinides, John Lees, Derrick W Crook. "Deacon: fast sequence filtering and contaminant depletion" *bioRxiv* 2025.06.09.658732, https://doi.org/10.1101/2025.06.09.658732

Please also consider citing the SimdMinimizers paper:

> Ragnar Groot Koerkamp, Igor Martayan. "SimdMinimizers: Computing random minimizers, *fast*" *bioRxiv* 2025.01.27.634998, https://doi.org/10.1101/2025.01.27.634998
