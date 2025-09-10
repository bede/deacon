![CI status](https://img.shields.io/github/actions/workflow/status/bede/deacon/test.yml?style=flat-square)
[![Crates.io version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)
[![Conda version](https://img.shields.io/conda/v/bioconda/deacon?style=flat-square&label=bioconda&color=blue)](https://anaconda.org/bioconda/deacon)
[![Crates.io downloads](https://img.shields.io/crates/d/deacon?color=orange&label=crates.io%20downloads&style=flat-square)](https://crates.io/crates/deacon)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/deacon.svg?style=flat-square&label=conda%20downloads&color=blue)](https://anaconda.org/bioconda/deacon)
[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

# Deacon

<div align="center"><img src="deacon.png" width="180" alt="Logo"></div>

Search and depletion of FASTA/FASTQ files and streams using accelerated minimizer matching. Default parameters balance sensitivity and specificity for the application of microbial metagenomic host depletion, for which a validated prebuilt index is available. Classification sensitivity, specificity and memory requirements may be tuned by varying *k*-mer length (`-k`), window size (`-w`), and the two match thresholds (`-a` and `-r`). Minimizer `k` and `w` are chosen at index time, while the match thresholds can be chosen at filter time. To be considered a match, sequences must meet both an absolute threshold (`-a`, default 2 minimizer hits) and a relative threshold (`-r`, default 0.01 or 1% of minimizers). Paired sequences are also supported: a match in either mate causes both mates in the pair to be retained or discarded; `deacon filter` retains only matches by default (search mode) and discards matches in `--deplete` mode. Deacon reports filtering performance during execution and optionally writes a JSON `--summary` upon completion. Sequences can optionally be renamed using `--rename` for privacy and smaller file sizes. Gzip, zst and xz compression formats are natively supported and detected by file extension.

Deacon is capable of filtering compressed long reads at >500Mbp/s and indexing a human genome in <30s (Apple M1). Filtering at >2Gbp/s is possible with uncompressed input. Peak memory usage during filtering is 5GB for the default panhuman index. Use Zstandard (zst) compression and/or pipe output to an external compressor such as `pigz` for best performance.

Benchmarks for panhuman host depletion of complex microbial metagenomes are described in a [preprint](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1). Among tested approaches, Deacon with the panhuman-1 (*k*=31, w=15) index exhibited the highest balanced accuracy for both long and short simulated reads. Deacon was however less specific than Hostile for short reads.

> [!IMPORTANT]
> Deacon is actively developed and unstable. Take note of software and index version(s) used in order to guarantee reproducibility of your results. Carefully review the CHANGELOG when updating. Version 0.7.0 introduced a new index container format that is incompatible with prior versions. Please report any problems you encounter by creating an issue or using the email address in my profile.

## Install

### cargo [![Crates.io version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)

```bash
cargo install deacon
```

### conda/mamba/pixi  [![Conda version](https://img.shields.io/conda/v/bioconda/deacon?style=flat-square&label=bioconda&color=blue)](https://anaconda.org/bioconda/deacon)

```bash
conda install -c bioconda deacon
```

## Usage

### Indexing

Use `deacon index build` to quickly build custom indexes. For human host depletion, the prebuilt validated panhuman index is recommended, available for download below from Zenodo or faster object storage. Object storage is provided by the [ModMedMicro research unit](https://www.expmedndm.ox.ac.uk/modernising-medical-microbiology) at the University of Oxford.

```shell
deacon index build chm13v2.fa > human.k31w15.idx

# Discard very low complexity minimizers
deacon index build -e 0.5 chm13v2.fa > human.k31w15e5.idx
```

#### Prebuilt indexes

|                           Name/URL                           |                         Composition                          | Minimizers  | Subtracted minimizers | Size  | Date    |
| :----------------------------------------------------------: | :----------------------------------------------------------: | ----------- | --------------------- | ----- | ------- |
| **panhuman-1 (*k*=31, *w*=15)** [Cloud](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/2/panhuman-1.k31w15.idx), [Zenodo](https://zenodo.org/records/15838532) | [HPRC Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index) ∪ [`CHM13v2.0`](https://www.ncbi.nlm.nih.gov/assembly/11828891) ∪ [`GRCh38.p14`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40) - bacteria (FDA-ARGOS) - viruses (RefSeq) | 409,913,780 | 20,781 (**0.0051%**)  | 3.7GB | 2025-07 |
| **panmouse-1a (k=31, w=15, e=0.5)** [Cloud](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/2/panmouse-1.k31w15e05.idx) | [`GRCm39`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27) ∪ [`PRJEB47108`](https://www.ebi.ac.uk/ena/browser/view/PRJEB47108?show=sequences) - bacteria (FDA-ARGOS) - viruses (RefSeq) | 548,331,948 | 8,246 (**0.0015%**)   | 4.6GB | 2025-08 |

### Filtering

The main command `deacon filter` accepts an index path followed by up to two query FASTA/FASTQ file paths, depending on whether query sequences originate from stdin, a single file, or paired input files. Paired queries are supported as either separate files or interleaved stdin, and written interleaved to either stdout or file, or else to separate paired output files. For paired reads, distinct minimizer hits originating from either mate are counted. By default, query sequences must meet both an absolute threshold of 2 minimizer hits (`-a 2`) and a relative threshold of 1% of minimizers (`-r 0.01`) to pass the filter. Filtering can be inverted for e.g. host depletion using the `--deplete` (`-d`) flag. Gzip, Zstandard, and xz compression formats are detected automatically by file extension. Use Zstandard compression rather than Gzip where possible for best performance.

**Examples**

```bash
# Keep only human sequences
deacon filter panhuman-1.k31w15.idx reads.fq.gz > filt.fq

# Host depletion using the panhuman-1 index and default thresholds
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz

# Max sensitivity with absolute threshold of 1 and no relative threshold
deacon filter -d -a 1 -r 0 panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz

# More specific 10% relative match threshold
deacon filter -d -r 0.1 panhuman-1.k31w15.idx reads.fq.gz > filt.fq.gz

# Stdin and stdout
zcat reads.fq.gz | deacon filter -d panhuman-1.k31w15.idx > filt.fq

# Faster Zstandard compression
deacon filter -d panhuman-1.k31w15.idx reads.fq.zst -o filt.fq.zst

# Fast gzip with pigz
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz | pigz > filt.fq.gz

# Paired reads
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz > filt12.fq
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz -o filt.r1.fq.gz -O filt.r2.fq.gz
zcat r12.fq.gz | deacon filter -d panhuman-1.k31w15.idx - - > filt12.fq

# Save summary JSON
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz -s summary.json

# Replace read headers with incrementing integers
deacon filter -d -R panhuman-1.k31w15.idx reads.fq.gz > filt.fq

# Only look for minimizer hits inside the first 1000bp per record
deacon filter -d -p 1000 panhuman-1.k31w15.idx reads.fq.gz > filt.fq

# Debug mode: see sequences with minimizer hits in stderr
deacon filter -d --debug panhuman-1.k31w15.idx reads.fq.gz > filt.fq
```



## Command line reference

### Filtering

```bash
$ deacon filter -h
Keep or discard DNA fastx records with sufficient minimizer hits to the index

Usage: deacon filter [OPTIONS] <INDEX> [INPUT] [INPUT2]

Arguments:
  <INDEX>   Path to minimizer index file
  [INPUT]   Optional path to fastx file (or - for stdin) [default: -]
  [INPUT2]  Optional path to second paired fastx file (or - for interleaved stdin)

Options:
  -o, --output <OUTPUT>
          Path to output fastx file (or - for stdout; detects .gz and .zst) [default: -]
  -O, --output2 <OUTPUT2>
          Optional path to second paired output fastx file (detects .gz and .zst)
  -a, --abs-threshold <ABS_THRESHOLD>
          Minimum absolute number of minimizer hits for a match [default: 2]
  -r, --rel-threshold <REL_THRESHOLD>
          Minimum relative proportion (0.0-1.0) of minimizer hits for a match [default: 0.01]
  -p, --prefix-length <PREFIX_LENGTH>
          Search only the first N nucleotides per sequence (0 = entire sequence) [default: 0]
  -d, --deplete
          Discard matching sequences (invert filtering behaviour)
  -R, --rename
          Replace sequence headers with incrementing numbers
  -s, --summary <SUMMARY>
          Path to JSON summary output file
  -t, --threads <THREADS>
          Number of execution threads (0 = auto) [default: 8]
      --compression-level <COMPRESSION_LEVEL>
          Output compression level (1-9 for gz & xz; 1-22 for zstd) [default: 2]
      --debug
          Output sequences with minimizer hits to stderr
  -q, --quiet
          Suppress progress reporting
  -h, --help
          Print help
```

### Indexing

```bash
$ deacon index -h
Create and compose minimizer indexes

Usage: deacon index <COMMAND>

Commands:
  build  Index minimizers contained within a fastx file
  info   Show index information
  union  Combine multiple minimizer indexes (A ∪ B…)
  diff   Subtract minimizers in one index from another (A - B)
  help   Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

```bash
$ deacon index build -h
Index minimizers contained within a fastx file

Usage: deacon index build [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Path to input fastx file (supports gz, zst and xz compression)

Options:
  -k <KMER_LENGTH>
          K-mer length used for indexing (1-32) [default: 31]
  -w <WINDOW_SIZE>
          Minimizer window size used for indexing [default: 15]
  -o, --output <OUTPUT>
          Path to output file (- for stdout) [default: -]
  -c, --capacity <CAPACITY_MILLIONS>
          Preallocated index capacity in millions of minimizers [default: 400]
  -t, --threads <THREADS>
          Number of execution threads (0 = auto) [default: 8]
  -q, --quiet
          Suppress sequence header output
  -e, --entropy-threshold <ENTROPY_THRESHOLD>
          Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
  -h, --help
          Print help
```

## Building custom indexes

Building custom Deacon indexes is quite fast. Nevertheless, when indexing many large genomes, it may be worthwhile separately indexing and subsequently combining indexes into one succinct index. Combine distinct minimizers from multiple indexes using `deacon index union`. Similarly, use `deacon index diff` to subtract the minimizers contained in one index from another. This can be helpful  for e.g. eliminating shared minimizers between the target and host genomes when building custom (non-human) indexes for host depletion.

- Use `deacon index union 1.idx 2.idx 3.idx… > 1+2+3.idx` to succinctly combine two (or more!) deacon indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 1.idx from 2.idx. Useful for masking out shared minimizer content between e.g. target and host genomes.
- In version `0.7.0` and above, `deacon index diff` also supports subtracting minimizers from an index using a fastx file or stream, e.g. `deacon index diff 1.idx 2.fa.gz > 1-2.idx` or `zcat *.fa.gz | deacon index diff 1.idx - > 1-2.idx`.

For best performance, set the `--capacity` argument of `deacon index build` to a number of minimizers in millions greater than that you expect your index to contain. Setting this too low will cause delays during indexing for hash table resizing.

## Filtering summary statistics

Use `-s summary.json` to save detailed filtering statistics:
```json
{
  "version": "deacon 0.9.0",
  "index": "panhuman-1.k31w15.idx",
  "input": "HG02334.1m.fastq.gz",
  "input2": null,
  "output": "-",
  "output2": null,
  "k": 31,
  "w": 15,
  "abs_threshold": 2,
  "rel_threshold": 0.01,
  "prefix_length": 0,
  "deplete": true,
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

## Features
There is an optional feature which can be enabled when building with `cargo build --features server`.

This enables the running of a server which pre-loads the index, allowing filtering to be remote rather than local. In most local use cases, this will inevitably be slower than local `filter`, but for cases where lots of small inputs would otherwise load the index for each input, this is expected to be faster.

> [!NOTE]
> Compiling with the server feature swaps the filtering engine for a slower one in order to enable faster client performance. For optimal local filtering, **do not** use this feature. Scientifically the results are identical between both modes.

### Testing
Unit tests covering all filter functionality (but with the server/client) are enabled when compiling with `--features server`. However, due to concurrency issues, the client tests can only be run single threaded, and will fail without this:
```bash
cargo test --features server -- --test-threads 1
```

This also tests the remote filtering engine in local mode to ensure that results are identical between the two filtering implementations.

### Server
Start up a server with a specific index loaded. Note that this by default runs in the foreground, so it may appear as if nothing is happening after a message about loading your index. To display incomming connection logs, set `RUST_LOG=trace` in your environment variables. 


#### Run on default port
Starts up the server on port `8888`
```bash
deacon server index.idx
```
#### Run on custom port
Starts up the server on port `12345`
```bash
deacon server index.idx -p 12345
```

### Client
Almost exactly identical to the `deacon filter` reference, but swapping index path for server address. 


#### Usage
```bash
Alternate version of Filter, swapping local compute for passing to a server which has the index pre-loaded. Will inevitably be slower than local filtering, but saves on index loading. Better used for cases of small input + large index

Requires "server" feature to be enabled at compile time.

Usage: deacon client [OPTIONS] <SERVER_ADDRESS> [INPUT] [INPUT2]

Arguments:
  <SERVER_ADDRESS>  Server address to connect to (including port)
  [INPUT]           Optional path to fastx file (or - for stdin) [default: -]
  [INPUT2]          Optional path to second paired fastx file (or - for interleaved stdin)

Options:
  -o, --output <OUTPUT>
          Path to output fastx file (or - for stdout; detects .gz and .zst) [default: -]
  -O, --output2 <OUTPUT2>
          Optional path to second paired output fastx file (detects .gz and .zst)
  -a, --abs-threshold <ABS_THRESHOLD>
          Minimum absolute number of minimizer hits for a match [default: 2]
  -r, --rel-threshold <REL_THRESHOLD>
          Minimum relative proportion (0.0-1.0) of minimizer hits for a match [default: 0.01]
  -p, --prefix-length <PREFIX_LENGTH>
          Search only the first N nucleotides per sequence (0 = entire sequence) [default: 0]
  -d, --deplete
          Discard matching sequences (invert filtering behaviour)
  -R, --rename
          Replace sequence headers with incrementing numbers
  -s, --summary <SUMMARY>
          Path to JSON summary output file
  -t, --threads <THREADS>
          Number of execution threads (0 = auto) [default: 8]
      --compression-level <COMPRESSION_LEVEL>
          Output compression level (1-9 for gz & xz; 1-22 for zstd) [default: 2]
      --debug
          Output sequences with minimizer hits to stderr
  -q, --quiet
          Suppress progress reporting
  -h, --help
          Print help
```

**Examples**

```bash
# Keep only sequences matching the index loaded in the server
deacon client http://0.0.0.0:8888 reads.fq.gz > filt.fq

# Host depletion using the index loaded in the server and default thresholds
deacon client -d http://0.0.0.0:8888 reads.fq.gz -o filt.fq.gz

# Max sensitivity with absolute threshold of 1 and no relative threshold
deacon client -d -a 1 -r 0 http://0.0.0.0:8888 reads.fq.gz -o filt.fq.gz

# More specific 10% relative match threshold
deacon client -d -r 0.1 http://0.0.0.0:8888 reads.fq.gz > filt.fq.gz

# Stdin and stdout
zcat reads.fq.gz | deacon client -d http://0.0.0.0:8888 > filt.fq

# Faster Zstandard compression
deacon client -d http://0.0.0.0:8888 reads.fq.zst -o filt.fq.zst

# Fast gzip with pigz
deacon client -d http://0.0.0.0:8888 reads.fq.gz | pigz > filt.fq.gz

# Paired reads
deacon client -d http://0.0.0.0:8888 r1.fq.gz r2.fq.gz > filt12.fq
deacon client -d http://0.0.0.0:8888 r1.fq.gz r2.fq.gz -o filt.r1.fq.gz -O filt.r2.fq.gz
zcat r12.fq.gz | deacon client -d http://0.0.0.0:8888 - - > filt12.fq

# Save summary JSON
deacon client -d http://0.0.0.0:8888 reads.fq.gz -o filt.fq.gz -s summary.json

# Replace read headers with incrementing integers
deacon client -d -R http://0.0.0.0:8888 reads.fq.gz > filt.fq

# Only look for minimizer hits inside the first 1000bp per record
deacon client -d -p 1000 http://0.0.0.0:8888 reads.fq.gz > filt.fq

# Debug mode: see sequences with minimizer hits in stderr
deacon client -d --debug http://0.0.0.0:8888 reads.fq.gz > filt.fq
```

## Citation

[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

>  Bede Constantinides, John Lees, Derrick W Crook. "Deacon: fast sequence filtering and contaminant depletion" *bioRxiv* 2025.06.09.658732, https://doi.org/10.1101/2025.06.09.658732 

Please also consider citing the SimdMinimizers paper:

> Ragnar Groot Koerkamp, Igor Martayan. "SimdMinimizers: Computing random minimizers, *fast*" *bioRxiv* 2025.01.27.634998, https://doi.org/10.1101/2025.01.27.634998 

