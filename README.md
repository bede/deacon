[![Crates.io version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)
[![Conda version](https://img.shields.io/conda/v/bioconda/deacon?style=flat-square&label=bioconda&color=blue)](https://anaconda.org/bioconda/deacon)
[![Crates.io downloads](https://img.shields.io/crates/d/deacon?color=orange&label=crates.io%20downloads&style=flat-square)](https://crates.io/crates/deacon)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/deacon.svg?style=flat-square&label=conda%20downloads&color=blue)](https://anaconda.org/bioconda/deacon)
[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

# Deacon

<div align="center"><img src="deacon.png" width="180" alt="Logo"></div>

Deacon filters DNA sequences in FASTA/Q files and streams using accelerated minimizer comparison with indexed query sequence(s), emitting either matching sequences (**search mode**), or sequences without matches (**deplete mode**). Sequences are matched when they share enough distinct minimizers with the indexed query to exceed specified absolute and relative thresholds. Query size has little impact on filtering speed, enabling ultrafast search and depletion with gene-, genome- and pangenome-scale queries using a laptop. Deacon filters uncompressed FASTA/Q at gigabases per second on recent AMD, Intel (`x86_64`), and Apple `arm64` systems. Deacon is built with panhuman host depletion in mind, delivering [leading classification accuracy](https://doi.org/10.1101/2025.06.09.658732) for this application with unrivalled speed and a 5GB memory footprint.

Default parameters are carefully chosen and easily changed. Classification sensitivity, specificity and memory requirements may be tuned by varying *k*-mer length (`-k`), window size (`-w`), absolute match threshold (`-a`) and relative match threshold (`-r`) . Minimizer `k` and `w` are chosen at index time, while the match thresholds can be chosen at filter time. Matching sequences are those that share enough distinct minimizers with the indexed query to exceed *both* the absolute threshold (`-a`, default 2 shared minimizers) and the relative threshold (`-r`, default 0.01 [1%] shared minimizers). Paired sequences are fully supported: a match in either mate causes both mates in the pair to be matched; `deacon filter` retains  matches by default (search mode) and discards matches in `--deplete` mode. Deacon reports live filtering performance during execution and optionally writes a JSON `--summary` upon completion. Sequences can optionally be renamed using `--rename` for privacy and smaller file sizes. Gzip, zst and xz compression formats are natively supported and detected by file extension. Other source formats can be converted to FASTA or FASTQ for consumption by Deacon via stdin.

Benchmarks for panhuman host depletion of complex microbial metagenomes are described in a [preprint](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1). Among tested approaches, Deacon with the panhuman-1 (*k*=31, w=15) index exhibited the highest balanced accuracy for both long and short simulated reads. Deacon was however less specific than Hostile for short reads.

> [!IMPORTANT]
> Deacon is actively developed. Take note of software and index version(s) used in order to guarantee reproducibility of your results. Carefully review the CHANGELOG when updating. Versions 0.7.0 and 0.11.0 introduced backwards incompatible index formats. Please report any problems you encounter by creating an issue or using the email address in my profile.

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

### Indexing search queries

Use `deacon index build` to index a search query. For panhuman host depletion, skip this step and download the validated `panhuman-1` index from the table below.

```shell
deacon index build chm13v2.fa > human.k31w15.idx

# Discard low complexity minimizers
deacon index build -e 0.5 chm13v2.fa > human.k31w15e5.idx
```

*N.B. Indexing a human genome takes <30s and uses 18GB of RAM. Filtering uses 5GB of RAM.*

#### Prebuilt indexes

|                          Name & URL                          |                         Composition                          | Minimizers  | Subtracted minimizers | Size  | Date    |
| :----------------------------------------------------------: | :----------------------------------------------------------: | ----------- | --------------------- | ----- | ------- |
| **`panhuman-1` (*k*=31, *w*=15)** [Cloud](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w15.idx), [Zenodo](https://zenodo.org/records/17288185) | [HPRC Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index) ∪ [`CHM13v2.0`](https://www.ncbi.nlm.nih.gov/assembly/11828891) ∪ [`GRCh38.p14`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40) - bacteria (FDA-ARGOS) - viruses (RefSeq) | 409,907,949 | 20,671 (**0.0050%**)  | 3.3GB | 2025-04 |
| **`panmouse-1` (k=31, w=15, e=0.5)** [Cloud](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panmouse-1.k31w15e05.idx) | [`GRCm39`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27) ∪ [`PRJEB47108`](https://www.ebi.ac.uk/ena/browser/view/PRJEB47108?show=sequences) - bacteria (FDA-ARGOS) - viruses (RefSeq) | 548,328,389 | 8,243 (**0.0015%**)   | 4.4GB | 2025-08 |

#### Index compatibility

Deacon `0.11.0` and above uses index format version 3. Using version 3 indexes with older Deacon versions and vice versa triggers an error. Prebuilt indexes in legacy formats are therefore archived in object storage for reproducibility. Should you wish to download indexes in legacy formats, replace the `/3/` in any prebuilt index download URL with either `/2/` or `/1/`  accordingly. Indexes are not however backported to Deacon versions predating their existence.

- Deacon **`0.11.0`** and above uses index format version **`3`**
- Deacon **`0.7.0`** through to **`0.10.0`** used index format version **`2`**

- Deacon **`0.1.0`** through to **`0.6.0`** used index format version **`1`**

### Filtering

The main command `deacon filter` accepts an index path followed by up to two FASTA/FASTQ file paths, depending on whether input sequences originate from stdin, a single file, or paired input files. Paired queries are supported as either separate files or interleaved stdin, and written interleaved to either stdout or file, or else to separate paired output files. For paired reads, distinct minimizer hits originating from either mate are counted. By default, input sequences must meet both an absolute threshold of 2 minimizer hits (`-a 2`) and a relative threshold of 1% of minimizers (`-r 0.01`) to pass the filter. Filtering can be inverted for e.g. host depletion using the `--deplete` (`-d`) flag. Gzip, Zstandard, and xz compression formats are detected automatically by file extension.

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
Build, compose and inspect minimizer indexes

Usage: deacon index <COMMAND>

Commands:
  build  Index minimizers contained within a fastx file
  info   Show index information
  union  Combine multiple minimizer indexes (A ∪ B…)
  diff   Subtract minimizers in one index from another (A - B)
  dump   Dump minimizer index to fasta
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
          K-mer length used for indexing (k+w-1 must be <= 96 and odd) [default: 31]
  -w <WINDOW_SIZE>
          Minimizer window size used for indexing [default: 15]
  -o, --output <OUTPUT>
          Path to output file (stdout if not specified)
  -t, --threads <THREADS>
          Number of execution threads (0 = auto) [default: 8]
  -q, --quiet
          Suppress sequence header output
  -e, --entropy-threshold <ENTROPY_THRESHOLD>
          Minimum scaled entropy threshold for k-mer filtering (0.0-1.0) [default: 0.0]
  -h, --help
          Print help
```

## Building custom indexes

Building custom Deacon indexes is fast. Nevertheless, when indexing many large genomes, it may be worthwhile separately indexing and subsequently combining indexes into one succinct index. Combine distinct minimizers from multiple indexes using `deacon index union`. Similarly, use `deacon index diff` to subtract the minimizers contained in one index from another. This can be helpful for e.g. eliminating shared minimizers between the target and host genomes when building custom (non-human) indexes for host depletion.

- Use `deacon index union 1.idx 2.idx 3.idx… > 1+2+3.idx` to succinctly combine two (or more!) deacon indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 1.idx from 2.idx. Useful for masking out shared minimizer content between e.g. target and host genomes.
- From version `0.7.0`, `deacon index diff` also supports subtracting minimizers from an index using a fastx file or stream, e.g. `deacon index diff 1.idx 2.fa.gz > 1-2.idx` or `zcat *.fa.gz | deacon index diff 1.idx - > 1-2.idx`.


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

## Server mode

From version `0.11.0`, it is possible to eliminate index loading overhead at the start of each filter operation by preloading the index in the memory of a local server process. Subsequent filtering commands with `--use-server` are executed by the server process using a UNIX socket. Having started a server process, the index of the first filtering command it receives persists in memory for the life of that server process, enabling subsequent filter commands to be served rapidly without hash set construction overhead.

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

## Citation

[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

>  Bede Constantinides, John Lees, Derrick W Crook. "Deacon: fast sequence filtering and contaminant depletion" *bioRxiv* 2025.06.09.658732, https://doi.org/10.1101/2025.06.09.658732 

Please also consider citing the SimdMinimizers paper:

> Ragnar Groot Koerkamp, Igor Martayan. "SimdMinimizers: Computing random minimizers, *fast*" *bioRxiv* 2025.01.27.634998, https://doi.org/10.1101/2025.01.27.634998 

