[![Crates.io version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)
[![Conda version](https://img.shields.io/conda/v/bioconda/deacon?style=flat-square&label=bioconda&color=blue)](https://anaconda.org/bioconda/deacon)
[![Crates.io downloads](https://img.shields.io/crates/d/deacon?color=orange&label=crates.io%20downloads&style=flat-square)](https://crates.io/crates/deacon)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/deacon.svg?style=flat-square&label=conda%20downloads&color=blue)](https://anaconda.org/bioconda/deacon)
[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

<div align="center"><img src="https://github.com/bede/deacon/blob/main/deacon.png?raw=true" width="180" alt="Logo"></div>

# Deacon

Deacon filters DNA sequences in FASTA/Q files and streams using SIMD-accelerated minimizer comparison with query sequence(s), emitting either matching sequences (**search mode**), or sequences without matches (**deplete mode**). Sequences match when they share enough distinct minimizers with the indexed query to exceed chosen absolute and relative thresholds. Query size has little impact on filtering speed, enabling ultrafast search and depletion with gene-, genome- and pangenome-scale queries using a laptop. Deacon filters uncompressed FASTA/Q at **gigabases per second** on recent AMD, Intel (`x86_64`), and Apple `arm64` systems. Built with panhuman host depletion in mind—yet broadly useful for searching large sequence collections—Deacon delivers [leading classification accuracy](https://doi.org/10.1101/2025.06.09.658732) for host depletion and unrivalled speed using 5GB of RAM.

Default parameters are carefully chosen but easily changed. Classification sensitivity, specificity and memory requirements may be tuned by varying *k*-mer length, window size, absolute match threshold and relative match threshold. Minimizer `k` and `w` are chosen at query index time, while the match thresholds can be chosen at filter time. Matching sequences are those that share enough distinct minimizers with the indexed query to exceed *both* the absolute threshold (default 2 shared minimizers) and the relative threshold (default 0.01 [1%] shared minimizers). For paired sequences, hits in either mate counts towards a single match threshold for the pair. Deacon reports filtering performance during execution and optionally writes a JSON upon completion. Sequences can optionally be renamed for privacy and smaller file sizes. Deacon natively handles gz, zst and xz compression formats, detected by file extension.

Benchmarks for panhuman host depletion of complex microbial metagenomes are described in a [preprint](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1). Deacon with the `panhuman-1` (*k*=31, w=15) index exhibited the highest balanced accuracy for both long and short simulated reads. Deacon was less specific only than Hostile for short reads.

## Use cases

- Depletion of human or other host genome sequences in FASTQ reads or streams.
- Ultrafast binary classification of genes, genomes or pangenomes in terabase genome catalogues like [AllTheBacteria](https://allthebacteria.org/) without tedious pre-indexing.

## Note
For use as a standalone CLI tool, please refer to the non-python package https://github.com/bede/deacon

## Install
```bash
pip install deacon
```

## Quickstart

### Ultrafast panhuman host depletion

```python
import deacon

# Download validated 3GB human pangenome index to the current directory
deacon.index_fetch()

# Deplete long reads
deacon.filter("panhuman-1.k31w15.idx", "reads.fq", "filt.fq", deplete=True)

# Deplete short paired reads
deacon.filter("panhuman-1.k31w15.idx", "reads.r1.fq.gz", "filt.r1.fq.gz", input2_path="reads.r2.fq.gz", output2_path="filt.r2.fq.gz", deplete=True)
```

### Ultrafast gene/genome/pangenome search

```python
import deacon

deacon.index_build("amr-genes.fa", "amr-genes.idx")

deacon.filter("amt-genes.idx", "AllTheBacteria.fa.zst", "hits.fa")
```

*N.B. Indexing a 3Gbp human genome takes ~30s using 18GB of RAM with default parameters. Filtering uses 5GB.*

## Prebuilt indexes

Prebuilt pangenome indexes are provided for human and mouse host classification and depletion. These can be downloaded using the links below, or with `deacon.index_fetch(index_name="panhuman-1")`.

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

#### Examples

```python
import deacon

# Keep only sequences matching a collection of genes
deacon.index_build("genes.fa", "genes.idx")
deacon.filter("genes.idx", "sequences.fa.gz", "matches.fa.gz")

# Host depletion using the panhuman-1 index
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fq.gz", deplete=True)

# High sensitivity host depletion with absolute threshold of 1 and no relative threshold
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fq.gz", deplete=True, abs_threshold=1, rel_threshold=0)

# High specificity 10% relative match threshold
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fq.gz", deplete=True, rel_threshold=0.1)

# Zstandard compression
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.zst", "filt.fq.zst", deplete=True)

# Save summary JSON
summary = deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fq.gz", deplete=True, summary_path="summary.json")
# `deacon.filter` also returns a summary object
print(summary.seqs_in, summary.seqs_out)

# Replace read headers with incrementing integers
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fq", deplete=True, rename=True)

# Only look for minimizer hits inside the first 1000bp per record
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fq", deplete=True, prefix_length=1000)

# Output FASTA regardless of input format (discards quality scores)
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fa", deplete=True, output_fasta=True)

# Debug mode: see sequences with minimizer hits in stderr
deacon.filter("panhuman-1.k31w15.idx", "reads.fq.gz", "filt.fq.gz", deplete=True, debug=True)
```

> [!NOTE]
>
> `deacon.filter()` uses 8 threads by default. Using more threads (e.g.  `deacon.filter(..., threads=16)`) can accelerate filtering given sufficient resources, especially with uncompressed sequences whose processing is not rate limited by decompression. Since version `0.13.0`, Deacon writes gzipped output files in parallel, providing particular practical benefit for gzipped paired reads. If output file(s) ending in `.gz` are detected, total `threads` are allocated 1:1 to compression and filtering tasks respectively. Gzip compression thread allocation can be overriden with the `compression-threads` kwarg.

### Indexing

```python
# Index one FASTA/FASTQ file
import deacon
deacon.index_build("genome.fa.gz", "genome.idx")
```

`deacon.index_build()` accepts either a FASTA or a FASTQ file enabling convenient indexing of compressed sequences in one or many files with a single step. Indexing a human genome takes a few seconds. Indexing uses 2-4x as much RAM as filtering. For indexing large collections approaching terabase scale—such as mammalian pangenomes—it may be practical to index genomes individually in parallel and later combine them using the `deacon.index_union()` set operation, described below.

#### Set operations

A differentiating feature of Deacon is the ease of combining, subtracting and intersecting minimizer indexes. For example, `deacon.index_diff()`can be used to subtract shared minimizers between target and host genomes when building custom indexes for host depletion.

- Use `deacon.index_union(["1.idx", "2.idx", "3.idx"], "1+2+3.idx")` to succinctly combine two or more indexes.
- Use `deacon.index_diff("1.idx", "2.idx", "1-2.idx")` to subtract minimizers in 2.idx from 1.idx. Useful for masking out shared minimizer content between e.g. target and host genomes.
- Use `deacon.index_intersect(["1.idx", "2.idx"], "1∩2.idx")` to find the intersection of minimizers in two or more indexes.

#### Inspecting indexes

- Use `deacon.index_info("1.idx")` to display index information including minimizer *k* and *w* parameters, number of minimizers, and index format version.
- Use `deacon.index_dump("1.idx", "1.fa")` to dump a minimizer index to FASTA.


## Filtering summary statistics

Use the `--summary_path=summary.json` flag on `deacon.filter()` to save detailed filtering statistics:
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

## Citation

[![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red?&style=flat-square)](https://doi.org/10.1101/2025.06.09.658732)

>  Bede Constantinides, John Lees, Derrick W Crook. "Deacon: fast sequence filtering and contaminant depletion" *bioRxiv* 2025.06.09.658732, https://doi.org/10.1101/2025.06.09.658732

Please also consider citing the SimdMinimizers paper:

> Ragnar Groot Koerkamp, Igor Martayan. "SimdMinimizers: Computing random minimizers, *fast*" *bioRxiv* 2025.01.27.634998, https://doi.org/10.1101/2025.01.27.634998

