[![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon) [![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon) ![CI status](https://img.shields.io/github/actions/workflow/status/bede/deacon/test.yml?style=flat-square) [![Downloads](https://img.shields.io/conda/dn/bioconda/deacon.svg?style=flat-square)](https://anaconda.org/bioconda/deacon) [![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red)](https://doi.org/10.1101/2025.06.09.658732)


# Deacon

<div align="center"><img src="deacon.png" width="180" alt="Logo"></div>

Fast minimizer-based search and depletion of FASTA/FASTQ files and streams. Default parameters balance sensitivity and specificity for microbial (meta)genomic host depletion, for which a validated prebuilt index is available. Classification sensitivity, specificity and required memory can be tuned by varying *k*-mer length (`-k`), minimizer window size (`-w`), and the number or proportion of required index matches (`-m`) per query. Minimizer `k` and `w`  are chosen at index time, while the match threshold `m` can be varied at filter time. `m` can be specified either as a minimum integer of minimizer matches (default is 2), else a minimum proportion of minimizer hits between 0.0 and 1.0.

Building on [simd-minimizers](https://github.com/rust-seq/simd-minimizers), Deacon is capable of filtering long reads at >250Mbp/s (Apple M4) and indexing a human genome in <30s. Short and/or paired reads are fully supported albeit more slowly. Peak memory usage during filtering is 5GB for the default panhuman index. Sequences can optionally be renamed for privacy and smaller file sizes. Deacon reports filtering performance during execution and optionally writes a JSON summary on completion.

Benchmarks for panhuman host depletion of complex microbial metagenomes are described in a [preprint](https://www.biorxiv.org/content/10.1101/2025.06.09.658732v1). Among tested approaches, Deacon with the panhuman-1 (*k*=31, w=15) index exhibited the highest balanced accuracy for both long and short simulated reads. Deacon was however less specific than Hostile for short reads.

> [!IMPORTANT]
> Deacon is not yet stable. 0.5.0 for instance introduced major CLI changes. Please carefully review the CHANGELOG when upgrading. Please report problems by creating an issue or using the email address in my profile

## Install

### conda/mamba/pixi  [![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon)

```bash
conda install -c bioconda deacon
```

### cargo [![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)

```bash
cargo install deacon
```

## Usage

### Indexing

Build indexes with `deacon index build`. For human host depletion, the prebuilt validated panhuman index is recommended, available for download below. Object storage is provided by the [ModMedMicro research unit](https://www.expmedndm.ox.ac.uk/modernising-medical-microbiology) at the University of Oxford.

```shell
deacon index build chm13v2.fa > human.k31w15.idx
```

#### Prebuilt indexes

|                           Name/URL                           |                         Composition                          | Minimizers  | Subtracted minimizers | Size  | Date    |
| :----------------------------------------------------------: | :----------------------------------------------------------: | ----------- | --------------------- | ----- | ------- |
| [**panhuman-1 (*k*=31, *w*=15)**](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/panhuman-1.k31w15.idx) | ([HPRC Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index) ∪ [CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) ∪ [GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40)) - bacteria (FDA-ARGOS)  - viruses (RefSeq) | 409,914,298 | 20,741 (**0.0051%**)  | 3.7GB | 2025-04 |

### Filtering

The command `deacon filter` accepts an index path followed by up to two query FASTA/FASTQ file paths, depending on whether query sequences originate from stdin, a single file, or paired input files. Paired queries are supported as either separate files or interleaved stdin, and written interleaved to either stdout or file, or else to separate paired output files. For paired reads, distinct minimizer hits originating from either mate are counted. By default, query sequences with two or more minimizer hits to the index (`-m 2`) pass the filter. Filtering can be inverted for e.g. host depletion using the `--deplete` (`-d`) flag. Gzip (.gz) and Zstandard (.zst) compression formats are detected automatically by file extension. Use Zstandard compression rather than Gzip where possible for best performance.

**Examples**

```bash
# Keep only human sequences
deacon filter panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz

# Host depletion using the panhuman-1 index
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz

# Host depletion using stdin and stdout
zcat reads.fq.gz | deacon filter -d panhuman-1.k31w15.idx > filt.fq

# Faster Zstandard compression
deacon filter -d panhuman-1.k31w15.idx reads.fq.zst -o filt.fq.zst

# More sensitive match threshold of at least 1 minimizer hit
deacon filter -d -m 1 panhuman-1.k31w15.idx reads.fq.gz > filt.fq.gz

# More specific match threshold of 50% minimizer hits (minimum 1)
deacon filter -d -m 0.5 panhuman-1.k31w15.idx reads.fq.gz > filt.fq.gz

# Deplete paired reads
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz > filt12.fq
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz -o filt.r1.fq.gz -O filt.r2.fq.gz
zcat r12.fq.gz | deacon filter -d panhuman-1.k31w15.idx - - > filt12.fq

# Save summary JSON
deacon filter -d panhuman-1.k31w15.idx reads.fq.gz -o filt.fq.gz -s summary.json
```

## Building custom indexes

Building custom Deacon indexes is quite fast. Nevertheless, when indexing many large genomes, it may be worthwhile separately indexing and subsequently combining indexes into one succinct index. Combine distinct minimizers from multiple indexes using `deacon index union`. Similarly, use `deacon index diff` to subtract the minimizers contained in one index from another. This can be helpful  for e.g. eliminating shared minimizers between the target and host genomes when building custom (non-human) indexes for host depletion.

- Use `deacon index union 1.idx 2.idx 3.idx… > 1+2+3.idx` to succinctly combine two (or more!) deacon indexes.
- Use `deacon index diff host.idx fungi.idx > 1-2.idx` to subtract minimizers in fungi.idx from host.idx. Useful for masking out shared minimizer content between e.g. target and host genomes.

For best performance, set the `--capacity` argument of `deacon index build` to a number of minimizers in millions greater than that you expect your index to contain. Setting this too low will cause delays during indexing for hash table resizing.

## Filtering summary statistics

Use `-s summary.json` to save detailed filtering statistics:
```json
{
  "version": "deacon 0.5.0",
  "index": "panhuman-1.k31w15.idx",
  "input1": "HG02334.1m.fastq.gz",
  "input2": null,
  "output": "-",
  "output2": null,
  "k": 31,
  "w": 21,
  "match_threshold": "2",
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

 [![biorXiv preprint](https://img.shields.io/badge/biorXiv-10.1101/2025.06.09.658732-red)](https://doi.org/10.1101/2025.06.09.658732)

>  Bede Constantinides, John Lees, Derrick W Crook. "Deacon: fast sequence filtering and contaminant depletion" *bioRxiv* 2025.06.09.658732, https://doi.org/10.1101/2025.06.09.658732 

Please also consider citing the SimdMinimizers paper:

> Ragnar Groot Koerkamp, Igor Martayan. "SimdMinimizers: Computing random minimizers, *fast*" *bioRxiv* 2025.01.27.634998, https://doi.org/10.1101/2025.01.27.634998 

