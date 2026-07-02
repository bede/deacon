# Deacon for Python

Python bindings for [Deacon](https://github.com/bede/deacon), enabling fast multithreaded DNA sequence filtering for e.g. host pangenome depletion using Python code. These bindings load an index once, allowing subsequent filtering runs with low latency. Deacon's complete functionality is currently only available using the [Rust/CLI version of Deacon](https://github.com/bede/deacon).

## Installation

```bash
uv install deacon
```

## Quickstart

```python
from deacon import Index

index = Index("panhuman-1.k31w15.idx")
stats = index.filter(
    fastq_path,
    deplete=True,
    rename=True,
    output=fastq_path.replace(".fastq.gz", ".clean.fastq.gz")
)
print(stats["seqs_in"], stats["seqs_out"])
```

## `Index()`

Load a minimizer index or probabilistic filter from disk. The resulting object may be reused across many `filter` calls.

```python
index = Index("panhuman-1.k31w15.idx", complexity_threshold=None)
```

Pass `complexity_threshold` (0.0–1.0, e.g. `0.9`) to discard low-complexity index minimizers once at load using kdust; the filtered set is then reused across all `filter` calls. Not supported for `bff` (binary fuse filter) indexes.

## `Index.fetch()`

Download a prebuilt index, then load and return it (a static method, so `Index.fetch(...)` returns an `Index`). `output` is the local path to save to; when omitted it defaults to `"{name}.k{k}w{w}.idx"` in the working directory. The index is downloaded on every call — there is no local cache, so an existing file at that path is overwritten.

```python
index = Index.fetch("panhuman-1", k=31, w=15, output=None, complexity_threshold=None)
```

## `Index.info()`

`Index.info(index_path)` Returns a `dict` of index metadata:

| Key | Meaning |
| --- | --- |
| `k` | *k*-mer length |
| `w` | minimizer window size |
| `format` | `exact-u64`, `exact-u128`, or `bff` (binary fuse filter) |
| `count` | number of keys in index (fingerprint slot count for `bff`) |

## `Index.filter()`

Filter a FASTA/FASTQ file or file pair against the index and return a `dict` of summary statistics. Auto-detects `.gz`/`.zst`/`.xz` compression on both input and output based on file extension. The Python GIL is released while filtering meaning calls benefit from multithreading. Refer to the [main Deacon readme](https://github.com/bede/deacon) for more detailed usage examples.

```python
def filter(
    fastq,                   # input path (FASTA/FASTQ, optionally .gz/.zst/.xz)
    fastq2=None,             # second mate for paired reads
    interleaved=False,       # treat fastq as interleaved paired reads (cannot combine with fastq2)
    deplete=False,           # False = search (keep matches); True = deplete (remove matches)
    rename=False,            # replace read names with sequential integers
    rename_random=False,     # replace read names with random strings
    output=None,             # output path; None writes to stdout
    output2=None,            # second output path for paired reads
    abs_threshold=2,         # min absolute minimizer hits to call a match
    rel_threshold=0.01,      # min proportion of minimizers hitting to call a match
    prefix_length=0,         # only use the first N bp of each read (0 = whole read)
    output_fasta=False,      # emit FASTA instead of FASTQ
    threads=8,               # worker threads for filtering
    compression_level=2,     # output compression level
    compression_threads=0,   # threads for output compression (0 = auto)
    debug=False,             # verbose per-read debug output
    quiet=True,              # suppress progress/log output on stderr
) -> dict
```

**Modes.** With `deplete=False` (the default, *search* mode) reads that match the index are kept; with `deplete=True` reads that match are removed (host depletion). A read is a match only when it clears **both** thresholds: at least `abs_threshold` minimizer hits **and** at least `rel_threshold` of its minimizers hitting the index.

**Output.** When `output` is `None` the filtered records are written to stdout. To count without keeping the filtered sequences, pass `output="/dev/null"`. Statistics are returned regardless.

**Return value.** A `dict` including the run configuration (`version`, `index`, `input`/`input2`, `output`/`output2`, `k`, `w`, `abs_threshold`, `rel_threshold`, `prefix_length`, `deplete`, `rename`, `rename_random`) and the results:

| Key | Meaning |
| --- | --- |
| `seqs_in`, `seqs_out`, `seqs_removed` | sequence counts |
| `seqs_out_proportion`, `seqs_removed_proportion` | sequence proportions |
| `bp_in`, `bp_out`, `bp_removed` | base-pair counts |
| `bp_out_proportion`, `bp_removed_proportion` | base-pair proportions |
| `time` | wall-clock seconds |
| `seqs_per_second`, `bp_per_second` | throughput (filtering only) |
| `seqs_per_second_total`, `bp_per_second_total` | throughput (including load/IO) |
