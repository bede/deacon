# Deacon for Python

Python bindings for [Deacon](https://github.com/bede/deacon), enabling fast multithreaded DNA sequence filtering for e.g. host pangenome depletion using Python code. These bindings load an index once, allowing subsequent filtering runs with low latency. Deacon's complete functionality is currently only available using the [Rust/CLI version of Deacon](https://github.com/bede/deacon).

## Installation

```bash
uv pip install deacon
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

`path` is positional-only and accepts a string or path-like object. `complexity_threshold` is keyword-only.

## `Index.fetch()`

Download a prebuilt index, then load and return it (a static method, so `Index.fetch(...)` returns an `Index`). `output` is the local path to save to; when omitted it defaults to `"{name}.k{k}w{w}.idx"` in the working directory. The index is downloaded on every call — there is no local cache, so an existing file at that path is overwritten.

```python
index = Index.fetch(
    name="panhuman-1",
    k=31,
    w=15,
    output=None,
    complexity_threshold=None,
)
```

All `fetch()` arguments are keyword-only.

## `Index.info()`

`index.info()` returns a `dict` of index metadata:

| Key | Meaning |
| --- | --- |
| `k` | *k*-mer length |
| `w` | minimizer window size |
| `format` | `exact-u64`, `exact-u128`, or `bff` (binary fuse filter) |
| `count` | number of minimizers/keys represented by the index |

## `Index.filter()`

Filter FASTA, FASTQ, or CBQ input against the index and return a `dict` of summary statistics. FASTA/FASTQ compression (`.gz`, `.zst`, or `.xz`) is detected automatically. A `.cbq` output writes CBQ; `.cba` writes quality-free CBQ. The Python GIL is released while filtering, so calls benefit from multithreading. Refer to the [main Deacon readme](https://github.com/bede/deacon) for more detailed usage examples.

```python
def filter(
    input,                   # positional-only FASTA, FASTQ, or CBQ path
    /,
    *,                       # every remaining argument is keyword-only
    input2=None,             # second FASTA/FASTQ mate
    interleaved=False,       # treat input as interleaved pairs (cannot combine with input2)
    check_pairs=False,       # validate paired read names
    deplete=False,           # False = search (keep matches); True = deplete (remove matches)
    rename=False,            # replace read names with sequential integers
    output=None,             # output path; None writes FASTA/FASTQ to stdout
    output2=None,            # second output path for paired reads
    inverse_output=None,     # path for discarded records; None discards them
    inverse_output2=None,    # second output path for discarded paired reads
    summary=None,            # optional JSON summary output path
    abs_threshold=2,         # min absolute minimizer hits to call a match
    rel_threshold=0.01,      # min proportion of minimizers hitting to call a match
    prefix_length=0,         # only use the first N bp of each read (0 = whole read)
    discard_quality=False,   # discard quality scores
    ordered=False,           # preserve input record ordering (deterministic, slightly slower)
    threads=8,               # worker threads for filtering
    compression_level=2,     # output compression level
    compression_threads=0,   # threads for output compression (0 = auto)
    cbq_block_size=16,        # CBQ output block size in MiB (1-1024)
    quiet=True,              # suppress progress/log output on stderr
    debug=False,             # verbose per-read debug output
) -> dict
```

**Modes.** With `deplete=False` (the default, *search* mode) reads that match the index are kept; with `deplete=True` reads that match are removed (host depletion). A read is a match only when it clears **both** thresholds: at least `abs_threshold` minimizer hits **and** at least `rel_threshold` of its minimizers hitting the index.

**Paired and CBQ input.** Use `input2=` for separate FASTA/FASTQ mates or `interleaved=True` for an interleaved stream. `check_pairs=True` validates Illumina CASAVA or `/1` and `/2` names. CBQ stores pairing internally, so CBQ input cannot be combined with `input2` or `interleaved`; paired CBQ output uses one `output` and cannot use `output2`.

**Output.** When `output` is `None` retained records are written to stdout. Set `inverse_output=` to write discarded records. Primary and inverse outputs must both be FASTX or both be CBQ/CBA; compression and paired-file layout may differ. For paired FASTX, set `inverse_output2=` to use separate mate files. Otherwise the inverse output is interleaved. To count without keeping retained sequences, pass `output="/dev/null"`. Pass `summary=` to write the same statistics returned by the call as JSON. For CBQ output, `cbq_block_size` is a lower bound. A larger input block size is preserved.

**Return value.** A `dict` including the run configuration (`version`, `index`, `input`/`input2`, `output`/`output2`, `inverse_output`/`inverse_output2`, `k`, `w`, `abs_threshold`, `rel_threshold`, `prefix_length`, `deplete`, `rename`, `ordered`, `check_pairs`) and the results:

| Key | Meaning |
| --- | --- |
| `seqs_in`, `seqs_out`, `seqs_removed` | sequence counts |
| `seqs_out_proportion`, `seqs_removed_proportion` | sequence proportions |
| `bp_in`, `bp_out`, `bp_removed` | base-pair counts |
| `bp_out_proportion`, `bp_removed_proportion` | base-pair proportions |
| `time` | wall-clock seconds |
| `seqs_per_second`, `bp_per_second` | throughput (filtering only) |
| `seqs_per_second_total`, `bp_per_second_total` | throughput (including load/IO) |
