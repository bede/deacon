# Deacon Python bindings

Enables index reuse between filtering sessions.

```python
from deacon import Index

index = Index("panhuman-1.k31w15.idx")
stats = [
    index.filter(fq, deplete=True, rename=True,
                 output=fq.replace(".fastq.gz", ".clean.fastq.gz"))
    for fq in input_fastqs
]
print(stats[0]["seqs_in"], stats[0]["seqs_out"])
```

`filter` accepts single or paired input (`fastq2=`, `output2=`), auto-detects `.gz`/`.zst`/
`.xz` output compression from the extension, and returns a `dict` of summary statistics.
`from deacon import filter` also exposes `filter(index, ...)` as a free function.
