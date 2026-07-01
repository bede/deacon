#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage: scripts/bench-filter-ab.sh <baseline-deacon-bin> <candidate-deacon-bin> [iters]

Runs native A/B filtering benchmarks with 8 filter threads against:
  index: data/panhuman-1.k31w15.idx
  reads: ../host-depletion-bench/data/rsviruses17900.fasta

Filtered output is always written to /dev/null.
USAGE
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
baseline_bin="$1"
candidate_bin="$2"
iters="${3:-3}"
index_path="$repo_root/data/panhuman-1.k31w15.idx"
reads_path="$repo_root/../host-depletion-bench/data/rsviruses17900.fasta"

if [[ ! -x "$baseline_bin" ]]; then
  echo "Baseline binary is not executable: $baseline_bin" >&2
  exit 1
fi
if [[ ! -x "$candidate_bin" ]]; then
  echo "Candidate binary is not executable: $candidate_bin" >&2
  exit 1
fi
if [[ ! -f "$index_path" ]]; then
  echo "Missing index: $index_path" >&2
  exit 1
fi
if [[ ! -f "$reads_path" ]]; then
  echo "Missing reads: $reads_path" >&2
  exit 1
fi

run_one() {
  local label="$1"
  local bin="$2"
  local iter="$3"
  local log
  log="$(mktemp "${TMPDIR:-/tmp}/deacon-bench-${label}-${iter}.XXXXXX")"

  echo "BENCH label=$label iter=$iter output=/dev/null"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    if ! /usr/bin/time -l "$bin" filter -d -t 8 -q "$index_path" "$reads_path" -o /dev/null 2>"$log"; then
      cat "$log" >&2
      rm -f "$log"
      return 1
    fi
    awk -v label="$label" -v iter="$iter" '
      /real/ { real=$1 }
      /maximum resident set size/ { rss=$1 }
      END { printf("RESULT label=%s iter=%s real_secs=%s max_rss_bytes=%s\n", label, iter, real, rss) }
    ' "$log"
  else
    if ! /usr/bin/time -v "$bin" filter -d -t 8 -q "$index_path" "$reads_path" -o /dev/null 2>"$log"; then
      cat "$log" >&2
      rm -f "$log"
      return 1
    fi
    awk -v label="$label" -v iter="$iter" '
      /Elapsed \(wall clock\) time/ { elapsed=$NF }
      /Maximum resident set size/ { rss_kb=$NF }
      END { printf("RESULT label=%s iter=%s elapsed=%s max_rss_kb=%s\n", label, iter, elapsed, rss_kb) }
    ' "$log"
  fi
  rm -f "$log"
}

echo "Benchmark input:"
echo "  index=$index_path"
echo "  reads=$reads_path"
echo "  output=/dev/null"
echo "  threads=8"
echo "  iters=$iters"

for ((i = 1; i <= iters; i += 2)); do
  run_one baseline "$baseline_bin" "$i"
  if (( i + 1 <= iters )); then
    run_one baseline "$baseline_bin" "$((i + 1))"
  fi
  run_one candidate "$candidate_bin" "$i"
  if (( i + 1 <= iters )); then
    run_one candidate "$candidate_bin" "$((i + 1))"
  fi
done
