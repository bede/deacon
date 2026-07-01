#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/deacon-wasm-paired-parity.XXXXXX")"
pkg_dir="$work_dir/pkg"
native_r1="$work_dir/native.r1.fastq"
native_r2="$work_dir/native.r2.fastq"
wasm_r1="$work_dir/wasm.r1.fastq"
wasm_r2="$work_dir/wasm.r2.fastq"
config_path="$repo_root/.cargo/config-wasm.toml"
index_path="$repo_root/data/mn908947.idx"
reads1_path="$repo_root/tests/data/test_small_1.fastq.gz"
reads2_path="$repo_root/tests/data/test_small_2.fastq.gz"

cleanup() {
  rm -rf "$work_dir"
}
trap cleanup EXIT

mkdir -p "$pkg_dir"
export LC_ALL=C
export LANG=C

echo "Running native paired filter..."
cargo run --release -- filter -d -t 1 "$index_path" "$reads1_path" "$reads2_path" \
  -o "$native_r1" -O "$native_r2"

echo "Building Node-target WASM bundle..."
export CARGO_BUILD_CONFIG="$config_path"
wasm-pack build "$repo_root/deacon-wasm" --target nodejs --no-opt --out-dir "$pkg_dir"

echo "Running WASM paired filter..."
node "$repo_root/scripts/wasm_paired_filter_to_file.mjs" \
  --pkg "$pkg_dir" \
  --index "$index_path" \
  --reads1 "$reads1_path" \
  --reads2 "$reads2_path" \
  --output1 "$wasm_r1" \
  --output2 "$wasm_r2" \
  --deplete

if ! cmp -s "$native_r1" "$wasm_r1" || ! cmp -s "$native_r2" "$wasm_r2"; then
  echo "Native/WASM paired output mismatch" >&2
  echo "  native R1: $(shasum -a 256 "$native_r1" | awk '{print $1}')" >&2
  echo "  wasm R1:   $(shasum -a 256 "$wasm_r1" | awk '{print $1}')" >&2
  echo "  native R2: $(shasum -a 256 "$native_r2" | awk '{print $1}')" >&2
  echo "  wasm R2:   $(shasum -a 256 "$wasm_r2" | awk '{print $1}')" >&2
  exit 1
fi

echo "Native/WASM paired parity check passed"
echo "  R1 bytes: $(wc -c < "$native_r1" | tr -d ' ')"
echo "  R2 bytes: $(wc -c < "$native_r2" | tr -d ' ')"

# Second pass: --rename + --fasta (-R -f) must also match byte-for-byte.
native_rf_r1="$work_dir/native.rf.r1.fasta"
native_rf_r2="$work_dir/native.rf.r2.fasta"
wasm_rf_r1="$work_dir/wasm.rf.r1.fasta"
wasm_rf_r2="$work_dir/wasm.rf.r2.fasta"
echo "Running native paired filter (-R -f)..."
cargo run --release -- filter -d -t 1 -R -f "$index_path" "$reads1_path" "$reads2_path" \
  -o "$native_rf_r1" -O "$native_rf_r2"
echo "Running WASM paired filter (--rename --fasta)..."
node "$repo_root/scripts/wasm_paired_filter_to_file.mjs" \
  --pkg "$pkg_dir" \
  --index "$index_path" \
  --reads1 "$reads1_path" \
  --reads2 "$reads2_path" \
  --output1 "$wasm_rf_r1" \
  --output2 "$wasm_rf_r2" \
  --deplete --rename --fasta

if ! cmp -s "$native_rf_r1" "$wasm_rf_r1" || ! cmp -s "$native_rf_r2" "$wasm_rf_r2"; then
  echo "Native/WASM paired output mismatch (rename+fasta)" >&2
  echo "  native R1: $(shasum -a 256 "$native_rf_r1" | awk '{print $1}')" >&2
  echo "  wasm R1:   $(shasum -a 256 "$wasm_rf_r1" | awk '{print $1}')" >&2
  echo "  native R2: $(shasum -a 256 "$native_rf_r2" | awk '{print $1}')" >&2
  echo "  wasm R2:   $(shasum -a 256 "$wasm_rf_r2" | awk '{print $1}')" >&2
  exit 1
fi
echo "Native/WASM paired rename+fasta parity check passed"
