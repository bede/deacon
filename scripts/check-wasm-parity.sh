#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/deacon-wasm-parity.XXXXXX")"
pkg_dir="$work_dir/pkg"
native_out="$work_dir/native.fastq"
wasm_out="$work_dir/wasm.fastq"
config_path="$repo_root/.cargo/config-wasm.toml"
index_path="$repo_root/data/test-human-mit.k31w15.idx"
reads_path="$repo_root/data/HG02334.100MB.fastq.gz"

cleanup() {
  rm -rf "$work_dir"
}
trap cleanup EXIT

mkdir -p "$pkg_dir"
export LC_ALL=C
export LANG=C

echo "Running native filter..."
# Use one native thread so output order matches the single-threaded WASM path.
cargo run --release -- filter -t 1 "$index_path" "$reads_path" -o "$native_out"

echo "Building Node-target WASM bundle..."
export CARGO_BUILD_CONFIG="$config_path"
wasm-pack build "$repo_root/deacon-wasm" --target nodejs --no-opt --out-dir "$pkg_dir"

echo "Running WASM filter..."
node "$repo_root/scripts/wasm_filter_to_file.mjs" \
  --pkg "$pkg_dir" \
  --index "$index_path" \
  --reads "$reads_path" \
  --output "$wasm_out"

if ! cmp -s "$native_out" "$wasm_out"; then
  echo "Native/WASM output mismatch" >&2
  echo "  native: $(shasum -a 256 "$native_out" | awk '{print $1}')" >&2
  echo "  wasm:   $(shasum -a 256 "$wasm_out" | awk '{print $1}')" >&2
  exit 1
fi

sha="$(shasum -a 256 "$native_out" | awk '{print $1}')"
size="$(wc -c < "$native_out" | tr -d ' ')"

echo "Native/WASM parity check passed"
echo "  bytes: $size"
echo "  sha256: $sha"

# Second pass: --rename + --fasta (-R -f) must also match byte-for-byte.
native_rf="$work_dir/native.rename-fasta.fasta"
wasm_rf="$work_dir/wasm.rename-fasta.fasta"
echo "Running native filter (-R -f)..."
cargo run --release -- filter -t 1 -R -f "$index_path" "$reads_path" -o "$native_rf"
echo "Running WASM filter (--rename --fasta)..."
node "$repo_root/scripts/wasm_filter_to_file.mjs" \
  --pkg "$pkg_dir" \
  --index "$index_path" \
  --reads "$reads_path" \
  --output "$wasm_rf" \
  --rename --fasta

if ! cmp -s "$native_rf" "$wasm_rf"; then
  echo "Native/WASM output mismatch (rename+fasta)" >&2
  echo "  native: $(shasum -a 256 "$native_rf" | awk '{print $1}')" >&2
  echo "  wasm:   $(shasum -a 256 "$wasm_rf" | awk '{print $1}')" >&2
  exit 1
fi
echo "Native/WASM rename+fasta parity check passed ($(wc -c < "$native_rf" | tr -d ' ') bytes)"
