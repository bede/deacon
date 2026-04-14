#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/deacon-wasm-smoke.XXXXXX")"
pkg_dir="$work_dir/pkg"
config_path="$repo_root/.cargo/config-wasm.toml"

cleanup() {
  rm -rf "$work_dir"
}
trap cleanup EXIT

mkdir -p "$pkg_dir"

export CARGO_BUILD_CONFIG="$config_path"

wasm-pack build "$repo_root/deacon-wasm" --target nodejs --no-opt --out-dir "$pkg_dir"

node "$repo_root/scripts/wasm_smoke_test.mjs" \
  --pkg "$pkg_dir" \
  --index "$repo_root/data/mn908947.idx" \
  --reads "$repo_root/data/mn908947.fastq"
