#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
wasm_dir="$repo_root/deacon-wasm"
pkg_dir="$wasm_dir/pkg"
config_path="$repo_root/.cargo/config-wasm.toml"

mkdir -p "$pkg_dir"

export CARGO_BUILD_CONFIG="$config_path"

cargo check --manifest-path "$wasm_dir/Cargo.toml" --target wasm32-unknown-unknown
wasm-pack build "$wasm_dir" --target web --out-dir "$pkg_dir"

for path in \
  "$pkg_dir/deacon_wasm.js" \
  "$pkg_dir/deacon_wasm_bg.wasm" \
  "$pkg_dir/deacon_wasm.d.ts" \
  "$pkg_dir/package.json"
do
  [[ -f "$path" ]] || {
    echo "Missing generated asset: $path" >&2
    exit 1
  }
done

echo "WASM build check passed"
echo "  pkg: $pkg_dir"
