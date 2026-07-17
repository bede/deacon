#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ "$(uname)" == "Darwin" ]]; then
    LIB=target/debug/libdeacon_py.dylib
else
    LIB=target/debug/libdeacon_py.so
fi

cargo build -p deacon-py --features python-stubs
cargo run -q -p deacon-py --features python-stubs --bin gen_stubs -- "$LIB"
