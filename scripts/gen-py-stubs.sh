#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# macOS needs -undefined dynamic_lookup to link the cdylib outside maturin
if [[ "$(uname)" == "Darwin" ]]; then
    export RUSTFLAGS="-C link-arg=-undefined -C link-arg=dynamic_lookup"
    LIB=target/debug/libdeacon_py.dylib
else
    LIB=target/debug/libdeacon_py.so
fi

cargo build -p deacon-py --features python-stubs
cargo run -q -p deacon-py --features python-stubs --bin gen_stubs -- "$LIB"
