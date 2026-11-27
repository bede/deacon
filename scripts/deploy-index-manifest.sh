#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
exec oci os object put -bn human-genome-bucket --file bucket/manifest.json \
  --name deacon/3/manifest.json --content-type application/json \
  --cache-control 'public, max-age=60' "$@"
