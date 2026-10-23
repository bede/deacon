#!/usr/bin/env bash
# Rebuilds the browser bundle first (--no-build to deploy the existing pkg/).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO/deacon-wasm"
PAGES_ROOT="${PAGES_ROOT:-/Users/bede/Research/git/bede.github.io}"
APP_DIR="deacon"
DEST="$PAGES_ROOT/$APP_DIR"

[[ "${1:-}" == "--no-build" ]] || "$REPO/scripts/build-wasm.sh"

SOURCE_HASH="$(git -C "$REPO" rev-parse --short HEAD)"
# We copy the working tree, not HEAD, so don't label uncommitted code with a clean hash
if [[ -n "$(git -C "$REPO" status --porcelain -- deacon-wasm)" ]]; then
  SOURCE_HASH="$SOURCE_HASH-dirty-$(date -u +%Y%m%d%H%M%S)"
fi

# Deployed files (relative to both SRC and DEST); tracked files outside this set are pruned
FILES=(
  index.html
  app.js
  worker.js
  protocol.js
  pairing.js
  pkg/deacon_wasm.js
  pkg/deacon_wasm_bg.wasm
)

mkdir -p "$DEST/pkg"
for f in "${FILES[@]}"; do
  cp "$SRC/$f" "$DEST/$f"
done

# Stamp every cache token with the source revision so one deployment is one module graph
perl -i -pe "s/(ASSET_VERSION = \"|[?]v=)\d{8}-\d+/\${1}$SOURCE_HASH/g" \
  "$DEST/index.html" "$DEST/app.js" "$DEST/worker.js"

echo "Deployed deacon-wasm $SOURCE_HASH to $DEST"

cd "$PAGES_ROOT"
while IFS= read -r tracked; do
  rel="${tracked#"$APP_DIR"/}"
  keep=0
  for f in "${FILES[@]}"; do [[ "$f" == "$rel" ]] && { keep=1; break; }; done
  [[ "$keep" == 0 ]] && git rm -q "$tracked"
done < <(git ls-files "$APP_DIR/")
for f in "${FILES[@]}"; do
  git add "$APP_DIR/$f"
done

if git diff --cached --quiet; then
  echo "No GitHub Pages changes to commit"
  exit 0
fi
git commit -m "deacon-wasm $SOURCE_HASH"
git push origin "$(git rev-parse --abbrev-ref HEAD)"
echo "Pushed to GitHub Pages"
