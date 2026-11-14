#!/bin/sh
# Download the covering certificate of the 13a upper bound [K2026] (data only, no code) into a new directory and
# check every SHA-256: ethan-keller/moser-worm-improved-bounds @ 9a33941 (the commit the paper and PR #314 cite),
# certificates/upper/: header.json, tree.json.gz, the 50 multiplier files and manifest.json (about 23 MB). The
# hashes in SHA256SUMS are those of the certificate's own manifest.json (plus the manifest's own hash).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d k2026up ] && [ -n "$(ls -A k2026up)" ]; then echo "refusing: $root/k2026up is not empty"; exit 1; fi
mkdir -p k2026up
base=https://raw.githubusercontent.com/ethan-keller/moser-worm-improved-bounds/9a3394130fe4f7b1df6ffdeffda16601795d9ce4/certificates/upper
for f in $(sed 's#.*k2026up/##' "$here/SHA256SUMS"); do
  curl -s -f -L -o "k2026up/$f" "$base/$f"
done
shasum -a 256 -c "$here/SHA256SUMS" | grep -v ': OK$' || true
shasum -a 256 -c "$here/SHA256SUMS" > /dev/null && echo "all 53 files: OK"
