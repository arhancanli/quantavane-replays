#!/bin/sh
# Download the integer-coordinate witnesses of Rohan Arun's Erdos problem 838 repository (data only, no code) into a new
# directory and check the SHA-256: rohanarun/erdos-838-convex-polygon-minima @ c301008: reproduction/witness10.txt,
# witness11.txt, witness12.txt, reproduction/real13.jsonl (five 13-point realizations), lean-sources/Witness13.lean
# (read only for its coordinate list), reproduction/fourteen/finite-fourteen/witness780-small.txt and witness780.txt.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d erdos838 ] && [ -n "$(ls -A erdos838)" ]; then echo "refusing: $root/erdos838 is not empty"; exit 1; fi
mkdir -p erdos838
base=https://raw.githubusercontent.com/rohanarun/erdos-838-convex-polygon-minima/c3010080162c8424880222f6b0f0ff7aa5800f9f
for f in reproduction/witness10.txt reproduction/witness11.txt reproduction/witness12.txt reproduction/real13.jsonl \
         lean-sources/Witness13.lean reproduction/fourteen/finite-fourteen/witness780-small.txt \
         reproduction/fourteen/finite-fourteen/witness780.txt; do
  curl -s -f -L -o "erdos838/$(basename "$f")" "$base/$f"
done
shasum -a 256 -c "$here/SHA256SUMS"
