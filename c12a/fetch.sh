#!/bin/sh
# Download the data replayed here (no code): RECIPE.md and the nine spec files and nine values files of
# devinokeefe/travelling-sales-person-constant-bounds at 69b864b9 (main on 6 Nov 2026, the repository of PR #222),
# into a new directory, and check their SHA-256 against SHA256SUMS (the values in the repository's SHA256SUMS.txt).
# The repository's checkers, controls and workflow are neither downloaded nor run.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d tsp ] && [ -n "$(ls -A tsp)" ]; then echo "refusing: $root/tsp is not empty"; exit 1; fi
mkdir -p tsp/data/specs tsp/data/values
url=https://raw.githubusercontent.com/devinokeefe/travelling-sales-person-constant-bounds/69b864b97262f2c717a3bd3393a60577404fdfae
curl -s -f -L -o tsp/RECIPE.md "$url/RECIPE.md"
for s in tiny rounding bhh k4 k5 k6 k7 k8 k9; do
  curl -s -f -L -o "tsp/data/specs/$s.txt" "$url/data/specs/$s.txt"
  curl -s -f -L -o "tsp/data/values/$s.txt" "$url/data/values/$s.txt"
done
shasum -a 256 -c "$here/SHA256SUMS"
