#!/bin/sh
# Download the two mixture certificates of SeanSagait/c1c-dual-bound-limit @ 0b6cbdd (data only, no code) into a new
# directory and check their SHA-256 (as in that repository's SHA256SUMS).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d mix ] && [ -n "$(ls -A mix)" ]; then echo "refusing: $root/mix is not empty"; exit 1; fi
mkdir -p mix
for f in mix_K4096_P2.json mix_K1024_P2.json; do
  curl -s -f -L -o "mix/$f" \
    "https://raw.githubusercontent.com/SeanSagait/c1c-dual-bound-limit/0b6cbddd0566387110ddd1e44a0b9c1f6e3ed765/$f"
done
(cd mix && shasum -a 256 -c "$here/SHA256SUMS")
