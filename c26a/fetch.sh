#!/bin/sh
# Download the data of #284 (26a) at the commit it pins into a new directory and check their SHA-256 (the six instance
# hashes are those WRITEUP.md states): WRITEUP.md and problems/bh_cube/FORMAT.md (read, not run) and the six instances.
# cycle.py, core.py and verify/ are not downloaded.  Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
B=https://raw.githubusercontent.com/DennisKoshta/c26a-bohnenblust-hille-cube-family/2557cdb23d08cd45414e4dc312b5a81d629e0b20
curl -s -f -L -o WRITEUP.md $B/WRITEUP.md
curl -s -f -L -o FORMAT.md $B/problems/bh_cube/FORMAT.md
for f in d2_n5 d3_n8 d4_n11 cycle_d5_n15 cycle_d6_n18 cycle_d7_n21; do
  curl -s -f -L -o $f.json $B/problems/bh_cube/instances/$f.json
done
shasum -a 256 -c "$here/SHA256SUMS"
