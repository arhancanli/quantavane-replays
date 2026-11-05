#!/bin/sh
# Download the one data file compared here (no code): [Syl2026]'s k = 2 cost histogram, at the commit the 3a page pins,
# into a new directory, and check its SHA-256 (the theorem note states it). Everything else is rebuilt from the note.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d syl2026 ] && [ -n "$(ls -A syl2026)" ]; then echo "refusing: $root/syl2026 is not empty"; exit 1; fi
mkdir -p syl2026
curl -s -f -L -o syl2026/k2_cost_hist.tsv \
  "https://raw.githubusercontent.com/xSyllaa/c3a-controlled-carry-k2/5ca7cb88820995bf7666096773fd3799f5714b18/data/k2_cost_hist.tsv"
shasum -a 256 -c "$here/SHA256SUMS"
