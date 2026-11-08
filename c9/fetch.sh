#!/bin/sh
# Download the independent set of [IRCR2026] (134,753 vertices of the strong power C_7^10; nathanielitty/lower-bounds-
# for-shannon-capacity at commit fe89660, c7/CC_10_7_134753.txt; the page links this repository) into its own new
# directory and check its SHA-256. Only the data file is fetched: the repository's verifier.py and construction
# scripts are neither read nor run. Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d ircr2026 ] && [ -n "$(ls -A ircr2026)" ]; then echo "refusing: $root/ircr2026 is not empty"; exit 1; fi
mkdir -p ircr2026
curl -s -f -L -o ircr2026/CC_10_7_134753.txt \
  "https://raw.githubusercontent.com/nathanielitty/lower-bounds-for-shannon-capacity/fe8966079013abaac21d1c08441500100f8e3b86/c7/CC_10_7_134753.txt"
shasum -a 256 -c "$here/SHA256SUMS"
