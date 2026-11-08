#!/bin/sh
# Download Althoefer's note (the [Al26] row of 84b; https://althofer.de/improved_constant_052.tex, the TeX source the
# page links) into its own new directory and check its SHA-256 (the file as fetched on 8 November 2026).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d al26 ] && [ -n "$(ls -A al26)" ]; then echo "refusing: $root/al26 is not empty"; exit 1; fi
mkdir -p al26
curl -s -f -L -o al26/improved_constant_052.tex "https://althofer.de/improved_constant_052.tex"
shasum -a 256 -c "$here/SHA256SUMS"
