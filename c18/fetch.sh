#!/bin/sh
# Download the LaTeX source of Green and Tao, Freiman's theorem in finite fields via extremal set theory
# (arXiv math/0703668, the e-print is one gzipped main.tex), into its own new directory, and check its SHA-256.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d gt2009 ] && [ -n "$(ls -A gt2009)" ]; then echo "refusing: $root/gt2009 is not empty"; exit 1; fi
mkdir -p gt2009
curl -s -f -L -o gt2009/main.tex.gz "https://arxiv.org/e-print/math/0703668v1"
gunzip gt2009/main.tex.gz
shasum -a 256 -c "$here/SHA256SUMS"
