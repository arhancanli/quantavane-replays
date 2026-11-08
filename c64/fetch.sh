#!/bin/sh
# Download the LaTeX source of Li and Yang, arXiv 2308.14859v2 (the e-print is one gzipped main.tex), into its own new
# directory, and check its SHA-256 against SHA256SUMS. Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d ly2023 ] && [ -n "$(ls -A ly2023)" ]; then echo "refusing: $root/ly2023 is not empty"; exit 1; fi
mkdir -p ly2023
curl -s -f -L -o ly2023/main.tex.gz "https://arxiv.org/e-print/2308.14859v2"
gunzip ly2023/main.tex.gz
shasum -a 256 -c "$here/SHA256SUMS"
