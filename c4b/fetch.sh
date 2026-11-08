#!/bin/sh
# Download the LaTeX source of Krachun, arXiv 2608.01325v1 (one main.tex), the README of Jones's fs-lower-bound at
# c2c0687 (the [Jon2026a] note's definitions and data; its verify.py and Lean files are not fetched), and the 4b page
# at teorth/optimizationproblems 87b4ebb, each into its own new directory, and check their SHA-256 against SHA256SUMS.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d kra2026 ] && [ -n "$(ls -A kra2026)" ]; then echo "refusing: $root/kra2026 is not empty"; exit 1; fi
mkdir -p kra2026
curl -s -f -L -o kra2026/eprint.tar.gz "https://arxiv.org/e-print/2608.01325v1"
(cd kra2026 && tar xzf eprint.tar.gz && rm eprint.tar.gz)
for d in jon2026a page; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o jon2026a/README.md "https://raw.githubusercontent.com/JD-Jones-ASES/fs-lower-bound/c2c06875ee4544e069e8da16564ee429d870a13e/README.md"
curl -s -f -L -o page/4b.md "https://raw.githubusercontent.com/teorth/optimizationproblems/87b4ebbd426ff8db597a2004cc45da17903bdcf2/constants/4b.md"
shasum -a 256 -c "$here/SHA256SUMS"
