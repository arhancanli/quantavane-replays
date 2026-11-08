#!/bin/sh
# Download the LaTeX sources of Bai, arXiv 2609.11276v2 (main.tex, sections/, appendix), and of Zeilberger and Zudilin,
# arXiv 1912.06345v2 (one main.tex), each into its own new directory, and check the SHA-256 of the files read against
# SHA256SUMS. Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d b2026 ] && [ -n "$(ls -A b2026)" ]; then echo "refusing: $root/b2026 is not empty"; exit 1; fi
mkdir -p b2026
curl -s -f -L -o b2026/eprint.tar.gz "https://arxiv.org/e-print/2609.11276v2"
(cd b2026 && tar xzf eprint.tar.gz && rm eprint.tar.gz)
if [ -d zz2020 ] && [ -n "$(ls -A zz2020)" ]; then echo "refusing: $root/zz2020 is not empty"; exit 1; fi
mkdir -p zz2020
curl -s -f -L -o zz2020/main.tex.gz "https://arxiv.org/e-print/1912.06345v2"
gunzip zz2020/main.tex.gz
shasum -a 256 -c "$here/SHA256SUMS"
