#!/bin/sh
# Download the LaTeX source of Norton, arXiv 2608.20292v3 (one paper.tex), into its own new directory, and check its
# SHA-256 against SHA256SUMS. Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d n2026 ] && [ -n "$(ls -A n2026)" ]; then echo "refusing: $root/n2026 is not empty"; exit 1; fi
mkdir -p n2026
curl -s -f -L -o n2026/eprint.tar.gz "https://arxiv.org/e-print/2608.20292v3"
(cd n2026 && tar xzf eprint.tar.gz && rm eprint.tar.gz)
shasum -a 256 -c "$here/SHA256SUMS"
