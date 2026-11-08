#!/bin/sh
# Download the LaTeX source of Keller, arXiv 2609.21968v1 (an e-print with main.tex and sections/), into its own new
# directory, and check the SHA-256 of the two files read against SHA256SUMS. Usage: ./fetch.sh [CERTS_ROOT]
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d k2026 ] && [ -n "$(ls -A k2026)" ]; then echo "refusing: $root/k2026 is not empty"; exit 1; fi
mkdir -p k2026
curl -s -f -L -o k2026/eprint.tar.gz "https://arxiv.org/e-print/2609.21968v1"
(cd k2026 && tar xzf eprint.tar.gz && rm eprint.tar.gz)
shasum -a 256 -c "$here/SHA256SUMS"
