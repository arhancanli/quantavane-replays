#!/bin/sh
# Download the LaTeX source of Buys, Polak and Zuiddam, arXiv 2607.29681v1 (main3.tex from the e-print), and the base
# data files of their Lean repository (spectra-research/shannon-capacity-lean @ aa21eeb: ShannonBounds/BaseC7Data.lean,
# BaseC11Data.lean, BaseC15Data.lean, BaseC19Data.lean; read only as data) into their own new directory, and check
# their SHA-256.  Usage: ./fetch_bpz.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d bpz2026 ] && [ -n "$(ls -A bpz2026)" ]; then echo "refusing: $root/bpz2026 is not empty"; exit 1; fi
mkdir -p bpz2026
curl -s -f -L -o bpz2026/eprint.tar.gz "https://arxiv.org/e-print/2607.29681v1"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("main3.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' bpz2026/eprint.tar.gz bpz2026/main3.tex
rm bpz2026/eprint.tar.gz
for n in 7 11 15 19; do
  curl -s -f -L -o bpz2026/BaseC${n}Data.lean \
    "https://raw.githubusercontent.com/spectra-research/shannon-capacity-lean/aa21eeb12b75b0413d3fa9fb4208b5d0bf2c4d65/ShannonBounds/BaseC${n}Data.lean"
done
shasum -a 256 -c "$here/SHA256SUMS_bpz"
