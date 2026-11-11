#!/bin/sh
# Download the LaTeX source of Zheng, arXiv 2506.01896v1 (sumdif.tex from the e-print: Theorem 1 and the table of its
# section 3), into a new directory, and check its SHA-256. z2025_check.py reads it only as text (the MATLAB code
# printed in it is not run). The [Ru96] row needs no download.  Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d z2025 ] && [ -n "$(ls -A z2025)" ]; then echo "refusing: $root/z2025 is not empty"; exit 1; fi
mkdir -p z2025
curl -s -f -L -o z2025/eprint.tar.gz "https://arxiv.org/e-print/2506.01896v1"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("sumdif.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' z2025/eprint.tar.gz z2025/sumdif.tex
rm z2025/eprint.tar.gz
shasum -a 256 -c "$here/SHA256SUMS"
