#!/bin/sh
# Download the LaTeX source of Kunisky, arXiv 2111.02974v1 (unsatisfiable-matrices.tex from the e-print: Theorem 1.1 and
# the recursion of its second proof), into a new directory, and check its SHA-256. kunisky_check.py reads it only as
# text (to confirm the statement and the recursion it implements).  Usage: ./fetch.sh [CERTS_ROOT]  (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d kun2023 ] && [ -n "$(ls -A kun2023)" ]; then echo "refusing: $root/kun2023 is not empty"; exit 1; fi
mkdir -p kun2023
curl -s -f -L -o kun2023/eprint.tar.gz "https://arxiv.org/e-print/2111.02974v1"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("unsatisfiable-matrices.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' kun2023/eprint.tar.gz kun2023/unsatisfiable-matrices.tex
rm kun2023/eprint.tar.gz
shasum -a 256 -c "$here/SHA256SUMS"
