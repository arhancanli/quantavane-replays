#!/bin/sh
# Download the LaTeX source of Ma and Tang, arXiv 2602.23282v2 (757_submit_v2.tex from the e-print: the 14-point
# (4,5)-set A_base of its Lemma 6.1 and the Sidon subset its proof prints), into a new directory, and check its
# SHA-256. mt26_check.py reads it only as text.  Usage: ./fetch_mt26.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d mt26 ] && [ -n "$(ls -A mt26)" ]; then echo "refusing: $root/mt26 is not empty"; exit 1; fi
mkdir -p mt26
curl -s -f -L -o mt26/eprint.tar.gz "https://arxiv.org/e-print/2602.23282v2"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("757_submit_v2.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' mt26/eprint.tar.gz mt26/757_submit_v2.tex
rm mt26/eprint.tar.gz
shasum -a 256 -c "$here/SHA256SUMS_mt26"
