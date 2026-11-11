#!/bin/sh
# Download the LaTeX sources of Bondarenko, arXiv 1305.2584v3 (borsuk.tex: Theorem 1, the configuration Gamma_f with
# p = 1/5, q = -1/15) and Jenrich (with Brouwer), arXiv 1308.0206v6 (BORSUK64.tex only: the construction of G_2(4) with
# the Hermitian form a1 b3^4 + a2 b2^4 + a3 b1^4, the sets B_1, B_2, B_3, C, claim (1) and the vector p), each into its
# own new directory, and check their SHA-256. bon_jb_check.py reads them only as text; the Pascal program G24CHK.PAS
# in the second e-print is not extracted.  Usage: ./fetch_bonjb.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in bon2014 jb2014; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o bon2014/eprint.gz "https://arxiv.org/e-print/1305.2584v3"
python3 -I -c 'import gzip, sys
open(sys.argv[2], "wb").write(gzip.open(sys.argv[1]).read())' bon2014/eprint.gz bon2014/borsuk.tex
rm bon2014/eprint.gz
curl -s -f -L -o jb2014/eprint.tar.gz "https://arxiv.org/e-print/1308.0206v6"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("BORSUK64.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' jb2014/eprint.tar.gz jb2014/BORSUK64.tex
rm jb2014/eprint.tar.gz
shasum -a 256 -c "$here/SHA256SUMS_bonjb"
