#!/bin/sh
# Download the LaTeX sources of Polak and Schrijver, arXiv 1808.07438v2 (c7_367_rev2.tex: the 367 words of its
# appendix), and of Mathew and Ostergard, arXiv 1504.01472v1 (arxiv.tex from the e-print: the generator and orbit
# representatives of its code for G(5,7) >= 350), each into its own new directory, and check their SHA-256.
# ps_mo_check.py reads both only as text.  Usage: ./fetch_psmo.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in ps2018 mo2017; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o ps2018/eprint.gz "https://arxiv.org/e-print/1808.07438v2"
python3 -I -c 'import gzip, sys
open(sys.argv[2], "wb").write(gzip.open(sys.argv[1]).read())' ps2018/eprint.gz ps2018/c7_367_rev2.tex
rm ps2018/eprint.gz
curl -s -f -L -o mo2017/eprint.tar.gz "https://arxiv.org/e-print/1504.01472v1"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("arxiv.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' mo2017/eprint.tar.gz mo2017/arxiv.tex
rm mo2017/eprint.tar.gz
shasum -a 256 -c "$here/SHA256SUMS_psmo"
