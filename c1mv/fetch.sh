#!/bin/sh
# Download, each into its own new directory, and check the SHA-256 of:
#   mv2009/sch0903b.tex   the LaTeX source of Matolcsi and Vinuesa, arXiv 0907.1379v2 (from the e-print): the step
#                         functions of its appendix are read as text;
#   ae/dm_mathematical_results.ipynb   google-deepmind/alphaevolve_results @ 4226acb: the heights of section B.3
#                         (the list in cell 65) are read with Python's ast module as literals; nothing is run.
# Then extract_mv.py copies the numbers, as printed, into JSON files for ../autocorr_check.py.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in mv2009 ae; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o mv2009/eprint.tar.gz "https://arxiv.org/e-print/0907.1379v2"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("sch0903b.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' mv2009/eprint.tar.gz mv2009/sch0903b.tex
rm mv2009/eprint.tar.gz
curl -s -f -L -o ae/dm_mathematical_results.ipynb "https://raw.githubusercontent.com/google-deepmind/alphaevolve_results/4226acbf237ff9ad10ba7673a2af127a2d8a5971/mathematical_results.ipynb"
shasum -a 256 -c "$here/SHA256SUMS"
python3 -I "$here/extract_mv.py" "$root"
