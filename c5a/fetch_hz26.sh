#!/bin/sh
# Download, each into its own new directory, and check the SHA-256 of:
#   hz26/main830.tex   the LaTeX source of Hou and Zhao, arXiv 2607.01169v3 (Theorem 2.1, the definitions of a, b and
#                      the covering inequalities, section 5: R, m, L, the eight mixing weights, eta, the stated a, b and
#                      least slack);
#   hz26/sidon_certificate_8kernel.py   HbZhao1/sidon-vector-smoothing @ ef04456: the certificate file the paper
#                      identifies by its SHA-256 (957a5afa...92a8, as the paper states). hz26_check.py reads only its
#                      integer data literals with Python's ast module (LAMBDA_INTS, P*_HALF_INTS, W*_INTS); nothing
#                      in it is run, and its checking code is not used.
# Usage: ./fetch_hz26.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d hz26 ] && [ -n "$(ls -A hz26)" ]; then echo "refusing: $root/hz26 is not empty"; exit 1; fi
mkdir -p hz26
curl -s -f -L -o hz26/eprint.tar.gz "https://arxiv.org/e-print/2607.01169v3"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("main830.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' hz26/eprint.tar.gz hz26/main830.tex
rm hz26/eprint.tar.gz
curl -s -f -L -o hz26/sidon_certificate_8kernel.py "https://raw.githubusercontent.com/HbZhao1/sidon-vector-smoothing/ef044564300e546f8832b31f5fba133fd192cc3a/sidon_certificate_8kernel.py"
shasum -a 256 -c "$here/SHA256SUMS_hz26"
