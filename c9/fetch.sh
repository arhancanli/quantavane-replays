#!/bin/sh
# Download the independent set of [IRCR2026] (134,753 vertices of the strong power C_7^10; nathanielitty/lower-bounds-
# for-shannon-capacity at commit fe89660, c7/CC_10_7_134753.txt; the page links this repository) and, for [Gao2026],
# the LaTeX source of arXiv 2607.27869v1 and its 367-point base set (xyz2606/recursive_construction_of_the_Shannon_
# capacity_of_C_7 at b13031b, inputs/R367.txt), each into its own new directory, and check their SHA-256. Only data
# and TeX are fetched: the repositories' verifiers and construction scripts are neither read nor run.
# Also [Tan2026] and PR #296, see below.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d ircr2026 ] && [ -n "$(ls -A ircr2026)" ]; then echo "refusing: $root/ircr2026 is not empty"; exit 1; fi
mkdir -p ircr2026
curl -s -f -L -o ircr2026/CC_10_7_134753.txt \
  "https://raw.githubusercontent.com/nathanielitty/lower-bounds-for-shannon-capacity/fe8966079013abaac21d1c08441500100f8e3b86/c7/CC_10_7_134753.txt"
if [ -d gao2026 ] && [ -n "$(ls -A gao2026)" ]; then echo "refusing: $root/gao2026 is not empty"; exit 1; fi
mkdir -p gao2026
curl -s -f -L -o gao2026/eprint.tar.gz "https://arxiv.org/e-print/2607.27869v1"
(cd gao2026 && tar xzf eprint.tar.gz && rm eprint.tar.gz)
curl -s -f -L -o gao2026/R367.txt \
  "https://raw.githubusercontent.com/xyz2606/recursive_construction_of_the_Shannon_capacity_of_C_7/b13031ba76e373fdff4bd0152ae79d1c26969900/inputs/R367.txt"
# [Tan2026] (arXiv 2608.30273v1): only main.tex and appendix.tex of the e-print are kept (TeX, read, never run); and
# PR #296's certificate.json (gomyway/c9-shannon-c7-lower-bound-certificate at b5d6c2c; data only: its check_c9.py and
# the vendored package are not downloaded)
if [ -d tan2026 ] && [ -n "$(ls -A tan2026)" ]; then echo "refusing: $root/tan2026 is not empty"; exit 1; fi
mkdir -p tan2026
curl -s -f -L -A "Mozilla/5.0" -o tan2026/eprint.tar.gz "https://arxiv.org/e-print/2608.30273v1"
python3 -I -c "
import tarfile
t = tarfile.open('tan2026/eprint.tar.gz')
for n in ('main.tex', 'appendix.tex'):
    m = t.getmember(n)
    assert m.isfile()
    open('tan2026/' + n, 'wb').write(t.extractfile(m).read())
"
rm tan2026/eprint.tar.gz
if [ -d lt2026 ] && [ -n "$(ls -A lt2026)" ]; then echo "refusing: $root/lt2026 is not empty"; exit 1; fi
mkdir -p lt2026
curl -s -f -L -o lt2026/certificate.json \
  "https://raw.githubusercontent.com/gomyway/c9-shannon-c7-lower-bound-certificate/b5d6c2c2c981f89f5c62f3db9d26dab045464f81/certificate.json"
shasum -a 256 -c "$here/SHA256SUMS"
