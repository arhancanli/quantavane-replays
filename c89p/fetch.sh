#!/bin/sh
# Download, each into its own new directory, and check the SHA-256 of:
#   ppss/graphs.zip, ppss/c4.graph6.txt   the graph data of Parczyk, Pokutta, Spiegel and Szabo (Zenodo record 6602512,
#              the current version of 10.5281/zenodo.6364588, which the paper cites): only graphs/c4.graph6.txt, the
#              768-vertex graph of their Theorem 1.1, is extracted (with Python's zipfile, as data);
#   ppss/main.tex   the LaTeX source of arXiv 2206.04036v3 (the statement c_4 <= 4551721 * 2^-24 * 3^-2 < 0.030145);
#   ezl/ind.tex     the LaTeX source of Even-Zohar and Linial, arXiv 1312.1205v2 (section 2: Thomason's product and the
#              nested blow-up, with the values 11411/373248, 3769/124416 and 1411/46592).
# The certificates.zip of the same record (flag-algebra lower bounds) is not downloaded.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in ppss ezl; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o ppss/graphs.zip "https://zenodo.org/api/records/6602512/files/graphs.zip/content"
python3 -I -c 'import sys, zipfile
z = zipfile.ZipFile(sys.argv[1]); open(sys.argv[2], "wb").write(z.read("graphs/c4.graph6.txt"))' ppss/graphs.zip ppss/c4.graph6.txt
curl -s -f -L -o ppss/eprint.gz "https://arxiv.org/e-print/2206.04036v3"
python3 -I -c 'import gzip, sys
open(sys.argv[2], "wb").write(gzip.open(sys.argv[1]).read())' ppss/eprint.gz ppss/main.tex
rm ppss/eprint.gz
curl -s -f -L -o ezl/eprint.tar.gz "https://arxiv.org/e-print/1312.1205v2"
python3 -I -c 'import sys, tarfile
t = tarfile.open(sys.argv[1]); m = t.getmember("ind.tex")
open(sys.argv[2], "wb").write(t.extractfile(m).read())' ezl/eprint.tar.gz ezl/ind.tex
rm ezl/eprint.tar.gz
shasum -a 256 -c "$here/SHA256SUMS"
