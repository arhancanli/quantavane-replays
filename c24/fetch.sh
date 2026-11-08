#!/bin/sh
# Download the LaTeX source of Akbas and Sra, arXiv 2609.27172v2 (the e-print is a gzipped tar), into its own new
# directory; copy out komlos-paper.tex alone (tarfile, member checked); check SHA-256 against SHA256SUMS.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d as2026 ] && [ -n "$(ls -A as2026)" ]; then echo "refusing: $root/as2026 is not empty"; exit 1; fi
mkdir -p as2026
curl -s -f -L -o as2026/eprint.tar.gz "https://arxiv.org/e-print/2609.27172v2"
python3 -I -c "import tarfile; t = tarfile.open('as2026/eprint.tar.gz'); m = t.getmember('komlos-paper.tex'); assert m.isfile(); open('as2026/komlos-paper.tex', 'wb').write(t.extractfile(m).read())"
shasum -a 256 -c "$here/SHA256SUMS"
