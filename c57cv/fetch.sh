#!/bin/sh
# Download the data of the 57c row [vS2026] (W. van Staden, U <= 0.65628589*) at the repository's commit a6f17b5 into a
# new directory, and check their SHA-256 against SHA256SUMS: the paper's LaTeX source (read, not run), the inradius cells
# (certify/inradius_leaves_G.npz) and the density coefficients (certify/density_G_n64.json). No program of the
# repository is downloaded.  Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
B=https://raw.githubusercontent.com/WynandVStaden/bloch-landau-slit-discs/a6f17b56e107522655e94eee4138621d3d7c6979
curl -s -f -L -o bloch_slit_discs.tex $B/paper/bloch_slit_discs.tex
curl -s -f -L -o inradius_leaves_G.npz $B/certify/inradius_leaves_G.npz
curl -s -f -L -o density_G_n64.json $B/certify/density_G_n64.json
shasum -a 256 -c "$here/SHA256SUMS"
