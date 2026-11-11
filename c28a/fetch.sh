#!/bin/sh
# Download the proof PDF of maaxgrin/borsuk-63-counterexample @ cdcdbeac (the text our checker was written from) into
# its own new directory and check its SHA-256. borsuk_check.py does not read it: the construction is rebuilt from the
# PDF's description (Brouwer's projective model of G_2(4)); this is for provenance only.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d gri2026 ] && [ -n "$(ls -A gri2026)" ]; then echo "refusing: $root/gri2026 is not empty"; exit 1; fi
mkdir -p gri2026
curl -s -f -L -o gri2026/borsuk_63_counterexample.pdf \
  "https://raw.githubusercontent.com/maaxgrin/borsuk-63-counterexample/cdcdbeac2e692b8641218c70ce9f414522e125e5/borsuk_63_counterexample.pdf"
shasum -a 256 -c "$here/SHA256SUMS"
