#!/bin/sh
# Download the sources this replay reads, each into its own new directory, and check the SHA-256 of the files used
# against SHA256SUMS: Wang et al., arXiv 2603.15617v2 (the appendix with the quintic, A.3); Gupta, Ndiaye, Norin and
# Wei, arXiv 2407.19026v2 (the current version) and v1 (the version the benchmark and the 17a page cite by theorem
# number); and, for the [Gri26] row, the proof note and the certificate's JSON of maaxgrin/ramsey-3792-bound at
# 0b26aa2 (data and prose only; its verifier is not downloaded). Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # get DIR ARXIV_ID: the e-print into a directory that must be new or empty, unpacked there
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/eprint.tar.gz" "https://arxiv.org/e-print/$2"
  (cd "$1" && tar xzf eprint.tar.gz && rm eprint.tar.gz)
}
get wmr2026 2603.15617v2
get gnnw_v2 2407.19026v2
get gnnw_v1 2407.19026v1
if [ -d gri26 ] && [ -n "$(ls -A gri26)" ]; then echo "refusing: $root/gri26 is not empty"; exit 1; fi
mkdir -p gri26
for f in ramsey_3792_certificate.json ramsey_3792_bound.tex; do
  curl -s -f -L -o "gri26/$f" "https://raw.githubusercontent.com/maaxgrin/ramsey-3792-bound/0b26aa202babc82f9f3150cfc32c54992997480e/$f"
done
shasum -a 256 -c "$here/SHA256SUMS"
