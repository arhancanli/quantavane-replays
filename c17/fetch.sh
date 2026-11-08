#!/bin/sh
# Download the LaTeX sources this replay reads, each e-print into its own new directory, and check the SHA-256 of
# the files used against SHA256SUMS: Wang et al., arXiv 2603.15617v2 (the appendix with the quintic, A.3), and
# Gupta, Ndiaye, Norin and Wei, arXiv 2407.19026v2 (the current version) and v1 (the version the benchmark and the
# 17a page cite by theorem number). Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
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
shasum -a 256 -c "$here/SHA256SUMS"
