#!/bin/sh
# Download the WT-triple witnesses of Carter, O'Bryant and Tao (data only, no code) and check their SHA-256:
# dcartermath/sidon-error @ ff5b032, folders asymptotic/ (upsilon = 10) and explicit/ (upsilon = 3).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in asymptotic explicit; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
  for f in G.txt H.txt p.txt upsilon.txt; do
    curl -s -f -L -o $d/$f https://raw.githubusercontent.com/dcartermath/sidon-error/ff5b032fc3745140549777a96a9ec3ba7f8dd954/$d/$f
  done
done
shasum -a 256 -c "$here/SHA256SUMS"
