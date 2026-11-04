#!/bin/sh
# Download the C44 certificate data files replayed here (no code), each source into its own new directory, at the
# commits the 44a page and PR #258 pin, and check their SHA-256 against SHA256SUMS.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
newdir() {  # a directory that must be new or empty
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
}
get() {  # get DIR FILE REPO COMMIT PATH
  mkdir -p "$(dirname "$1/$2")"
  curl -s -f -L -o "$1/$2" "https://raw.githubusercontent.com/$3/$4/$5"
}
# [P2026], 4.3576549226* (44a page; PR #258): the new small certificate and Epi's data as redistributed
P=rpinho/c44-certificate; PC=f3fb25dfb88de3da8b1a0b30f03fdb51b78ccef9
newdir pin
get pin certificates/c44_candidate.json $P $PC certificates/c44_candidate.json
get pin sources/wellens-certificate.json $P $PC sources/wellens-certificate.json
get pin sources/wellens-result.json $P $PC sources/wellens-result.json
# [Epi2026] (44a page): the same data at BakeLab's pinned commit, for a byte comparison
E=BakeLab/Epi-results; EC=97fa6a420f417458e11723c8bceb70f1529e447a
newdir epi
get epi certificate.json $E $EC mathematics/wellens/certificate.json
get epi result.json $E $EC mathematics/wellens/result.json
shasum -a 256 -c "$here/SHA256SUMS"
