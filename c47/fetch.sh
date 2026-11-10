#!/bin/sh
# Download the C47 measures replayed here (data files only, no code), each into a new directory, at the commits the
# pull requests pin, and check their SHA-256 against SHA256SUMS.
# - [Tom2026]: lowerbound/c2lower_1p70272.json of tejasvi/square-maximal-weak-type @ 3af09fc (tag v5). Zenodo
#   10.5281/zenodo.23716078 archives the same commit as square-maximal-weak-type-3af09fc359ce.zip (sha256 a91bbc8e...);
#   the copy inside the zip has the same SHA-256 as this file.
# - PR #321: half59-measure.json and the certificate half59.json of rpinho/math-certificates @ 3d4cd55 (c47/round1).
# - PR #341: construction.json is only in Zenodo 10.5281/zenodo.24022791, which is restricted (files: HTTP 403). If you
#   have the file, put it in CERTS/pr341/construction.json; run_all.sh checks its SHA-256 (5f7a8d9d...310b, as the PR
#   states) before using it.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # get DIR FILE REPO COMMIT PATH: one file into a directory that must be new or empty
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "https://raw.githubusercontent.com/$3/$4/$5"
}
# [Tom2026], 1.7027237335 (Theorem 1.3): the 50-atom measure
get tom2026 c2lower_1p70272.json tejasvi/square-maximal-weak-type 3af09fc359ce262f4a7f63d34b0d97e6d4d95d1a \
    lowerbound/c2lower_1p70272.json
# PR #321, 1.70311634* (59 atoms): the measure in Tomar's input format, then the 559-witness certificate
get pr321 half59-measure.json rpinho/math-certificates 3d4cd5574f894d7df2159dfd29d9dd99461cbe1f \
    c47/round1/certificates/half59-measure.json
curl -s -f -L -o pr321/half59.json \
    "https://raw.githubusercontent.com/rpinho/math-certificates/3d4cd5574f894d7df2159dfd29d9dd99461cbe1f/c47/round1/certificates/half59.json"
# PR #341, 1.7032090736* (63 atoms): try the restricted record; a refusal is reported, not an error
if [ ! -f pr341/construction.json ]; then
  mkdir -p pr341
  if curl -s -f -L -o pr341/construction.json \
      "https://zenodo.org/api/records/24022791/files/construction.json/content" 2>/dev/null; then
    echo "pr341: construction.json downloaded"
  else
    rm -f pr341/construction.json
    echo "pr341: Zenodo 10.5281/zenodo.24022791 is restricted; construction.json was not downloaded"
  fi
fi
# PR #350 (rpinho): the round-2 certificate (measure and witness rectangles, one file) at commit 7273e97
if [ ! -f pr350/final.json ]; then
  mkdir -p pr350
  curl -s -f -L -o pr350/final.json \
    "https://raw.githubusercontent.com/rpinho/math-certificates/7273e97ad18bad296571ce55b1378459fe10fc19/c47/round2/certificates/final.json"
fi
# PR #350 (rpinho), update of 10 Nov (round 8, 440 atoms and 7806 witnesses) at commit 85c4680
if [ ! -f pr350r8/final.json ]; then
  mkdir -p pr350r8
  curl -s -f -L -o pr350r8/final.json \
    "https://raw.githubusercontent.com/rpinho/math-certificates/85c4680cc57130d4e9f3f0f48643e51775cdb286/c47/round8/certificates/final.json"
fi
if [ -f pr341/construction.json ]; then
  shasum -a 256 -c "$here/SHA256SUMS"
else
  grep -v pr341/ "$here/SHA256SUMS" | shasum -a 256 -c -
fi
