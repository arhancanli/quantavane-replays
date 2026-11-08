#!/bin/sh
# Download the [F2026] package from Zenodo (10.5281/zenodo.23714415, open access), check the archive's SHA-256 against
# the one PR #329 states (29c07851...2f158), and extract only the certificate data files: record.json and
# all-groups.json (the replayed bounds), standard.json and standard-five-page.json (self-tests), into a new directory;
# their own SHA-256 are then checked (as the 3b row and the package's CERTIFICATE_INDEX.md state). The package's
# programs are not extracted or run.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ ! -f f2026/record.json ] || [ ! -f f2026/all-groups.json ] || [ ! -f f2026/standard-five-page.json ]; then
  if [ -d f2026 ] && [ -n "$(ls -A f2026)" ]; then echo "refusing: $root/f2026 is not empty"; exit 1; fi
  mkdir -p f2026 zip
  [ -f zip/source-and-certificates.zip ] || curl -s -f -L -o zip/source-and-certificates.zip \
    "https://zenodo.org/api/records/23714415/files/source-and-certificates.zip/content"
  echo "29c0785139d6a6a9f2bbc3c569db9211df1ab1e2da1b29238a52cab5d41f2158  zip/source-and-certificates.zip" | shasum -a 256 -c -
  unzip -q -j zip/source-and-certificates.zip anc/certificates/record.json anc/certificates/all-groups.json \
    anc/certificates/standard.json anc/certificates/standard-five-page.json -d f2026
fi
# [MI2026] (PR #92), the 3b lower bound: Zenodo 10.5281/zenodo.20794135, then only c3b_pr92/certificate_3b_13pt.json
if [ ! -f mi2026/certificate_3b_13pt.json ]; then
  if [ -d mi2026 ] && [ -n "$(ls -A mi2026)" ]; then echo "refusing: $root/mi2026 is not empty"; exit 1; fi
  mkdir -p mi2026
  curl -s -f -L -o mi2026/sumdiff_certs_upload.zip "https://zenodo.org/records/20794135/files/sumdiff_certs_upload.zip?download=1"
  echo "8b0ae7571bbb9e584ed8f6a4e40a9d63031c9aba1526b660d6f56b12dfabfb87  mi2026/sumdiff_certs_upload.zip" | shasum -a 256 -c -
  unzip -p mi2026/sumdiff_certs_upload.zip c3b_pr92/certificate_3b_13pt.json > mi2026/certificate_3b_13pt.json
fi
shasum -a 256 -c "$here/SHA256SUMS"
