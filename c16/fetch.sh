#!/bin/sh
# Download the C16 certificate data file replayed here (no code) into its own new directory, at the commit PR #270
# pins, and check its SHA-256 against SHA256SUMS.
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
# PR #270 update of 4 Nov: C16 <= 2.15622552815422 (353 delta-intervals with their mu)
get vid certificate.json vidigoat/optimization-records 5da2dc1f4c5294a3f3bcb8fb1596580dfbaecec6 \
    tao/certificates/c16_pinned/certificate.json
shasum -a 256 -c "$here/SHA256SUMS"
