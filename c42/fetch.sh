#!/bin/sh
# Download the two data files replayed here (no code): Griego's certificate (sebastian-griego/turan-c42-certificate,
# v1.0.0 = the commit the 42a page pins) and Pinho's 256-step profile (rpinho/c42-certificate at the commit the page
# pins), into a new directory, and check their SHA-256 (Griego's against his .sha256 file, Pinho's against the value
# the page prints). Their checkers are neither downloaded nor run.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in griego pinho; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o griego/turan42_certificate.json \
  "https://raw.githubusercontent.com/sebastian-griego/turan-c42-certificate/c8ddce14d9a5e898406d5dc6b8d08bb8a39507c7/certificate/turan42_certificate.json"
curl -s -f -L -o pinho/c42_256.json \
  "https://raw.githubusercontent.com/rpinho/c42-certificate/a8a9bc309bd1992413fa07595bd60c5e3fa1b65b/certificates/c42_256.json"
shasum -a 256 -c "$here/SHA256SUMS"
