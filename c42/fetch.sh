#!/bin/sh
# Download the data files replayed here (no code), each into its own new directory, and check their SHA-256:
# Griego's certificate (sebastian-griego/turan-c42-certificate, v1.0.0 = the commit the 42a page pins; against his
# .sha256 file), Pinho's 256-step profile (rpinho/c42-certificate at the commit the page pins; the value the page
# prints) and Koshta's tau = 41/200 certificate (DennisKoshta/c42-turan-power-sum-certificate at commit 8e1797e0; the
# value its note and SHA256SUMS state). Their checkers are neither downloaded nor run. A file already present is kept
# and checked again; a download goes only into a new or empty directory.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # directory, file, url
  [ -f "$1/$2" ] && return 0
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "$3"
}
get griego turan42_certificate.json \
  "https://raw.githubusercontent.com/sebastian-griego/turan-c42-certificate/c8ddce14d9a5e898406d5dc6b8d08bb8a39507c7/certificate/turan42_certificate.json"
get pinho c42_256.json \
  "https://raw.githubusercontent.com/rpinho/c42-certificate/a8a9bc309bd1992413fa07595bd60c5e3fa1b65b/certificates/c42_256.json"
get koshta c42_tau41_200_M64_0.6873358063.json \
  "https://raw.githubusercontent.com/DennisKoshta/c42-turan-power-sum-certificate/8e1797e0d35454ad18cfebbbeb5858bf202b0748/certificates/c42_tau41_200_M64_0.6873358063.json"
shasum -a 256 -c "$here/SHA256SUMS"
