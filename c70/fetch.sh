#!/bin/sh
# Download the C70 certificate data file replayed here (no code) into its own new directory, at the commit the 70a
# page and PR #259 pin, and check its SHA-256 against SHA256SUMS.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pin ] && [ -n "$(ls -A pin)" ]; then echo "refusing: $root/pin is not empty"; exit 1; fi
mkdir -p pin/certificates
# [Pin2026], 1.20194336* (70a page; PR #259)
curl -s -f -L -o pin/certificates/c70_candidate.json \
  "https://raw.githubusercontent.com/rpinho/c70-certificate/680e074bb50929734bfd76484d716cd3e85827f0/certificates/c70_candidate.json"
shasum -a 256 -c "$here/SHA256SUMS"
