#!/bin/sh
# Download the one data file replayed here (no code): the certificate of PR #335 (ArseniyVarlamov/sidon-c5a-certificate,
# the commit the pull request links), into a new directory, and check its SHA-256. PROOF.md is read, not executed.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d var2026 ] && [ -n "$(ls -A var2026)" ]; then echo "refusing: $root/var2026 is not empty"; exit 1; fi
mkdir -p var2026
curl -s -f -L -o var2026/certificate.json \
  "https://raw.githubusercontent.com/ArseniyVarlamov/sidon-c5a-certificate/091844ea70a50bd042d6687ca0004ec7c9f1849f/certificate.json"
shasum -a 256 -c "$here/SHA256SUMS"
