#!/bin/sh
# Download the one data file replayed here (no code): the certificate of PR #337 (ArseniyVarlamov/dual-alpha-witness,
# the commit the pull request links), into a new directory, and check its SHA-256. README.md there is read, not run.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d var2026 ] && [ -n "$(ls -A var2026)" ]; then echo "refusing: $root/var2026 is not empty"; exit 1; fi
mkdir -p var2026
curl -s -f -L -o var2026/certificate.json \
  "https://raw.githubusercontent.com/ArseniyVarlamov/dual-alpha-witness/d2868faf31d83e27198868125e2b739637af7004/certificate.json"
shasum -a 256 -c "$here/SHA256SUMS"
