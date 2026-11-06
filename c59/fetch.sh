#!/bin/sh
# Download the one data file replayed here (no code): the certificate of PR #297 (gomyway/c59-bohr-bidisc-certificate,
# the commit the pull request links), into a new directory, and check its SHA-256 against the value the PR and the
# repository's README print. Its checker t7428_c59.py is neither downloaded nor run.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d gomyway ] && [ -n "$(ls -A gomyway)" ]; then echo "refusing: $root/gomyway is not empty"; exit 1; fi
mkdir -p gomyway
curl -s -f -L -o gomyway/c59_cert.json \
  "https://raw.githubusercontent.com/gomyway/c59-bohr-bidisc-certificate/a014dc9455815c2875507017ebcd303b8da27b63/c59_cert.json"
shasum -a 256 -c "$here/SHA256SUMS"
