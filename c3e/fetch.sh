#!/bin/sh
# Download the two data files replayed here (no code): the C_3e finite certificates in mizopl/ghr-c3a-odd-base-finite-bound
# (c3e_finite_certificate/, at the commit of 6 Nov 2026 06:39 UTC), into a new directory, and check their SHA-256 against
# the values the folder's README prints. They are compared with our results only; the README's Theorems C1-C3 are read,
# and its check_c3e.py is neither downloaded nor run.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d mizopl ] && [ -n "$(ls -A mizopl)" ]; then echo "refusing: $root/mizopl is not empty"; exit 1; fi
mkdir -p mizopl
base="https://raw.githubusercontent.com/mizopl/ghr-c3a-odd-base-finite-bound/2fa6c6a6c708c0b5e729ac4bad085e7d6a040b5d/c3e_finite_certificate"
for f in c3e_finite_certificate_80_161.json c3e_finite_certificate_160_321.json; do
  curl -s -f -L -o "mizopl/$f" "$base/$f"
done
shasum -a 256 -c "$here/SHA256SUMS"
