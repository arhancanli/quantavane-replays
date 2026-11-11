#!/bin/sh
# Download the 10c page of teorth/optimizationproblems at main 87b4ebb (the four sign-matrix certificates are printed
# on it) into its own new directory, and check its SHA-256. The page is read only as data by c10c_check.py.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d page ] && [ -n "$(ls -A page)" ]; then echo "refusing: $root/page is not empty"; exit 1; fi
mkdir -p page
curl -s -f -L -o page/10c.md "https://raw.githubusercontent.com/teorth/optimizationproblems/87b4ebbd426ff8db597a2004cc45da17903bdcf2/constants/10c.md"
shasum -a 256 -c "$here/SHA256SUMS"
