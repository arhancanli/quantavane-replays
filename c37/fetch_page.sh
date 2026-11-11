#!/bin/sh
# Download the 37a page of teorth/optimizationproblems at main 87b4ebb (the Kushilevitz and not-all-equal polynomials
# are printed in its comments) into its own new directory, and check its SHA-256. kushilevitz_check.py reads it only
# as data.
# Usage: ./fetch_page.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d page ] && [ -n "$(ls -A page)" ]; then echo "refusing: $root/page is not empty"; exit 1; fi
mkdir -p page
curl -s -f -L -o page/37a.md "https://raw.githubusercontent.com/teorth/optimizationproblems/87b4ebbd426ff8db597a2004cc45da17903bdcf2/constants/37a.md"
shasum -a 256 -c "$here/SHA256SUMS_page"
