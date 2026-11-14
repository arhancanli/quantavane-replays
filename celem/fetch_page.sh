#!/bin/sh
# Download the four constant pages whose rows are replayed here (teorth/optimizationproblems at main 87b4ebb, the
# commit of 27 October 2026 that is still main on 14 November) and check their SHA-256. The programs do not read the
# pages: each object (the cube, the simplex, f(x) = x) is written into the programs from the row's text, which
# RESULT.md quotes; the pages are pinned so that the rows can be compared with what was replayed.
# Usage: ./fetch_page.sh [ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pages ] && [ -n "$(ls -A pages)" ]; then echo "refusing: $root/pages is not empty"; exit 1; fi
mkdir -p pages
for p in 20b 25a 39a 11a; do
  curl -s -f -L -o pages/$p.md \
    https://raw.githubusercontent.com/teorth/optimizationproblems/87b4ebbd426ff8db597a2004cc45da17903bdcf2/constants/$p.md
done
(cd pages && shasum -a 256 -c "$here/SHA256SUMS")
