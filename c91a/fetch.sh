#!/bin/sh
# Download the one data file this replay uses: [J2026]'s certificate log (the rational x and t for k = 2), at the
# commit PR #273 pins (leandrejack-afk/chvatal-sankoff-bound @ 574220e, edit-distance-lower/); check its SHA-256.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d j2026 ] && [ -n "$(ls -A j2026)" ]; then echo "refusing: $root/j2026 is not empty"; exit 1; fi
mkdir -p j2026
curl -s -f -L -o j2026/certify_rule5_R32.log \
  "https://raw.githubusercontent.com/leandrejack-afk/chvatal-sankoff-bound/574220ee90b332181ba8f8eccf500c380c0f78bb/edit-distance-lower/certify_rule5_R32.log"
cd "$root" && shasum -a 256 -c "$here/SHA256SUMS"
