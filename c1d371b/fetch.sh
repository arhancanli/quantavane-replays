#!/bin/sh
# Download the eighth 1d step function of PR #371 (data only, no code) into a new directory and check its SHA-256:
# intervalplace/frontier @ 1ab20c7, c1d-signed-autocorrelation/C1d_327680_1.4546676818.txt (as stated in #371).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pr371b ] && [ -n "$(ls -A pr371b)" ]; then echo "refusing: $root/pr371b is not empty"; exit 1; fi
mkdir -p pr371b
curl -s -f -L -o pr371b/C1d_327680_1.4546676818.txt \
  https://raw.githubusercontent.com/intervalplace/frontier/1ab20c7ca4376cb7824b6d9fcc8408e620bd716d/c1d-signed-autocorrelation/C1d_327680_1.4546676818.txt
shasum -a 256 -c "$here/SHA256SUMS"
