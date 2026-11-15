#!/bin/sh
# Download the 1d step function of PR #371 (data only, no code) into a new directory and check its SHA-256:
# intervalplace/frontier @ 316b76c, c1d-signed-autocorrelation/C1d_327680_1.4546733294.txt (3.7 MB).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pr371 ] && [ -n "$(ls -A pr371)" ]; then echo "refusing: $root/pr371 is not empty"; exit 1; fi
mkdir -p pr371
curl -s -f -L -o pr371/C1d_327680_1.4546733294.txt \
  https://raw.githubusercontent.com/intervalplace/frontier/316b76ce6a161e4cf4f96a69faef484fb11a30b1/c1d-signed-autocorrelation/C1d_327680_1.4546733294.txt
shasum -a 256 -c "$here/SHA256SUMS"
