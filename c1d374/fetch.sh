#!/bin/sh
# Download the 1d step function of PR #374 (data only, no code) into a new directory and check its SHA-256:
# cweill/c1d-upper-bound @ 7889e93, c1d_n102400_1.4513064007.json (2.0 MB).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pr374 ] && [ -n "$(ls -A pr374)" ]; then echo "refusing: $root/pr374 is not empty"; exit 1; fi
mkdir -p pr374
curl -s -f -L -o pr374/c1d_n102400_1.4513064007.json \
  https://raw.githubusercontent.com/cweill/c1d-upper-bound/7889e938c7ceff997bddf643f018c06719b867b2/c1d_n102400_1.4513064007.json
shasum -a 256 -c "$here/SHA256SUMS"
