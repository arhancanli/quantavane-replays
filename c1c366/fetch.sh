#!/bin/sh
# Download the C_1c certificate of PR #366 (data only, no code) into a new directory and check its SHA-256.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pr366 ] && [ -n "$(ls -A pr366)" ]; then echo "refusing: $root/pr366 is not empty"; exit 1; fi
mkdir -p pr366
# PR #366, C_1c >= 0.4103735 (intervalplace/frontier @ c0bcdfb, c1c-barnard-steinerberger; SHA-256 stated in the PR)
curl -s -f -L -o pr366/K98304_0.4103735.json \
  https://raw.githubusercontent.com/intervalplace/frontier/c0bcdfb7c36fa20db3d3e32ac13963d790e6a33c/c1c-barnard-steinerberger/K98304_0.4103735.json
(cd pr366 && shasum -a 256 -c "$here/SHA256SUMS")
