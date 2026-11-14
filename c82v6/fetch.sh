#!/bin/sh
# Download the C82 certificate of PR #280's update of 14 November 2026 (v6; data only, no code) into a new directory
# and check its SHA-256 (the one stated in the pull request comment).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pr280f ] && [ -n "$(ls -A pr280f)" ]; then echo "refusing: $root/pr280f is not empty"; exit 1; fi
mkdir -p pr280f
# vidigoat/optimization-records @ e3bae2e (14 November 2026), tao/certificates/c82/v6
curl -s -f -L -o pr280f/certificate.json \
  https://raw.githubusercontent.com/vidigoat/optimization-records/e3bae2ea40c8b1f764a1ba08bcd21fe513b2c45d/tao/certificates/c82/v6/certificate.json
(cd pr280f && shasum -a 256 -c "$here/SHA256SUMS")
