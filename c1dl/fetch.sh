#!/bin/sh
# Download the statement of the C_1d >= 1.152 bound (#270, vidigoat, 10 Nov 2026, second update): the README of
# tao/certificates/c1d_lower at the commit the pull request links, into a new directory, and check its SHA-256.
# c1d_lower_check.py does not read it: the replay is written from that README's proof; the file is pinned here so the
# text we worked from can be compared. The folder's check_c1d_lower.py and sanity_random.py are not downloaded.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
curl -s -f -L -o c1d_lower_README.md https://raw.githubusercontent.com/vidigoat/optimization-records/a5fcdef1dc0e1ece83736e771491962652b2a1ec/tao/certificates/c1d_lower/README.md
shasum -a 256 -c "$here/SHA256SUMS"
