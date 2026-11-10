#!/bin/sh
# Download the statement of the 41b bound C41b <= 2 sqrt 2 - 1* (#358, vidigoat): the README of tao/certificates/c41b at
# the commit the pull request links, into a new directory, and check its SHA-256. sofa_check.py does not read it: the
# check is written from that README's proof; the file is pinned so the text we worked from can be compared. The folder's
# check_sofa.py is not downloaded.  Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
curl -s -f -L -o c41b_README.md https://raw.githubusercontent.com/vidigoat/optimization-records/cf620ea3fe673782b4e8b6d414012437d6165590/tao/certificates/c41b/README.md
shasum -a 256 -c "$here/SHA256SUMS"
