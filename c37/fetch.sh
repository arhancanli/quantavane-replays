#!/bin/sh
# Download the three certificate files of PR #283 (data only; the packet's checker, generator, verifier and tools are
# not fetched) at the commit the PR links, into a new directory, and check their SHA-256 against SHA256SUMS.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d kosh ] && [ -n "$(ls -A kosh)" ]; then echo "refusing: $root/kosh is not empty"; exit 1; fi
mkdir -p kosh
B=https://raw.githubusercontent.com/DennisKoshta/c37-degree-sensitivity-certificate/472e2d95c85fe681826e70c90bc60a05ab5c0bc3/certificates
for f in d10_4.cert.xz d7_3.cert d4_2.cert; do curl -s -f -L -o kosh/$f $B/$f; done
shasum -a 256 -c "$here/SHA256SUMS"
