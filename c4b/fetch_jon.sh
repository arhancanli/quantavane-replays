#!/bin/sh
# Download the data files of Jones's square certificate (JD-Jones-ASES/nk-lean @ 55acbf0, certificates/square/:
# certificate.json, data/odd215.txt, data/odd437.txt, data/binary_policy.json) into their own new directory and check
# their SHA-256 (the odd437.txt hash is the one the repository's README states; the others are the files at that
# commit). Data only: no Lean file or program of the repository is fetched.  Usage: ./fetch_jon.sh [CERTS_ROOT]
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d jon2026b ] && [ -n "$(ls -A jon2026b)" ]; then echo "refusing: $root/jon2026b is not empty"; exit 1; fi
mkdir -p jon2026b
base=https://raw.githubusercontent.com/JD-Jones-ASES/nk-lean/55acbf0/certificates/square
curl -s -f -L -o jon2026b/certificate.json "$base/certificate.json"
for f in odd215.txt odd437.txt binary_policy.json; do curl -s -f -L -o jon2026b/$f "$base/data/$f"; done
shasum -a 256 -c "$here/SHA256SUMS_jon"
