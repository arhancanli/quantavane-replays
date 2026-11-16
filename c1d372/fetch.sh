#!/bin/sh
# Download the 1d step function of PR #372 (data only, no code) into a new directory and check its SHA-256:
# ArthurAnache/Hackathon-AI-Reasearch @ c0f8d18, records/1d/autonomous/c006_31bbf628_int56.json (0.8 MB).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pr372 ] && [ -n "$(ls -A pr372)" ]; then echo "refusing: $root/pr372 is not empty"; exit 1; fi
mkdir -p pr372
curl -s -f -L -o pr372/c006_31bbf628_int56.json \
  https://raw.githubusercontent.com/ArthurAnache/Hackathon-AI-Reasearch/c0f8d1837e65d9194350edb7c029897cf031908d/records/1d/autonomous/c006_31bbf628_int56.json
shasum -a 256 -c "$here/SHA256SUMS"
