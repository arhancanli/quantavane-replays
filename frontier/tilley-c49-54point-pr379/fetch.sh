#!/bin/sh
# For reference only (the check needs no download: the construction is on page 2 of the note and is written into
# both programs): fetch Tilley's Zenodo record 24405229 (the note and the proof archive) into certs/ and check SHA-256.
set -e
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/certs" && cd "$here/certs"
curl -s -f -L -o C49-proof.zip 'https://zenodo.org/records/24405229/files/C49-proof.zip?download=1'
sleep 10
curl -s -f -L -o tilley-2026-sunflower-capacity.pdf 'https://zenodo.org/records/24405229/files/tilley-2026-sunflower-capacity.pdf?download=1'
shasum -a 256 -c "$here/SHA256SUMS"
