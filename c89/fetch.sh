#!/bin/sh
# Download the C89 data files replayed here (no code), each into a new directory, at the commits PRs #279, #307 and
# #277 pin, and the PPSS control graphs from Zenodo; then check their SHA-256 against SHA256SUMS. PR #331's
# construction is replayed from the eleven n_w stated in the pull request: its certificate record (Zenodo
# 10.5281/zenodo.23974391) is restricted, so there is nothing to fetch for it.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # get DIR FILE REPO COMMIT PATH: one file into a directory that must be new or empty
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "https://raw.githubusercontent.com/$3/$4/$5"
}
# PR #279, [SK2026] 0.7422548616*2^-14 (JeetuSK0808/tao-constants-certificates @ 98b6ee3)
get sk2026 construction.json JeetuSK0808/tao-constants-certificates 98b6ee34e8f51f82be85517e4b7f0e7e7af393b8 \
    c89c-k6-ramsey-multiplicity/construction.json
# PR #307, [Pin2026] C89b < 0.001648970526* (rpinho/math-certificates @ 1d38853, c89b/round4/certificates): the
# frozen full matrix and the two-mask definition
get pr307 c5_cayley256_fine.json.gz rpinho/math-certificates 1d38853f31112da873c632fb3d3bd05519eb832a \
    c89b/round4/certificates/c5_cayley256_fine.json.gz
curl -s -f -L -o pr307/c5_cayley256_definition.json \
    "https://raw.githubusercontent.com/rpinho/math-certificates/1d38853f31112da873c632fb3d3bd05519eb832a/c89b/round4/certificates/c5_cayley256_definition.json"
# PR #277, [SK2026] C89b <= 0.0016870052* (JeetuSK0808/tao-constants-certificates @ 82a6e4d): the 192-block matrix
get pr277 c89b_stepgraphon_m192_D20.json.gz JeetuSK0808/tao-constants-certificates \
    82a6e4def9c63e9a7254bc2393f61a872891c480 c89b-k5-ramsey-multiplicity/c89b_stepgraphon_m192_D20.json.gz
# control: the Parczyk-Pokutta-Spiegel-Szabo graphs (Zenodo record 6602512, CC BY 4.0); c5.graph6.txt is read inside
if [ -d ppss ] && [ -n "$(ls -A ppss)" ]; then echo "refusing: $root/ppss is not empty"; exit 1; fi
mkdir -p ppss && curl -s -f -L -o ppss/graphs.zip "https://zenodo.org/api/records/6602512/files/graphs.zip/content"
cd "$root" && shasum -a 256 -c "$here/SHA256SUMS"
