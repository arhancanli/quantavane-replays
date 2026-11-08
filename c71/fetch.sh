#!/bin/sh
# Download the data of the three 71a functions replayed here (no code), each into its own new directory, and check
# their SHA-256 against SHA256SUMS. The [MI2026b] archive also holds programs: only its certificate JSON is extracted.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
new() { if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi; mkdir -p "$1"; }
# [MI2026] (PR #94): the truth table of g14, gist 463464q435q43/d642802c... at revision db9eb64
new mi2026
curl -s -f -L -o mi2026/c71_input_table_n14.json \
  "https://gist.githubusercontent.com/463464q435q43/d642802c9a6d69bc96685006756a6fae/raw/db9eb64ba25d9ad8b25065d2f9313ed2d6d45a7c/c71_input_table_n14.json"
# [MI2026b] (PR #124): Zenodo 10.5281/zenodo.21113935, then only fei_c71_n17/certs/ridge_n17_g4_best2_iter14.json
new mi2026b
curl -s -f -L -o mi2026b/fei_c71_n17_bundle.zip "https://zenodo.org/records/21113935/files/fei_c71_n17_bundle.zip?download=1"
unzip -p mi2026b/fei_c71_n17_bundle.zip fei_c71_n17/certs/ridge_n17_g4_best2_iter14.json > mi2026b/ridge_n17_g4_best2_iter14.json
# [Num2026] (PR #130): Zenodo 10.5281/zenodo.21497769, the artifact JSON
new num2026
curl -s -f -L -o num2026/fei_c71_n18_artifact.json "https://zenodo.org/records/21497769/files/fei_c71_n18_artifact.json?download=1"
shasum -a 256 -c "$here/SHA256SUMS"
