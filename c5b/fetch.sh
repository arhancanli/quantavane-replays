#!/bin/sh
# Download the thirteen data files replayed here (no code): the (4,5)-sets of [SK2026] on the 5b page (PR #253;
# JeetuSK0808/tao-constants-certificates at the commit the page links), into a new directory, and check their SHA-256
# against the hashes the repository states. Its verify.py and check_all.sh are not downloaded.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d sk2026 ] && [ -n "$(ls -A sk2026)" ]; then echo "refusing: $root/sk2026 is not empty"; exit 1; fi
mkdir -p sk2026
base=https://raw.githubusercontent.com/JeetuSK0808/tao-constants-certificates/afa2e1010cd5f36cfbad213dc83b9568fe77f934/c5b-erdos757-small-n/certificates
for n in 6 7 8 9 10 11 12 13 14 15 17 19 21; do
  curl -s -f -L -o sk2026/n$n.txt "$base/n$n.txt"
done
shasum -a 256 -c "$here/SHA256SUMS"
