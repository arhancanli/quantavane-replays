#!/bin/sh
# Download the 49a certificate of PR #370 (data only, no code), each file into a new directory, and check the SHA-256:
# rpinho/math-certificates @ a50021b, c49/round1/certificates/best_1.5678643845728069793620.json (the 252,303-set
# block) and best_delta.json (its stated difference from Koshta's block), and, for program A's second path,
# Koshta's block itself (DennisKoshta/sunflower3-49a-certificate @ c4ee817, as in ../c49/fetch.sh; xz unpacked with
# Python's lzma module).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in pr370 kos; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
base=https://raw.githubusercontent.com/rpinho/math-certificates/a50021b04687719f4a2ea27609029982c980573e/c49/round1/certificates
curl -s -f -L -o pr370/best_1.5678643845728069793620.json "$base/best_1.5678643845728069793620.json"
curl -s -f -L -o pr370/best_delta.json "$base/best_delta.json"
curl -s -f -L -o kos/two_pendant_extra_1.5678599.txt.xz \
  https://raw.githubusercontent.com/DennisKoshta/sunflower3-49a-certificate/c4ee81737f9f6a4ae9793586a27004553a12d5c0/two_pendant_extra_1.5678599.txt.xz
python3 -I -c "import lzma, sys; open(sys.argv[2], 'wb').write(lzma.open(sys.argv[1]).read())" \
  kos/two_pendant_extra_1.5678599.txt.xz kos/two_pendant_extra_1.5678599.txt
shasum -a 256 -c "$here/SHA256SUMS"
