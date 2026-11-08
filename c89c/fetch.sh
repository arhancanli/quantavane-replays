#!/bin/sh
# Download the 89c certificate and graph list of PR #308 (rpinho/math-certificates @ 603b40c, c89c/round3: data only,
# no code) into a new directory, decompress the certificate with the standard library, and check the SHA-256 of all
# three files against SHA256SUMS. Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d c89c ] && [ -n "$(ls -A c89c)" ]; then echo "refusing: $root/c89c is not empty"; exit 1; fi
mkdir -p c89c
P=https://raw.githubusercontent.com/rpinho/math-certificates/603b40c6fee0f0c8a5b6cfa2b07f6f96545fb857/c89c/round3
curl -s -f -L -o c89c/c6_ordered_type8.cert.xz "$P/certificates/c6_ordered_type8.cert.xz"
curl -s -f -L -o c89c/graph8.g6 "$P/sources/graph8.g6"
python3 -I -c "import lzma, sys; open(sys.argv[2], 'wb').write(lzma.open(sys.argv[1]).read())" \
    c89c/c6_ordered_type8.cert.xz c89c/c6_ordered_type8.cert
shasum -a 256 -c "$here/SHA256SUMS"
