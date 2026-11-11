#!/bin/sh
# Download the data replayed here: the Zenodo record of [SK2026] on the 17b page (PR #252; doi:10.5281/zenodo.23552014,
# r3_33_ramsey_certificate.zip), into a new directory; check its SHA-256; extract only the two certificate files (JSON
# data: connection set, edge list, group table, stabilizer maps). The bundle's checker/, verification_main/ and src/
# programs are not extracted, read or run.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d sk2026 ] && [ -n "$(ls -A sk2026)" ]; then echo "refusing: $root/sk2026 is not empty"; exit 1; fi
mkdir -p sk2026
curl -s -f -L -o sk2026/r3_33_ramsey_certificate.zip \
  "https://zenodo.org/records/23552014/files/r3_33_ramsey_certificate.zip?download=1"
python3 -I -c "
import zipfile, sys
z = zipfile.ZipFile('sk2026/r3_33_ramsey_certificate.zip')
for name in ('R3_33_n240_Z240.json', 'R3_34_n240_Z240.json'):
    open('sk2026/' + name, 'wb').write(z.read('r3_33_bundle/certificates/' + name))
"
shasum -a 256 -c "$here/SHA256SUMS"
