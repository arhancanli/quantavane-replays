#!/bin/sh
# Download the data of the three 3c entropy certificates replayed here (no code), each into its own new directory,
# and check their SHA-256 against SHA256SUMS. MI2026's archive also holds programs: only its 3c JSON is extracted.
# The 26-point law of [G2026] is the table in the description of PR #70; it is written out in c3c_check.py.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
new() { if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi; mkdir -p "$1"; }
# [L2026]: gist CoolRmal/5368357c... at the pinned revision 62621d2 (PR #185)
new l2026
curl -s -f -L -o l2026/certificate_3c_147pt.json \
  "https://gist.githubusercontent.com/CoolRmal/5368357cd781d7e5c676c9d68ad24d22/raw/62621d2d30df0a4405dce760cde8a5b59af51281/certificate_3c_147pt.json"
# [J2026]: Zenodo 10.5281/zenodo.23289305 (PR #221)
new j2026
curl -s -f -L -o j2026/certificate.json "https://zenodo.org/records/23289305/files/certificate.json?download=1"
# [MI2026]: Zenodo 10.5281/zenodo.20794135 (PR #93); the archive, then only c3c_pr93/certificate_3c_v3_95pt.json
new mi2026
curl -s -f -L -o mi2026/sumdiff_certs_upload.zip "https://zenodo.org/records/20794135/files/sumdiff_certs_upload.zip?download=1"
unzip -p mi2026/sumdiff_certs_upload.zip c3c_pr93/certificate_3c_v3_95pt.json > mi2026/certificate_3c_v3_95pt.json
shasum -a 256 -c "$here/SHA256SUMS"
