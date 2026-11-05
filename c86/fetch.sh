#!/bin/sh
# Download the C86 certificate data files replayed here (no code), each into its own new directory, at the commits
# the 86a page and the pull requests pin, and check their SHA-256 against SHA256SUMS.
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
# [Kos2026], 1.8084990* (86a page; PR #292)
get kos c86_A19d7_robust_1.8084990.json DennisKoshta/c86-sss-trace-certificate \
    90c57242b2f43aa942257154be6731c29898b537 certificate/c86_A19d7_robust_1.8084990.json
# [Pin26], 1.8117251* (86a page; PR #242)
get pin cert_n1200_A11.json rpinho/c86-trace-certificate 116a58ff99681226ce6299d4a02eeff1e7d3c5c7 cert_n1200_A11.json
# [Pat2026] row of PR #270, 1.8077686* (the file has the same git blob at fb80010, the PR body, and 7c6586f, the row)
get vid c86_deg8_h400_1.807768536863.json vidigoat/optimization-records fb800103581d5c63fa096805f73eff48fce32c98 \
    tao/certificates/c86/c86_deg8_h400_1.807768536863.json
# PR #298, 1.80560940* (piecewise-linear density on 5,375 nodes; auxiliaries with nonreal roots)
get p298 certificate.json rpinho/math-certificates 0d821621bab3d3d9408015f2679b89d5d5e9b379 c86/adaptive-cone/certificate.json
# PR #270 update of 5 Nov (vidigoat c86_v2), 1.8055311* (the format of #298; 31 auxiliaries, 5,734 knots)
get vid2 certificate.json vidigoat/optimization-records fd4373c76ec2fa9969501fa9a9663773e8668aa0 \
    tao/certificates/c86_v2/certificate.json
# PR #270 update of 5 Nov, night (vidigoat c86_v3), 1.8054956* (#298's format; 33 auxiliaries, 6,999 knots)
get vid3 certificate.json vidigoat/optimization-records 9e2b86ca420ab00390cf9b81b4002da0a0248001 \
    tao/certificates/c86_v3/certificate.json
# arXiv:2609.29298 (D. Niedbala Giraudin), the lower bound 1.80220* of PR #315: the certificate, an ancillary file of v1
if [ -d ng ] && [ -n "$(ls -A ng)" ]; then echo "refusing: $root/ng is not empty"; exit 1; fi
mkdir -p ng && curl -s -f -L -o ng/CERTIFICAT_C86.json "https://arxiv.org/src/2609.29298v1/anc/CERTIFICAT_C86.json"
shasum -a 256 -c "$here/SHA256SUMS"
