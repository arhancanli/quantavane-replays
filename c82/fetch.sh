#!/bin/sh
# Download the C82 certificate data files replayed here (no code), each into its own new directory, at the commits
# the 82a page and the pull requests pin, and check their SHA-256 against SHA256SUMS.
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
# [Pin26], 0.2534449* (82a page: certificate SHA-256 b33626f3..., commit 69d8964; PR #242)
get pin certificate.json rpinho/c82-height-certificate 69d89640dffd8bf3e2067dab769a54c7b8e8c490 certificate.json
# [Epi26], 0.2535744749068227* (82a page; the only commit of the folder; result.json states the artifact SHA-256)
get epi construction.json BakeLab/Epi-results 97fa6a420f417458e11723c8bceb70f1529e447a \
    records/zhang-zagier/zhang_zagier_c82/construction.json
# PR #280, 0.2534350402* (vidigoat/optimization-records @ 32898ac, tao/certificates/c82)
get pr280 certificate.json vidigoat/optimization-records 32898ac1dc0d75bf7f6d8a4e0a779eea629467a4 \
    tao/certificates/c82/certificate.json
# PR #280 updated 4 Nov, 0.2524177045* (vidigoat/optimization-records @ 3a7a84f, tao/certificates/c82/v2)
get pr280b certificate.json vidigoat/optimization-records 3a7a84fec9c64ae426f06a5f52597c84b3453ba0 \
    tao/certificates/c82/v2/certificate.json
# PR #280 updated 5 Nov, 0.2520591386* (vidigoat/optimization-records @ 88c39d5, tao/certificates/c82/v3)
get pr280c certificate.json vidigoat/optimization-records 88c39d53109a808a0e0315ed75d5868e582b8e57 \
    tao/certificates/c82/v3/certificate.json
# PR #301, 0.253091834* (rpinho/c82-height-certificate @ 6d4e80c, degree44-denominators-v2) and its first filing
get pr301 certificate.json rpinho/c82-height-certificate 6d4e80cb0b375c767df91eae7985ae3e8bef7203 \
    degree44-denominators-v2/certificate.json
get pr301a certificate.json rpinho/c82-height-certificate 6d4e80cb0b375c767df91eae7985ae3e8bef7203 \
    degree44-denominators/certificate.json
# positive control: Doche 2001b's polynomials and weights as [Pin26] records them; red control file of [Pin26]
get doche doche2001.json rpinho/c82-height-certificate 69d89640dffd8bf3e2067dab769a54c7b8e8c490 controls/doche2001.json
get pinm mutated_v3_should_fail.json rpinho/c82-height-certificate 69d89640dffd8bf3e2067dab769a54c7b8e8c490 \
    controls/mutated_v3_should_fail.json
shasum -a 256 -c "$here/SHA256SUMS"
