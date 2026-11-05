#!/bin/sh
# Download the C68 certificate data file replayed here (no code) into its own new directory, at the commit the 68a
# page pins, and check its SHA-256 against SHA256SUMS.
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
# [Pin2026], C68 >= 0.432* (68a page; added in PR #239)
get pin certificate-jung.json rpinho/korenblum-c68-certificate 9d48de64987911f5a07a397323345be8b97770a8 certificate-jung.json
# PR #274 (vidigoat): C68 >= 0.446 and the supplementary 0.44605 (flux and half-diameter lemmas)
get vid0446 certificate-c0446.json vidigoat/optimization-records 1a76f9dd053d51cce531b57b2ad5a68e3205ca75 \
    tao/certificates/c68/certificate-c0446.json
get vid044605 certificate-c044605.json vidigoat/optimization-records 1a76f9dd053d51cce531b57b2ad5a68e3205ca75 \
    tao/certificates/c68/certificate-c044605.json
shasum -a 256 -c "$here/SHA256SUMS"
