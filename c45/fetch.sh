#!/bin/sh
# Download the base histograms of vidigoat's C45 packet (data only; its generators and evaluators are not fetched), each
# into its own new directory, at the commit PR #281 links, and check their SHA-256 against SHA256SUMS. They are read
# only to compare them with the histograms we compute.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # get DIR FILE: one file into a directory that must be new or empty
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "https://raw.githubusercontent.com/vidigoat/optimization-records/eff83aa564d9cc74e98d6fbae308b877ad692a4a/tao/certificates/c45/data/$2"
}
get griego base_griego13.txt
get yoo base_yoo14.txt
get pin21 base_final.txt
shasum -a 256 -c "$here/SHA256SUMS"
