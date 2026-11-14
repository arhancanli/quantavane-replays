#!/bin/sh
# Download the step-function certificates replayed here (data only, no code), each into a new directory, and check
# their SHA-256: the two 1a files of PR #367 (the row's file at 9edc297 and the file of the PR description at
# 8c9fa72), the 1a ancillary file of [Cum2026] that PR #367 records (xz, unpacked with Python's lzma module; the
# unpacked file is named by its own SHA-256), and the 1b file of PR #368.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # get DIR FILE URL: one file into a directory that must be new or empty
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "$3"
}
# PR #367 (1a), the row in its files: intervalplace/frontier @ 9edc297, c1-first-autocorrelation
get pr367 C1_262144_1.5026648332.txt \
  https://raw.githubusercontent.com/intervalplace/frontier/9edc297b5d1c7a6e22d7d7d00c7278923ca68fc8/c1-first-autocorrelation/C1_262144_1.5026648332.txt
# PR #367 (1a), the certificate its description and title still name: intervalplace/frontier @ 8c9fa72
get pr367d C1_262144_1.5026648336.txt \
  https://raw.githubusercontent.com/intervalplace/frontier/8c9fa72af12796c00a82066ede262655ffc3abd3/c1-first-autocorrelation/C1_262144_1.5026648336.txt
# [Cum2026] (recorded by PR #367): the ancillary file of L. S. Cumes's paper, as served on 14 November 2026 (no commit)
get cum2026 49c83eb6.xz \
  https://qwntl.com/paper-new-constructions/anc/49c83eb6f9f82e244578fab15d056757886679700bd8c8c653f0b100ae882f2d.xz
python3 -I -c "import lzma, sys; open(sys.argv[2], 'wb').write(lzma.open(sys.argv[1]).read())" \
  cum2026/49c83eb6.xz cum2026/49c83eb6f9f82e244578fab15d056757886679700bd8c8c653f0b100ae882f2d
# PR #368 (1b): intervalplace/frontier @ 6d6c4ba, c1b-erdos-minimum-overlap
get pr368 C1b_15360_0.3808564712.json \
  https://raw.githubusercontent.com/intervalplace/frontier/6d6c4bac7b04cc224c1bb6016f594e5b785e3f48/c1b-erdos-minimum-overlap/C1b_15360_0.3808564712.json
shasum -a 256 -c "$here/SHA256SUMS"
