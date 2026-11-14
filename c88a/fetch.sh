#!/bin/sh
# Download the admissible tuples behind the rows of the 88a page (data only), each source into its own new directory,
# and check every SHA-256 against SHA256SUMS. Nothing downloaded is executed.
#   sutherland/   five tuples from the tables the page links (math.mit.edu/~primegaps); a web site without commits,
#                 so the hashes pin the files as retrieved on 14 November 2026
#   primegaps186/ Challenge.lean of openai/PrimeGaps186 @ 61340d0 (the 40-tuple of [OAI26]; read as text)
#   primebound/   data/certificate.json of hanxinzhang/prime-bound @ 2319570 (the 49-tuple of [Zha26]; read as JSON)
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in sutherland primegaps186 primebound; do
  if [ -d "$d" ] && [ -n "$(ls -A "$d")" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p "$d"
done
for t in 40_186 49_240 50_246 105_600 632_4680; do
  curl -s -f -L -o "sutherland/admissible_$t.txt" "https://math.mit.edu/~primegaps/tuples/admissible_$t.txt"
done
curl -s -f -L -o primegaps186/Challenge.lean \
  https://raw.githubusercontent.com/openai/PrimeGaps186/61340d0b74163003b32756bb16e91d9209a5e330/Challenge.lean
curl -s -f -L -o primebound/certificate.json \
  https://raw.githubusercontent.com/hanxinzhang/prime-bound/2319570b0056676537190512f20af80d3a7b72a7/data/certificate.json
shasum -a 256 -c "$here/SHA256SUMS"
