#!/bin/sh
# Download the one data file replayed here (no code): the four certificates of Emmerich's report (arXiv 2606.03419;
# Sawin's published example, the greedy, the evolution-strategy and the discrete-recombination certificates) at the
# repository commit 2279d0a, into its own new directory, and check its SHA-256 against SHA256SUMS.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d emm ] && [ -n "$(ls -A emm)" ]; then echo "refusing: $root/emm is not empty"; exit 1; fi
mkdir -p emm
curl -s -f -L -o emm/certificate_verification_all_v23.json \
  "https://raw.githubusercontent.com/emmerichmtm/UnitDistanceProblemOptimizationOfSawinsLowerBound/2279d0a277502a8bbfc770630724f34a72f0c73a/certificate_verification_all_v23.json"
shasum -a 256 -c "$here/SHA256SUMS"
