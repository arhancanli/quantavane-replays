#!/bin/sh
# Download the two data files replayed here (no code), each into its own new directory, and check their SHA-256
# against SHA256SUMS: the four certificates of Emmerich's report (arXiv 2606.03419; Sawin's published example, the
# greedy, the evolution-strategy and the discrete-recombination certificates) at the repository commit 2279d0a, and
# the data file mlewko's post in Erdos Problems forum thread 90 links (21 May 2026: the sets T and S, the weights k(p)
# and an inertness witness table; Google Drive file 1ygqjOSatMDYe439YMP8VJA7ZT5XVDBvG).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in emm epf90; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o emm/certificate_verification_all_v23.json \
  "https://raw.githubusercontent.com/emmerichmtm/UnitDistanceProblemOptimizationOfSawinsLowerBound/2279d0a277502a8bbfc770630724f34a72f0c73a/certificate_verification_all_v23.json"
curl -s -f -L -o epf90/mlewko_sawin_data.txt \
  "https://drive.google.com/uc?export=download&id=1ygqjOSatMDYe439YMP8VJA7ZT5XVDBvG"
shasum -a 256 -c "$here/SHA256SUMS"
