#!/bin/sh
# Download, each into its own new directory, and check the SHA-256 of:
#   page/3d.md   the 3d page of teorth/optimizationproblems at main 87b4ebb (the [FrPi73] set is printed in its row);
#   ae/sums_differences_problems.ipynb   google-deepmind/alphaevolve_repository_of_problems at 8f44745 (the
#                [GGSWT2025] set is the list in its first code cell; c3d_check.py reads it with Python's ast module
#                as a literal and runs nothing from the notebook);
#   pewe/n57.pdf   Penman and Wells, Integers 13 (2013) #A57 (provenance only: c3d_check.py transcribes its family
#                Q_j and does not read the PDF).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in page ae pewe; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
curl -s -f -L -o page/3d.md "https://raw.githubusercontent.com/teorth/optimizationproblems/87b4ebbd426ff8db597a2004cc45da17903bdcf2/constants/3d.md"
curl -s -f -L -o ae/sums_differences_problems.ipynb "https://raw.githubusercontent.com/google-deepmind/alphaevolve_repository_of_problems/8f447457957deac61e28bf1676746f0753b3b2f8/experiments/sums_differences_problems/sums_differences_problems.ipynb"
curl -s -f -L -o pewe/n57.pdf "http://math.colgate.edu/~integers/n57/n57.pdf"
shasum -a 256 -c "$here/SHA256SUMS"
