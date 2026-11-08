#!/bin/sh
# Download the data files compared here (no code): [Syl2026]'s k = 2 cost histogram and [M2026]'s data, at the commits the 3a page pins,
# into new directories, and check their SHA-256. Everything else is rebuilt from the theorem statements.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d syl2026 ] && [ -n "$(ls -A syl2026)" ]; then echo "refusing: $root/syl2026 is not empty"; exit 1; fi
mkdir -p syl2026
curl -s -f -L -o syl2026/k2_cost_hist.tsv \
  "https://raw.githubusercontent.com/xSyllaa/c3a-controlled-carry-k2/5ca7cb88820995bf7666096773fd3799f5714b18/data/k2_cost_hist.tsv"
# [Syl2026]'s k = 3 theorem note, certificate summary and cost histogram (compared only; the histogram is rebuilt by
# c3a_hist3.c).
if [ -d syl2026k3 ] && [ -n "$(ls -A syl2026k3)" ]; then echo "refusing: $root/syl2026k3 is not empty"; exit 1; fi
mkdir -p syl2026k3
for f in proof/THEOREM_k.md data/k3_cost_hist.tsv certificates/cert_k3_primary_ell_1e24.json; do
  curl -s -f -L -o "syl2026k3/$(basename $f)" \
    "https://raw.githubusercontent.com/xSyllaa/c3a-controlled-carry-k2/5ca7cb88820995bf7666096773fd3799f5714b18/k3/$f"
done
# [M2026]'s data files (its finite certificates, mask, sum support and pruned-difference table) at the commit the 3a page
# pins, compared with what we rebuild; none of its code is fetched.
if [ -d m2026 ] && [ -n "$(ls -A m2026)" ]; then echo "refusing: $root/m2026 is not empty"; exit 1; fi
mkdir -p m2026
for f in finite_certificate_17086_27995.json finite_certificate_17032_27995.json kleinwaks_finite_certificate_17032_34065.json \
         mask_17086.txt sum_support_17086.txt pruned_differences_17086_27995.tsv; do
  curl -s -f -L -o "m2026/$f" \
    "https://raw.githubusercontent.com/mizopl/ghr-c3a-odd-base-finite-bound/464ceda9058798f0070ab92d35cd1c51865eb0ec/object/$f"
done
# [GGSWT2025]'s two sets U_1, U_2 (the 3a rows 1.1479 and 1.1584): the AlphaEvolve notebook at commit 8f44745, read as
# data (ae_check.py extracts the two list literals with Python's ast module; nothing in it is run)
if [ ! -f ae/sums_differences_problems.ipynb ]; then
  if [ -d ae ] && [ -n "$(ls -A ae)" ]; then echo "refusing: $root/ae is not empty"; exit 1; fi
  mkdir -p ae
  curl -s -f -L -o ae/sums_differences_problems.ipynb \
    "https://raw.githubusercontent.com/google-deepmind/alphaevolve_repository_of_problems/8f44745/experiments/sums_differences_problems/sums_differences_problems.ipynb"
fi
shasum -a 256 -c "$here/SHA256SUMS"
