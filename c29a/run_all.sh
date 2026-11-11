#!/bin/sh
# Replay of the 29a row [KZ1873]: the 40 roots of D_5 form a kissing configuration in R^5 (exact integers), with the
# page's cross-polytope row as an anchor and red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I d5_check.py --controls > results/kz1873.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/kz1873.txt
exit $status
