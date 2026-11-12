#!/bin/sh
# Replay of the 4a row [P1970]: a 20-point cap in F_3^4 (found by our own search from the affine frame), checked by two
# routes, its Cartesian square a cap in F_3^8, completeness, and 20^(1/4) > 2.1146 exactly; red controls; no download.
# Also the [EG2016] row C_4a <= 2.756 (eg2016_check.py: the base of the monomial count M_n, exact at x = 5931/10000).
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about five seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I cap_check.py --controls > results/p1970.txt 2>&1; status=$?
python3 -I eg2016_check.py --controls > results/eg2016.txt 2>&1 || status=1
grep -H "RESULT\|controls all" results/p1970.txt results/eg2016.txt
exit $status
