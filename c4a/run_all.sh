#!/bin/sh
# Replay of the 4a row [P1970]: a 20-point cap in F_3^4 (found by our own search from the affine frame), checked by two
# routes, its Cartesian square a cap in F_3^8, completeness, and 20^(1/4) > 2.1146 exactly; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I cap_check.py --controls > results/p1970.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/p1970.txt
exit $status
