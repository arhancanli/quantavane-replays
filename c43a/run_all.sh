#!/bin/sh
# Replay of the 43a row [GP1968]: the equilateral triangle, L_S = sqrt 3 (Fermat point, exact in Q(sqrt 3)) and L_M = 2,
# so C_43 <= sqrt(3)/2; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; a few seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I steiner_check.py --controls > results/gp1968.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/gp1968.txt
exit $status
