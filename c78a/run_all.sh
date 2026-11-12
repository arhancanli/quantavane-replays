#!/bin/sh
# Replay of the 78a row [FP2011]: straight-line thrackle drawings of the odd cycles C_5..C_101 (star polygons on a
# parabola), every pair of edges checked with integer orientation tests; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I thrackle_check.py --controls > results/fp2011.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/fp2011.txt
exit $status
