#!/bin/sh
# Replay of the 27b rows [R1959] (K_8) and [G1980] (the join of K_6 and C_5): each split into two planar graphs, the
# embeddings verified (darts, rotation system, Euler's formula), chromatic numbers 8 and 9 by exhaustion; controls.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about a second; no download).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I biplanar_check.py --controls > results/r1959_g1980.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/r1959_g1980.txt
exit $status
