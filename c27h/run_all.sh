#!/bin/sh
# Replay of the 27a row [Had1945], [CR2017]: the hexagon 7-colouring of the plane at side 2/5, exact (rationals), with
# the colour sublattice enumerated; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I hexagon_check.py --controls > results/had1945.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/had1945.txt
exit $status
