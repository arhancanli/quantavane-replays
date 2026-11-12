#!/bin/sh
# Replay of the 76a row [BS2007]: 3d^2 distinct lines on the smooth Fermat surface of degree d, d = 3..12, exact in
# cyclotomic fields; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about 15 seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I fermat_lines_check.py --controls > results/bs2007.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/bs2007.txt
exit $status
