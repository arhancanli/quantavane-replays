#!/bin/sh
# Replay of the 15a row [S1969]: Strassen's seven products, the bilinear identity over any ring, the recursion
# (7^k multiplications for 2^k x 2^k, k <= 6) and log_2 7; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I strassen_check.py --controls > results/s1969.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/s1969.txt
exit $status
