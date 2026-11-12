#!/bin/sh
# Replay of the 53a row [GG2006]: e1^(n-1) e2^(n-1) e3^(n-1) is zero-sum free in C_n^3 (n = 2..9, two routes for
# n <= 5), so D(C_n^3) >= 3(n-1) + 1; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I gg2006_check.py --controls > results/gg2006.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/gg2006.txt
exit $status
