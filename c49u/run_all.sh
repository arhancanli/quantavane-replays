#!/bin/sh
# Replay of the counting step of the 49a row [NS2017]: sum_{k <= n/3} C(n, k) <= (27/4)^(n/3), n = 1..900, exact;
# red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; a few seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I ns2017_check.py --controls > results/ns2017.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/ns2017.txt
exit $status
