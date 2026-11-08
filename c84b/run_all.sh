#!/bin/sh
# Replay the 84b row 1.9993* ([Al26]: the note's chain at its parameters) and PR #275's ceiling over every choice of
# s, X, Y, with red controls; one output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 25 seconds.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results
python3 -I althoefer_check.py "$root/al26" --controls > results/al26.txt 2>&1; status=$?
grep -H "RESULT\|so c >=\|ceiling:\|rho(s) >=" results/al26.txt
exit $status
