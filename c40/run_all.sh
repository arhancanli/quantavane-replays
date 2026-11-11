#!/bin/sh
# Replay of the upper-bound row of the 40a page ([Leh1933]/[BDM2007], Lehmer's polynomial):
# M(l) enclosed exactly by two routes (Sturm sequences and bisection over the rationals), with controls.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about a second.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1   # paths in the output stay relative to this folder
root=${1:-certs}
mkdir -p results
python3 -I lehmer_check.py "$root" --controls > results/c40.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/c40.txt
exit $status
