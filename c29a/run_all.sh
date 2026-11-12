#!/bin/sh
# Replays of the 29a rows: [KZ1873] (the 40 roots of D_5, exact integers) and the two linear-programming rows
# [OS1979] 46.345 and [Cox1963] 48 (our own Delsarte certificates, exact rationals, Sturm sequences and Bernstein
# coefficients), each with red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; a few seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I d5_check.py --controls > results/kz1873.txt 2>&1; s1=$?
python3 -I os1979_check.py --controls > results/os1979.txt 2>&1; s2=$?
grep -H "RESULT\|controls all" results/kz1873.txt results/os1979.txt
[ $s1 = 0 ] && [ $s2 = 0 ]
