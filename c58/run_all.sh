#!/bin/sh
# Replay of the 58a lower-bound row [Kan2021] (A_Zar >= 5: the A = 4 version fails at d = 54 and d = 150), exactly,
# with controls; no download. Usage: ./run_all.sh   (Python 3.9+, standard library only; a few seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I zaremba_check.py --controls > results/kan2021.txt 2>&1; status=$?
grep -H "RESULT\|anchor\|controls all" results/kan2021.txt
exit $status
