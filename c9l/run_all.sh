#!/bin/sh
# Replay of the 9a row [BMRRST1971]: a 343-word linear code in Z_7^5 independent in the strong power C_7^(x5)
# (our own, found by exhaustive search; the paper is not accessible to us), checked by two routes; red controls.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I bmrrst_check.py --controls > results/bmrrst1971.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/bmrrst1971.txt
exit $status
