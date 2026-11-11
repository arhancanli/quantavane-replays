#!/bin/sh
# Replay of the 26b lower-bound row [DMPSS2014] (B_R,m >= 2^(1-1/m), so C_26b >= 2): the m-linear forms T_m built by
# the paper's recursion, m = 2..7, exhaustive norms for m <= 4, with controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about half a minute).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I dmpss_check.py --controls > results/dmpss2014.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/dmpss2014.txt
exit $status
