#!/bin/sh
# Replay of the 36a row [dLOV2014] (lower): the D_4 lattice packing, density pi^2/16 exactly; and the conversions of
# the three upper-bound rows with the README cell of PR #364; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I d4_check.py --controls > results/dlov2014.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/dlov2014.txt
exit $status
