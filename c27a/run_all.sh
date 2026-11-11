#!/bin/sh
# Replay of the 27a row [MM1961]: the Moser spindle, exact coordinates in Q(sqrt 3, sqrt 11), 11 unit distances,
# chromatic number 4 by exhaustion; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I spindle_check.py --controls > results/mm1961.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/mm1961.txt
exit $status
