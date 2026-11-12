#!/bin/sh
# Replay of the 89a row [Erd1962]: the random colouring gives 1/32 (all colourings of K_n enumerated for n <= 6);
# red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; a few seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I erdos1962_check.py --controls > results/erd1962.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/erd1962.txt
exit $status
