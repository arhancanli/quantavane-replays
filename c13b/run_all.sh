#!/bin/sh
# Replay of the two 13b rows [Pal1920]: the regular hexagon of width 1 (area sqrt(3)/2) and Pal's truncation by two
# corners of the circumscribed dodecagon (area 2 - 2/sqrt 3), exact in Q(sqrt 3); red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I pal_check.py --controls > results/pal1920.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/pal1920.txt
exit $status
