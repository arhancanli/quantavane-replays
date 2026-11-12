#!/bin/sh
# Replay of the 41a row [Hammersley1968]: Hammersley's sofa, the motion (exact path identities in s = tan(t/2), numeric
# anchor over 2001 angles) and the area pi/2 + 2/pi (exact, pi as a symbol); red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about ten seconds).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I hammersley_check.py --controls > results/hammersley1968.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/hammersley1968.txt
exit $status
