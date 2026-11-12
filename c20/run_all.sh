#!/bin/sh
# Replay of the 20a and 20c rows [Pat2026]: the one-dimensional examples E - 1 (Var(X^2) = 8) and the uniform law on
# [-sqrt 3, sqrt 3] (psi >= sqrt 3), exact; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I pat2026_check.py --controls > results/pat2026.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/pat2026.txt
exit $status
