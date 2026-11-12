#!/bin/sh
# Replay of the 69a row [Tao2022]: z^n - 1, lambda_0 = 1, n = 2..60, exact; red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I sendov_example_check.py --controls > results/tao2022.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/tao2022.txt
exit $status
