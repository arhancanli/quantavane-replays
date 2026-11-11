#!/bin/sh
# Replay of three lower-bound rows of the 3d page ([FrPi73], [PeWe13], [GGSWT2025]): sum and difference sets
# counted exactly, the Penman-Wells family checked for 1 <= j <= 120, with controls.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1   # paths in the output stay relative to this folder
root=${1:-certs}
mkdir -p results
python3 -I c3d_check.py "$root" --controls > results/c3d.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/c3d.txt
exit $status
