#!/bin/sh
# Replay of the 9a upper-bound rows [L1979] (Lovasz: theta(C_7) = 7 cos(pi/7)/(1 + cos(pi/7)) ~ 3.3177, exact spectrum
# and the matrix certificate) and [S1956] (Shannon: the fractional clique cover, 7/2); red controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I theta_check.py --controls > results/l1979_s1956.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/l1979_s1956.txt
exit $status
