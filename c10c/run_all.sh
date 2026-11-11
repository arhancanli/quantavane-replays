#!/bin/sh
# Replay of the four certificate rows of the 10c page ([G2026] #74, [X2026] #86, [L2026] #129, [H2026] #170): each
# matrix read from the page, its discrepancy computed exactly by two exhaustive routes, with controls.
# Also the classical row [Band2024] (band2024_check.py).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c10c_check.py" "$root" --controls > "$here/results/c10c.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/c10c.txt"
# the classical row [Band2024] (the 2 x 2 sign matrix printed in its row): disc 2 exactly, so sqrt 2; controls
python3 -I "$here/band2024_check.py" "$root" --controls > "$here/results/band2024.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/band2024.txt"
exit $status
