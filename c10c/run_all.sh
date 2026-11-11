#!/bin/sh
# Replay of the four certificate rows of the 10c page ([G2026] #74, [X2026] #86, [L2026] #129, [H2026] #170): each
# matrix read from the page, its discrepancy computed exactly by two exhaustive routes, with controls.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c10c_check.py" "$root" --controls > "$here/results/c10c.txt" 2>&1; status=$?
grep -H "RESULT\|controls all" "$here/results/c10c.txt"
exit $status
