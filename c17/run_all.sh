#!/bin/sh
# Replay of the 17a row [WMR2026], C_17 <= 3.6960839126* (PR #313): the printed value, the benchmark's acceptance
# rule at lambda = 1, and the conditions of Gupta-Ndiaye-Norin-Wei's theorem at lambda = 1 with the region a bound
# certifies (both orders of k and l), for three bounds; red controls. One output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 5 seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/ramsey_check.py" "$root" --controls > "$here/results/wmr2026.txt" 2>&1; status=$?
grep -H "RESULT\|controls all" "$here/results/wmr2026.txt"
exit $status
