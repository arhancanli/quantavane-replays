#!/bin/sh
# Replay the exact values f(n) of the 5b page ([SK2026], PR #253) from its thirteen (4,5)-sets; red controls included.
# Also the upper-bound row [MT26] C_5b <= 4/7 from the 14-point set A_base of Ma and Tang (./fetch_mt26.sh).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about ten seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c5b_check.py" "$root" --controls > "$here/results/sk2026.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/sk2026.txt"
python3 -I "$here/mt26_check.py" "$root" --controls > "$here/results/mt26.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/mt26.txt"
exit $status
