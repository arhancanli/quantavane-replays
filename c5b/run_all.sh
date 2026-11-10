#!/bin/sh
# Replay the exact values f(n) of the 5b page ([SK2026], PR #253) from its thirteen (4,5)-sets; red controls included.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about ten seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c5b_check.py" "$root" --controls > "$here/results/sk2026.txt" 2>&1
status=$?
grep -H "RESULT" "$here/results/sk2026.txt"
exit $status
