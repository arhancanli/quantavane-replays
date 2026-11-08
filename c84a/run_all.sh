#!/bin/sh
# Replay the 84a row 1.0152* (PR #315) on Sawin's criterion, with red controls; one output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/sawin_check.py" "$root" --controls > "$here/results/emm2026.txt" 2>&1; status=$?
grep -H "RESULT\|run time" "$here/results/emm2026.txt"
exit $status
