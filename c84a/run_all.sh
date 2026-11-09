#!/bin/sh
# Replay the 84a rows on Sawin's criterion, with red controls: the row 1.0152* (PR #315; results/emm2026.txt), and Sawin's
# own example [S2026] 1.014 with the three finite-data rows [EPF90] 1.03184..., [MO511514] 1.03188... and 1.03190... (results/forum2026.txt).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/sawin_check.py" "$root" --controls > "$here/results/emm2026.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/emm2026.txt"
python3 -I "$here/forum_check.py" "$root" --controls > "$here/results/forum2026.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/forum2026.txt"
exit $status
