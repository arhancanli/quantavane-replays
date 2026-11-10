#!/bin/sh
# Replay #284 (26a): the six certificates and our own f_d from WRITEUP.md Section 3, with red controls;
# output in results/koshta_c26a.txt. Usage: ./run_all.sh [CERTS_ROOT]   (Python 3.9+, standard library only; seconds)
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/bh_check.py" "$root" --controls > "$here/results/koshta_c26a.txt" 2>&1
status=$?
grep -H "RESULT\|refused\|term by term" "$here/results/koshta_c26a.txt"
exit $status
