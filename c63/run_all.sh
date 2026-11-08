#!/bin/sh
# Replay the 63a row theta* = 0.3144831759741 (PR #313) and the 64a row 2 theta* (PR #315) from Definition 1.1 of
# Li and Yang (arXiv 2308.14859v2), with red controls; one output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a second.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/theta_check.py" "$root/ly2023" --controls > "$here/results/ly2023.txt" 2>&1; status=$?
grep -H "RESULT\|printed" "$here/results/ly2023.txt" | grep -v ":\[ctl"
exit $status
