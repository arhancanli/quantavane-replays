#!/bin/sh
# Replay the 5a bound C_5a <= 0.9428095* through the certificate of PR #335; red controls included.
# Also the row [HZ26] C_5a <= 0.9435 (Hou and Zhao's eight-pair kernel-boundary system; needs ./fetch_hz26.sh).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c5a_check.py" var2026 "$root" --controls > "$here/results/var2026.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/var2026.txt" | tr '\n' ' '; echo
python3 -I "$here/hz26_check.py" "$root" --controls > "$here/results/hz26.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/hz26.txt"
exit $status
