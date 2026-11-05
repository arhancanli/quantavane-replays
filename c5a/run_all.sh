#!/bin/sh
# Replay the 5a bound C_5a <= 0.9428095* through the certificate of PR #335; red controls included.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c5a_check.py" var2026 "$root" --controls > "$here/results/var2026.txt" 2>&1
status=$?
grep -H "RESULT\|run time" "$here/results/var2026.txt" | tr '\n' ' '; echo
exit $status
