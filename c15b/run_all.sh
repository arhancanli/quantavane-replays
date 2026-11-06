#!/bin/sh
# Replay the 15b bounds alpha > 0.46652506* ([Var2026]) and alpha > 0.465* ([OAI2026]) of PR #337; red controls included.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c15b_check.py" var2026 "$root" --controls > "$here/results/var2026.txt" 2>&1
status=$?
grep -H "RESULT\|run time" "$here/results/var2026.txt" | tr '\n' ' '; echo
exit $status
