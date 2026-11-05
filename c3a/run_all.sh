#!/bin/sh
# Replay the 3a row [Syl2026] k = 2 (C_3a > 1.19281555*) exactly from its theorem note; red controls included.
# Usage: ./run_all.sh [CERTS_ROOT] [--workers=N]   (default: ./certs). Python 3.9+, standard library only, and cc.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
shift 2>/dev/null
mkdir -p "$here/results"
python3 -I "$here/c3a_check.py" syl2 "$root" --controls "$@" > "$here/results/syl2.txt" 2>&1
status=$?
grep -H "RESULT\|run time" "$here/results/syl2.txt" | tr '\n' ' '; echo
exit $status
