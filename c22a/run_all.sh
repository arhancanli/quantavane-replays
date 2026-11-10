#!/bin/sh
# Replay the finite part of [V2026] C_22a >= 2.124 (#356) with red controls; output in results/v2026.txt.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about ten seconds)
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/results"
python3 -I "$here/c22a_check.py" --controls > "$here/results/v2026.txt" 2>&1
status=$?
grep -H "RESULT\|control" "$here/results/v2026.txt"
exit $status
