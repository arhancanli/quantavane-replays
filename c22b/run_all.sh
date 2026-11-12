#!/bin/sh
# Replay of [Claude2026] C_22b <= 7.2012719* (PR #267) with red controls; output in results/claude2026.txt.
# Usage: ./run_all.sh [CERTS_ROOT]   (Python 3.9+, standard library only; about ten minutes on one core.) If
# CERTS_ROOT holds the files of fetch.sh, design.json and the q strings are compared with the script's own copies.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
src=""
[ -f "$root/design.json" ] && src="--sources $root"
python3 -I -u "$here/c22b_check.py" --controls $src > "$here/results/claude2026.txt" 2>&1
status=$?
grep -H "RESULT\|controls all" "$here/results/claude2026.txt"
exit $status
