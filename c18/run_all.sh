#!/bin/sh
# Replay of the 18a row [GT2009], C_18 >= 1.46601 (added in PR #202): the Hamming-ball ratio of Green and Tao's
# Section 5 enclosed in interval arithmetic of our own, with controls. One output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/hamming_check.py" "$root" --controls > "$here/results/gt2009.txt" 2>&1; status=$?
grep -H "RESULT\|controls all" "$here/results/gt2009.txt"
exit $status
