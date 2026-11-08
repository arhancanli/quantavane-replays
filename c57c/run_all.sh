#!/bin/sh
# Replay of the finite certificate (Lemma 6.1) behind the 57c row [W2026], U >= 0.612242 (Wikstrom, arXiv:2610.34580v1,
# Theorem 1.2; recorded in PR #311), from the paper's statements alone (no download), with the red controls.
# Usage: ./run_all.sh. Python 3.9+, standard library only; about 75 minutes, under 30 MB.
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/results"
python3 -I -u "$here/w2026_check.py" --controls > "$here/results/w2026.txt" 2>&1; status=$?
grep -H "RESULT" "$here/results/w2026.txt"
exit $status
