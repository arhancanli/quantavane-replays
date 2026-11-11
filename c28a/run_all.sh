#!/bin/sh
# Replay of the 28a row [Gri2026], C_28 <= 63 (PR #84): the G_2(4) graph built from Brouwer's projective model, its
# strongly regular parameters and clique number, the partition and degree data, Lemmas 3-6 exactly, with controls.
# Also the rows [Bon2014] (C_28 <= 65) and [JB2014] (C_28 <= 64) on G_2(4) built a second way, from Jenrich's paper
# (bon_jb_check.py; needs ./fetch_bonjb.sh for the two LaTeX sources, read as text).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 15 seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/borsuk_check.py" --controls > "$here/results/gri2026.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/gri2026.txt"
python3 -I "$here/bon_jb_check.py" "$root" --controls > "$here/results/bon2014_jb2014.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/bon2014_jb2014.txt"
exit $status
