#!/bin/sh
# Replay of the 28a row [Gri2026], C_28 <= 63 (PR #84): the G_2(4) graph built from Brouwer's projective model, its
# strongly regular parameters and clique number, the partition and degree data, Lemmas 3-6 exactly, with controls.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; no download needed; about 10 seconds)
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/results"
python3 -I "$here/borsuk_check.py" --controls > "$here/results/gri2026.txt" 2>&1; status=$?
grep -H "RESULT\|controls all" "$here/results/gri2026.txt"
exit $status
