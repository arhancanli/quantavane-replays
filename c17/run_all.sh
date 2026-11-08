#!/bin/sh
# Replays of two 17a rows inside the theorem they rely on (Gupta-Ndiaye-Norin-Wei), with the region a bound
# certifies (both orders of k and l): [WMR2026] C_17 <= 3.6960839126* (PR #313; lambda = 1) and [Gri26]
# C_17 <= 3.791993699438612 (PR #127; Lemma 14 of v1 against Lemma 15 of v2, lambda = 0.01); red controls.
# One output file per row in results/. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs).
# Python 3.9+, standard library only; about 15 seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/ramsey_check.py" "$root" --controls > "$here/results/wmr2026.txt" 2>&1 || status=1
python3 -I "$here/ramsey_check.py" "$root" --controls --row=gri26 > "$here/results/gri26.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/wmr2026.txt" "$here/results/gri26.txt"
exit $status
