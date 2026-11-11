#!/bin/sh
# Replay R(3,33) >= 241 ([SK2026] on the 17b page) from its certificate data: structure, triangle-freeness, the
# multiplier automorphisms and their orbits, red controls, and the exact search for an independent 33-set.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). A C compiler and Python 3.9+. The K = 33 search is the heavy
# part (single core, a few MB; about an hour or more), one orbit per process and resumable through
# results/orbits_K33.log; run it under nice.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results" "$here/bin"
cc -O2 -o "$here/bin/c17b_indep" "$here/c17b_indep.c" || exit 1
polite=; [ -n "$POLITE" ] && polite=--polite   # POLITE=1: wait while the machine is busy (see c17b_check.py)
python3 -I "$here/c17b_check.py" "$root" "$here/bin" --controls --main $polite > "$here/results/sk2026_r3_33.txt" 2>&1
status=$?
grep -H "RESULT" "$here/results/sk2026_r3_33.txt" | tail -1
exit $status
