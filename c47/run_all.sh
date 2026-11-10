#!/bin/sh
# Replay the C47 lower bounds of [Lin2026] (PR #192, Lin's lattice), [Tom2026] (Theorem 1.3, 50 atoms) and PR #321
# (59 atoms) exactly, and PR #341 (63 atoms) when its restricted construction.json is present; one output file per row
# in results/, with the red controls (and for the last three the full family of boxes), after the self-tests (Lin's
# Phi, a single atom, Tomar's bare lattice and brute force from the definition).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about two minutes.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c47_check.py" selftest "$root" > "$here/results/selftest.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/selftest.txt"
# the lower-bound row [Lin2026] (PR #192), from its closed forms (no file); red controls
python3 -I "$here/c47_check.py" lin2026 --controls > "$here/results/lin2026.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/lin2026.txt"
for row in tom2026 pr321 pr341 pr350 pr350r8; do
  python3 -I "$here/c47_check.py" $row "$root" --controls --full > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
exit $status
