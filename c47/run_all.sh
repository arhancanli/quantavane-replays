#!/bin/sh
# Replay the C47 lower bounds of [Tom2026] (Theorem 1.3, 50 atoms) and PR #321 (59 atoms) exactly, and PR #341
# (63 atoms) when its restricted construction.json is present; one output file per row in results/, with the red
# controls and the full family of boxes, after the self-tests (Lin's Phi, a single atom, Tomar's bare lattice and
# brute force from the definition).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about two minutes.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c47_check.py" selftest "$root" > "$here/results/selftest.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/selftest.txt"
for row in tom2026 pr321 pr341; do
  python3 -I "$here/c47_check.py" $row "$root" --controls --full > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
exit $status
