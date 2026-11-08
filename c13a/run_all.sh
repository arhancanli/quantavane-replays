#!/bin/sh
# Replay of the 13a upper-bound row [K2026] (PR #314): Keller's quadrilateral K from the side data in the paper's
# LaTeX source, exact rational arithmetic, with red controls; one output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a second.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/quad_check.py" "$root" --controls > "$here/results/k2026.txt" 2>&1; status=$?
grep -H "RESULT\|controls all" "$here/results/k2026.txt"
exit $status
