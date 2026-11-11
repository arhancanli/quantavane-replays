#!/bin/sh
# Replay of the 2a row [C2007]: p(z) = z and A = (0 2; 0 0) from the page; ||A|| = 2 and W(A) the closed unit disk,
# exactly; red controls. Needs ./fetch.sh first (the page at main 87b4ebb).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a second.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results
python3 -I crouzeix_check.py "$root" --controls > results/c2007.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/c2007.txt
exit $status
