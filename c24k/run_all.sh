#!/bin/sh
# Replay of the 24a row [Kun2023] C_24 >= 1 + sqrt 2: Kunisky's recursive matrices A_k, exact in Z[sqrt 2], unit columns
# and exhaustive discrepancies for k <= 4, with red controls. Needs ./fetch.sh first (the paper's LaTeX source).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about ten seconds.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results
python3 -I kunisky_check.py "$root" --controls > results/kun2023.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/kun2023.txt
exit $status
