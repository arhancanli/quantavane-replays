#!/bin/sh
# Replay the 24a row C24 <= 6.9013* (Akbas and Sra, Theorem 2.4 of arXiv 2609.27172v2; PR #250) from the equation in
# the paper's source, two routes, with red controls; one output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a second.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results
python3 -I komlos_check.py "$root/as2026" --controls > results/as2026.txt 2>&1; status=$?
grep -H "RESULT\|enclosure C_cos\|bound C_cos\|C_cos in \[" results/as2026.txt
exit $status
