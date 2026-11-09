#!/bin/sh
# Replay of PR #283's statement (37a): no Boolean function of degree 4 has sensitivity 10, from the definition and the
# certificate files alone (fetched by fetch.sh; read as data), with red controls; one output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; a few seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c37_check.py" "$root/kosh" --controls > "$here/results/kos2026.txt" 2>&1; status=$?
grep -H "RESULT\|run time" "$here/results/kos2026.txt"
exit $status
