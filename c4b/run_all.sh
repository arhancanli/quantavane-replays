#!/bin/sh
# Replays of two 4b rows added in PR #187: [Kra2026] C_4b >= 0.752796455874514... (the ten Paley chains, alpha_*) and
# [Jon2026a] C_4b >= 0.753741541837329405... (the two square-DAG blocks on 235 and 299, the eleven-block alpha_inf),
# in interval arithmetic, with the constructions checked on small instances by brute force; controls.
# One output file per row in results/. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/chain_check.py" "$root" --controls > "$here/results/kra2026.txt" 2>&1 || status=1
python3 -I "$here/chain_check.py" "$root" --controls --row=jon2026a > "$here/results/jon2026a.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/kra2026.txt" "$here/results/jon2026a.txt"
exit $status
