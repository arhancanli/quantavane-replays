#!/bin/sh
# Replay both [F2026] certificates (the record 1.82634836... for groups without 2-torsion, the 3b row of PR #329, and
# the all-groups 1.82641618...) with their red controls, after the self-tests on the paper's two standard-form
# certificates (2762/1509 and 152/83, a chain configuration). Output in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c3b_check.py" selftest "$root" > "$here/results/selftest.txt" 2>&1 || status=1
python3 -I "$here/c3b_check.py" both "$root" --controls > "$here/results/c3b_check.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/selftest.txt" "$here/results/c3b_check.txt"
exit $status
