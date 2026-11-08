#!/bin/sh
# Replays of the two latest 7a rows (PR #214): [B2026] mu(pi) < 7.101862832357 (Bai's sigma = r + C and tau = -s - C:
# the periodic set E, omega by two routes, C; zeta and r; y_c and s) and [ZZ2020] 7.103205334137... (Zeilberger and
# Zudilin's Part II: the prime saving psi(2/3) - psi(1/2) by two routes, the indicial cubic and its saddle-point
# route, the two rates), interval arithmetic of our own, with controls. One output file per row in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 20 seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/pi_check.py" "$root" --controls > "$here/results/b2026.txt" 2>&1 || status=1
python3 -I "$here/pi_check.py" "$root" --controls --row=zz2020 > "$here/results/zz2020.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/b2026.txt" "$here/results/zz2020.txt"
exit $status
