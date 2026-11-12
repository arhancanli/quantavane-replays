#!/bin/sh
# Replay the C51 lower bound of [R2026b] (PR #184) with its red controls, and the two parameter pairs it is compared
# with ([S2026]'s and He and Tang's), after the self-tests; one output file per row in results/. Also the [CH1964] row
# B >= 4/7 (ch1964_check.py: the Clunie-Hayman series at K = 6, exact Sturm certificate; no download).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
for row in selftest r2026b s2026 ht2026; do
  python3 -I "$here/c51_check.py" $row "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
python3 -I "$here/ch1964_check.py" --controls > "$here/results/ch1964.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/ch1964.txt" | tr '\n' ' '; echo
exit $status
