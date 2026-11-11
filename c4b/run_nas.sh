#!/bin/sh
# Replays of the 4b rows [Nas2026a] and [Jon2026b] (PR #187) and [Nas2026b] (PR #219): every finite hypothesis of the
# interval-moment criterion (Naslund's Theorem 3.3; Jones's NK.powerLowerBound_of_intervalMoments), from the papers'
# TeX and witness files and from Jones's data files, with controls. One output file per row.
# Usage: ./run_nas.sh [CERTS_ROOT]   (default: ./certs; run ./fetch_nas.sh and ./fetch_jon.sh first). Python 3.9+,
# standard library only; about 7 minutes and under 300 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
for row in nas2026a nas2026b jon2026b; do
  python3 -I "$here/naslund_check.py" $row "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|controls all" "$here/results/$row.txt"
done
exit $status
