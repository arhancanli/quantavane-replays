#!/bin/sh
# Replay the four entropy lower bounds of the 3c page ([G2026], [MI2026], [L2026], [J2026]) from their laws, fetched by
# fetch.sh as data (the [G2026] table is in c3c_check.py); one output file per row in results/, red controls included.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; a few seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
for row in g2026 mi2026 l2026 j2026; do
  python3 -I "$here/c3c_check.py" $row "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT" "$here/results/$row.txt"
done
exit $status
