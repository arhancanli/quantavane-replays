#!/bin/sh
# Replay the three explicit 71a lower bounds ([MI2026] n = 14, [MI2026b] n = 17, [Num2026] n = 18) from the data fetched
# by fetch.sh; one output file per row in results/, red controls included.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
for row in mi2026 mi2026b num2026; do
  python3 -I "$here/c71_check.py" $row "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT" "$here/results/$row.txt"
done
exit $status
