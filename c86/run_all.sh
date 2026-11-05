#!/bin/sh
# Replay the asterisked 86a upper bounds (and the open PRs #270, #298) from the files fetched by fetch.sh; one output
# file per row in results/, red controls included. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
for row in kos pin vid p298 vid2 vid3; do      # [Kos2026], [Pin26], PR #270, PR #298: every condition, red controls included
  python3 -I "$here/c86_check.py" "$row" "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
exit $status
