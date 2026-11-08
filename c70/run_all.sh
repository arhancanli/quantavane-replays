#!/bin/sh
# Replay the asterisked 70a lower bound [Pin2026] from the file fetched by fetch.sh, and the [Pat2026] row 1.201927 (the
# full Steiner sum at n = 300000, a/r = 13/1000, with the same interval primitives), red controls included; output in
# results/pin.txt and results/pat2026.txt. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c70_check.py" "$root" --controls > "$here/results/pin.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/pin.txt" | tr '\n' ' '; echo
python3 -I "$here/pat_check.py" "$here/c70_check.py" "$root/pin/certificates/c70_candidate.json" --controls \
    > "$here/results/pat2026.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/pat2026.txt"
exit $status
