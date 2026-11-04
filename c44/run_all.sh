#!/bin/sh
# Replay the asterisked 44a upper bound [P2026] from the files fetched by fetch.sh, red controls included; output in
# results/pin.txt. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c44_check.py" "$root" --controls > "$here/results/pin.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/pin.txt" | tr '\n' ' '; echo
exit $status
