#!/bin/sh
# Replay the quantitative hypotheses of the 68a row [Pin2026] (C68 >= 54/125) from the file fetched by fetch.sh, red
# controls included. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c68_check.py" "$root/pin/certificate-jung.json" --controls > "$here/results/pin.txt" 2>&1
status=$?
grep -H "RESULT\|run time" "$here/results/pin.txt"
exit $status
