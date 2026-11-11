#!/bin/sh
# Replay the asterisked 44a upper bound [P2026] from the files fetched by fetch.sh, red controls included; output in
# results/pin.txt. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/c44_check.py" "$root" --controls > "$here/results/pin.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/pin.txt" | tr '\n' ' '; echo
# the two [CHS2020] lower-bound rows (read-once decision trees; the functions Xi_d), d = 1..8, exactly; no download;
# about a minute: results/chs2020.txt
python3 -I "$here/chs2020_check.py" --controls > "$here/results/chs2020.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/chs2020.txt"
exit $status
