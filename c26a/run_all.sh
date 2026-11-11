#!/bin/sh
# Replay #284 (26a): the six certificates and our own f_d from WRITEUP.md Section 3, with red controls;
# output in results/koshta_c26a.txt. Usage: ./run_all.sh [CERTS_ROOT]   (Python 3.9+, standard library only; seconds)
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/bh_check.py" "$root" --controls > "$here/results/koshta_c26a.txt" 2>&1
status=$?
grep -H "RESULT\|refused\|term by term" "$here/results/koshta_c26a.txt"
# the [ADGP2025] row (the degree-d address function, d = 1..5, exactly; no download): results/adgp2025_address.txt
python3 -I "$here/address_check.py" --controls > "$here/results/adgp2025_address.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/adgp2025_address.txt"
exit $status
