#!/bin/sh
# Replay the quantitative hypotheses of the 68a row [Pin2026] (C68 >= 54/125) and of PR #274's two certificates
# (0.446, 0.44605) and of PR #299's contact-row certificate (0.447, row_check.py) from the files fetched by fetch.sh,
# red controls included. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
for row in "pin certificate-jung.json" "vid0446 certificate-c0446.json" "vid044605 certificate-c044605.json"; do
  set -- $row
  python3 -I "$here/c68_check.py" "$root/$1/$2" --controls > "$here/results/$1.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$1.txt"
done
# PR #299 (contact-row model, 0.447): measure, row cover and red controls (a few minutes)
python3 -I "$here/row_check.py" "$root" --controls > "$here/results/pin0447.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/pin0447.txt"
# the [Wik2026] row (0.4263; added with PR #239), from the paper's Table 1 alone (no file); red controls (two minutes)
python3 -I "$here/wik_check.py" --controls > "$here/results/wik2026.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/wik2026.txt"
exit $status
