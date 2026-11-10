#!/bin/sh
# Replay the contour inequality behind C_1d >= 1.152* (#270) by two methods, with red controls;
# output in results/vid_c1d_lower.txt. Usage: ./run_all.sh   (Python 3.9+, standard library only; about 7 minutes)
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/results"
python3 -I "$here/c1d_lower_check.py" --controls > "$here/results/vid_c1d_lower.txt" 2>&1
status=$?
grep -H "RESULT\|refused\|HOLDS" "$here/results/vid_c1d_lower.txt"
exit $status
