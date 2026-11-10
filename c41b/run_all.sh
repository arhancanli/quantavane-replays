#!/bin/sh
# Check the 41b bound C41b <= 2 sqrt 2 - 1* (#358) with red controls; output in results/vid_c41b.txt.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; under a minute)
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/results"
python3 -I "$here/sofa_check.py" --controls > "$here/results/vid_c41b.txt" 2>&1
status=$?
grep -H "RESULT\|fails\|True\|exceptions" "$here/results/vid_c41b.txt" | head -20
exit $status
