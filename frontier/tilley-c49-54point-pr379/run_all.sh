#!/bin/sh
# Tilley, PR #379 (49a): the 54-point block, checked by program A (Python, two-level box method on our 4 November
# checker ../../c49/sunflower_check.py) and program B (JavaScript, written separately), then the red controls.
# Python 3.9+ and node 18+; a few seconds.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
status=0
python3 -I tilley_a.py ../../c49/sunflower_check.py > results/A.txt 2>&1 && echo "A PASS" || { echo "A FAIL"; status=1; }
node tilley_b.mjs > results/B.txt 2>&1 && echo "B PASS" || { echo "B FAIL"; status=1; }
for c in plant-cross plant-contain drop-A1; do
  if python3 -I tilley_a.py ../../c49/sunflower_check.py --control $c > results/control_A_$c.txt 2>&1; then echo "control A $c NOT refused"; status=1; else echo "control A $c refused"; fi
done
for c in plant-sunflower drop-chair swap12; do
  if node tilley_b.mjs $c > results/control_B_$c.txt 2>&1; then echo "control B $c NOT refused"; status=1; else echo "control B $c refused"; fi
done
[ $status = 0 ] && echo "ALL PASS" || echo FAIL
exit $status
