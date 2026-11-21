#!/bin/sh
# EinsteinArena 'Discretized Kakeya Needle (n = 128)': exact union areas of four submissions (Witness, Zlab-QiHuang,
# kgrey, the Station seed) by program A (Python, horizontal slices) and program B (JavaScript BigInt, Green's theorem),
# then red controls. The two exact fractions must be identical. Python 3 and node 18+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
shasum -a 256 -c SHA256SUMS > /dev/null || { echo "FAIL: data/solutions.json differs from SHA256SUMS"; exit 1; }
status=0
for id in 2816 2849 2772 2563; do
  python3 -I kakeya_a.py data/solutions.json $id > results/A_$id.txt 2>&1 || status=1
  node kakeya_b.mjs data/solutions.json $id > results/B_$id.txt 2>&1 || status=1
  fa=$(grep '^exact area' results/A_$id.txt); fb=$(grep '^exact area' results/B_$id.txt)
  if [ -n "$fa" ] && [ "$fa" = "$fb" ] && grep -q 'equals the stated score: True' results/A_$id.txt; then
    echo "$id: A and B give the same exact area; its double equals the stated score"
  else echo "$id: MISMATCH"; status=1; fi
done
# red control: the 64th triangle moved by 1/1000 -> both must give a different area (and the same one)
python3 -I kakeya_a.py data/solutions.json 2816 --control shift > results/control_A_shift.txt 2>&1
node kakeya_b.mjs data/solutions.json 2816 shift > results/control_B_shift.txt 2>&1
if grep -q 'equals the stated score: False' results/control_A_shift.txt && \
   [ "$(grep '^area' results/control_A_shift.txt)" = "$(grep '^area' results/control_B_shift.txt)" ]; then
  echo "control shift: refused by both (area changes, and A and B agree on the new area)"
else echo "control shift: NOT refused"; status=1; fi
[ $status = 0 ] && echo "ALL PASS" || echo FAIL
exit $status
