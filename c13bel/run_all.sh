#!/bin/sh
# 13b, the lower bound [Elek1994] 0.8257 (a disk of diameter 1 and an equilateral triangle of side 1 in one convex set):
# two separately written programs. c13bel_a.py: the argument (caps, convexity of L, the 120-degree symmetry, the two
# cases) with its numbers in exact rational intervals, min = pi/8 + sqrt(3)/4 = 0.82571178359...; c13bel_b.py: branch
# and bound over the triangle's position in 40-digit decimals, no convexity used, every cell certified >= 0.8257.
# Then three red controls each, all refused. Nothing is downloaded; no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; about two minutes on one core.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c13bel.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c13bel_a.py" > "$tmp/a" 2>&1
for c in claim d1 side; do python3 -I "$here/c13bel_a.py" --control $c > "$tmp/a.$c" 2>&1; done
python3 -I "$here/c13bel_b.py" > "$tmp/b" 2>&1
for c in claim radius no-ii; do python3 -I "$here/c13bel_b.py" --control $c > "$tmp/b.$c" 2>&1; done
{
  echo "== program A (c13bel_a.py)"; cat "$tmp/a"
  echo "== program B (c13bel_b.py)"; cat "$tmp/b"
  echo "== red controls, program A: the claim 0.8258; d1 = 3/4; a triangle of side 0.9"
  for c in claim d1 side; do grep -h "printed lower\|case (a)\|VERDICT" "$tmp/a.$c"; done
  echo "== red controls, program B: the claim 0.8258; a disk of radius 0.49; bound (ii) removed"
  for c in claim radius no-ii; do grep -h "CERTIFIED\|certified\|VERDICT" "$tmp/b.$c"; done
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (A: min over placements = pi/8 + sqrt(3)/4 = 0.82571178359..., attained with the centres together; B: every placement certified >= 0.8257 by branch and bound; the printed 0.8257 is a truncation; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
