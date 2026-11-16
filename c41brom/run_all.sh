#!/bin/sh
# 41b, the lower bound [Romik2016] 1.64495521: the area of Romik's ambidextrous sofa in closed form (arXiv 1606.08111:
# X + arctan Y with x^2 (x + 3) = 8 and x (4x^2 + 3) = 1; eq. (5) in radicals). Two separately written programs:
# c41brom_a.py (the cubics' real roots by bisection in exact rationals, arctan by halving and series) and
# c41brom_b.py (the radical form in 60-digit decimals, arctan by Euler's series). Three red controls each.
# Nothing is downloaded; no program of any source was read or run. Usage: ./run_all.sh   Python 3.9+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c41brom.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c41brom_a.py" > "$tmp/a" 2>&1
python3 -I "$here/c41brom_b.py" > "$tmp/b" 2>&1
for c in claim Y X; do python3 -I "$here/c41brom_a.py" --control $c > "$tmp/a.$c" 2>&1; done
for c in claim radical atan; do python3 -I "$here/c41brom_b.py" --control $c > "$tmp/b.$c" 2>&1; done
{
  echo "== program A (c41brom_a.py)"; cat "$tmp/a"
  echo "== program B (c41brom_b.py)"; cat "$tmp/b"
  echo "== red controls, A: the claim 1.64495522; 4x^3 + 3x + 1; x^2 (x + 3) = 9"
  for c in claim Y X; do grep -h "printed\|VERDICT" "$tmp/a.$c"; done
  echo "== red controls, B: the claim 1.64495522; 3 + 2 sqrt 3 for 3 + 2 sqrt 2; the arctan term dropped"
  for c in claim radical atan; do grep -h "printed\|VERDICT" "$tmp/b.$c"; done
  echo "== agreement: the area (25 digits), A against B"
  ka=$(grep "^key" "$tmp/a"); kb=$(grep "^key" "$tmp/b")
  if [ -n "$ka" ] && [ "$ka" = "$kb" ]; then echo "agreement: True"; echo "$ka"; else echo "agreement: False ($ka / $kb)"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: area = 1.64495521842544085166880934...; the cubic form and the radical form agree; the printed 1.64495521 is a truncation; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
