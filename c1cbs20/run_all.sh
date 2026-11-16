#!/bin/sh
# 1c, the upper bound [BS2020] C_1c <= 1/(2(1 + theta_0)), theta_0 = -min_x sin(x)/x, from its definition. Two separately
# written programs: c1cbs20_a.py (the critical point x_0 of sin x/x by bisection on sin x - x cos x, exact rational
# intervals, theta_0 = 1/sqrt(1 + x_0^2)) and c1cbs20_b.py (branch and bound on sin x/x over [pi/2, 4 pi] with a
# second-order Taylor bound, 80-digit decimals moved outward). Then three red controls each, all refused.
# Nothing is downloaded; no program of any source was read or run. Usage: ./run_all.sh   Python 3.9+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c1cbs20.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c1cbs20_a.py" > "$tmp/a" 2>&1 &
python3 -I "$here/c1cbs20_b.py" > "$tmp/b" 2>&1
for c in theta bound range; do python3 -I "$here/c1cbs20_b.py" --control $c > "$tmp/b.$c" 2>&1; done
wait
for c in theta bound max; do python3 -I "$here/c1cbs20_a.py" --control $c > "$tmp/a.$c" 2>&1 & done
wait
{
  echo "== program A (c1cbs20_a.py)"; cat "$tmp/a"
  echo "== program B (c1cbs20_b.py)"; cat "$tmp/b"
  echo "== red controls, program A: theta_0 claimed as 0.2172337...; bound claimed as 0.4107674; theta_0 from the maximum"
  for c in theta bound max; do grep -h "page\|undefined\|VERDICT" "$tmp/a.$c"; done
  echo "== red controls, program B: theta_0 claimed as 0.2172337...; bound claimed as 0.4107674; search only on [pi/2, pi]"
  for c in theta bound range; do grep -h "page\|VERDICT" "$tmp/b.$c"; done
  echo "== agreement: theta_0 and the bound (30 digits), A against B"
  ka=$(grep "^key" "$tmp/a"); kb=$(grep "^key" "$tmp/b")
  if [ -n "$ka" ] && [ "$ka" = "$kb" ]; then echo "agreement: True"; echo "$ka"; else echo "agreement: False ($ka / $kb)"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: theta_0 = 0.21723362821122165740827932556..., the page's 0.2172336... a truncation; 1/(2(1 + theta_0)) = 0.41076748818940533767047553314..., so the page's 0.4107675... is a correct rounding upward, valid as an upper bound; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
