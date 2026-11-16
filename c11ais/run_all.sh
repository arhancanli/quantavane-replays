#!/bin/sh
# 11a, the upper bound [IS2024] C_11a <= pi/2 - delta, delta ~ 0.00013: delta from Ivanisvili and Stone's integral
# (arXiv 2407.04835, Section 2.4) by two separately written exact routes in the variable s = arcsin(2p - 1):
# c11ais_a.py (reduction formulas for int cos^n, exact rationals in the basis 1, sqrt 2, s*, pi) and c11ais_b.py
# (cosine-polynomial expansion and Chebyshev values of sin(2j s*), 60-digit decimals). Three red controls each.
# Nothing is downloaded; no program of any source was read or run. Usage: ./run_all.sh   Python 3.9+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c11ais.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c11ais_a.py" > "$tmp/a" 2>&1
python3 -I "$here/c11ais_b.py" > "$tmp/b" 2>&1
for c in den cross claim; do python3 -I "$here/c11ais_a.py" --control $c > "$tmp/a.$c" 2>&1; python3 -I "$here/c11ais_b.py" --control $c > "$tmp/b.$c" 2>&1; done
{
  echo "== program A (c11ais_a.py)"; cat "$tmp/a"
  echo "== program B (c11ais_b.py)"; cat "$tmp/b"
  echo "== red controls (A then B): 44 for 45 in the denominator; the crossover at sin^2 = 1/4; the claim delta ~ 0.000133"
  for c in den cross claim; do grep -h "^delta in\|row's\|VERDICT" "$tmp/a.$c" "$tmp/b.$c"; done
  echo "== agreement: delta (20 digits), A against B"
  ka=$(grep "^key" "$tmp/a"); kb=$(grep "^key" "$tmp/b")
  if [ -n "$ka" ] && [ "$ka" = "$kb" ]; then echo "agreement: True"; echo "$ka"; else echo "agreement: False ($ka / $kb)"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: delta = 0.00013279130299166442911386578...; pi/2 - delta = 1.5706635354919049548...; the row's delta ~ 0.00013 is a correct rounding; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
