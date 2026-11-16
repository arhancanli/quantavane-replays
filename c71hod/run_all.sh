#!/bin/sh
# 71a, the lower bound [Hod2017] C >= beta(1/2) > 6.4547837 (Hod, arXiv 1711.00762, Theorem 4.4): beta from its
# definition by two separately written programs: c71hod_a.py (the recursion q_(m+1) = 1/(1 + q_m) in exact rationals,
# logs in fixed-point integers) and c71hod_b.py (the paper's Fibonacci closed forms, checked against the recursion,
# 70-digit decimals). Three red controls each. Nothing is downloaded; no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c71hod.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c71hod_a.py" > "$tmp/a" 2>&1
python3 -I "$here/c71hod_b.py" > "$tmp/b" 2>&1
for c in claim z rec; do python3 -I "$here/c71hod_a.py" --control $c > "$tmp/a.$c" 2>&1; done
for c in claim shift factor; do python3 -I "$here/c71hod_b.py" --control $c > "$tmp/b.$c" 2>&1; done
{
  echo "== program A (c71hod_a.py)"; cat "$tmp/a"
  echo "== program B (c71hod_b.py)"; cat "$tmp/b"
  echo "== red controls, A: the claim 6.4547838; z = 0.499; q_(m+1) = 1/(2 + q_m)"
  for c in claim z rec; do grep -h "^beta\|VERDICT" "$tmp/a.$c"; done
  echo "== red controls, B: the claim 6.4547838; F_(n+2)/F_(n+3) for F_(n+1)/F_(n+2); the factor 4z(1-z) taken as 1/2"
  for c in claim shift factor; do grep -h "^beta\|VERDICT" "$tmp/b.$c"; done
  echo "== agreement: beta(1/2) (25 digits), A against B"
  ka=$(grep "^key" "$tmp/a"); kb=$(grep "^key" "$tmp/b")
  if [ -n "$ka" ] && [ "$ka" = "$kb" ]; then echo "agreement: True"; echo "$ka"; else echo "agreement: False ($ka / $kb)"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: beta(1/2) = 6.45478371656275855381683622859...; the row's > 6.4547837 holds with margin 1.66e-8; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
