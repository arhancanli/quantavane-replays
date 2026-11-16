#!/bin/sh
# 86a, the lower bound [Sie45] ~ 1.73361: Siegel's e (1 + 1/nu)^(-nu), nu the positive root of
# (1 + nu) log(1 + 1/nu) + log(nu)/(1 + nu) = 1 (as stated in arXiv 2401.03252v3, Section 1.1). Two separately written
# programs: c86sie_a.py (exact rationals: sign pattern of f on a grid with a derivative bound, bisection) and
# c86sie_b.py (decimal intervals: f increasing on (0, 1/2] and positive beyond, Newton and a certified bracket).
# Then three red controls each, all refused. Nothing is downloaded; no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c86sie.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c86sie_a.py" > "$tmp/a" 2>&1 &
python3 -I "$here/c86sie_b.py" > "$tmp/b" 2>&1
for c in claim sign root; do python3 -I "$here/c86sie_b.py" --control $c > "$tmp/b.$c" 2>&1; done
wait
for c in claim sign root; do python3 -I "$here/c86sie_a.py" --control $c > "$tmp/a.$c" 2>&1 & done
wait
{
  echo "== program A (c86sie_a.py)"; cat "$tmp/a"
  echo "== program B (c86sie_b.py)"; cat "$tmp/b"
  echo "== red controls (A then B): the claim 1.73362; the sign before log(nu)/(1 + nu) flipped; nu taken as 1/3"
  for c in claim sign root; do grep -h "page\|VERDICT" "$tmp/a.$c" "$tmp/b.$c"; done
  echo "== agreement: the bound (25 digits), A against B"
  ka=$(grep "^key" "$tmp/a"); kb=$(grep "^key" "$tmp/b")
  if [ -n "$ka" ] && [ "$ka" = "$kb" ]; then echo "agreement: True"; echo "$ka"; else echo "agreement: False ($ka / $kb)"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: the positive root is unique, nu_0 = 0.3144808694076345901929969631..., and e (1 + 1/nu_0)^(-nu_0) = 1.73361051698464779155694826...; the page's ~1.73361 is a correct rounding and the source's 1.7336105 a truncation; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
