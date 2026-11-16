#!/bin/sh
# 5a, the upper bound [COT26] C_5a <= 0.94304 (Carter, O'Bryant and Tao, arXiv 2610.36207v1): the WT-triple witnesses
# (dcartermath/sidon-error @ ff5b032, asymptotic/ for (10, 0.944814, 0.941259) and explicit/ for (3, 1.630563, 0.547610))
# checked against Definition 2.2 by two separately written exact programs: c5acot_a.py (the convolution at the knots by
# the index formula) and c5acot_b.py (a generic step-function engine, with linearity between knots checked at every
# midpoint). Then three red controls each. Needs ./fetch.sh first (data only; no program of the source read or run).
# Usage: ./run_all.sh [CERTS_ROOT]   Python 3.9+; about half a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/c5acot.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c5acot_a.py" "$root" > "$tmp/a" 2>&1 &
python3 -I "$here/c5acot_b.py" "$root" > "$tmp/b" 2>&1
wait
for c in claim middle mirror; do python3 -I "$here/c5acot_a.py" "$root" --control $c > "$tmp/a.$c" 2>&1; done
for c in claim prob height; do python3 -I "$here/c5acot_b.py" "$root" --control $c > "$tmp/b.$c" 2>&1; done
{
  echo "== program A (c5acot_a.py)"; cat "$tmp/a"
  echo "== program B (c5acot_b.py)"; cat "$tmp/b"
  echo "== red controls, A: the claim 0.94303; a middle-third table height set to 0.999; the mirroring ignored"
  for c in claim middle mirror; do grep -h "(i)\|(vi)\|the row\|VERDICT" "$tmp/a.$c"; done
  echo "== red controls, B: the claim 0.94303; 1/1000 of probability moved from atom 1 to atom 0; a window height raised"
  for c in claim prob height; do grep -h "shapes\|(vi)\|the row\|VERDICT" "$tmp/b.$c"; done
  echo "== agreement: alpha and beta of both witnesses (20 digits), A against B"
  ka=$(grep "key" "$tmp/a" | sed 's/^ *//'); kb=$(grep "key" "$tmp/b" | sed 's/^ *//')
  if [ -n "$ka" ] && [ "$ka" = "$kb" ]; then echo "agreement: True"; echo "$ka"; else echo "agreement: False ($ka / $kb)"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both witnesses satisfy every condition of Definition 2.2 exactly; (10, alpha, beta) with alpha = 0.94481304378915..., beta = 0.94125802624178... (the fraction printed in Appendix A.3), 2 sqrt(alpha beta) = 1.88606771963... <= 1.88607, threshold 482.34 <= 483, sqrt(alpha beta) = 0.94303385981... <= 0.94304; (3, 1.63056262452..., 0.54760956513...) gives 1.88988009115... <= 1.8898801 and 29.9987 < 30; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
