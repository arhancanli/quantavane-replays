#!/bin/sh
# 7a, the upper bound [ZZ2020] mu(pi) <= 7.103205334137...: Zeilberger and Zudilin's final formula recomputed by two
# separately written programs: c7azz_a.py (exact rational intervals: the real root of their cubic (18), |N_1|^2 from
# the product of the roots, logs by atanh series) and c7azz_b.py (60-digit decimals: the saddle points y_j of
# g(y) = y (y^2 + 6y + 25)^2/(y - 25)^3 and N_j = g(y_j), checked against the cubic). Three red controls each.
# Nothing is downloaded; no program of any source was read or run. Usage: ./run_all.sh   Python 3.9+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c7azz.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c7azz_a.py" > "$tmp/a" 2>&1
python3 -I "$here/c7azz_b.py" > "$tmp/b" 2>&1
for c in claim const log2; do python3 -I "$here/c7azz_a.py" --control $c > "$tmp/a.$c" 2>&1; done
for c in claim g log2; do python3 -I "$here/c7azz_b.py" --control $c > "$tmp/b.$c" 2>&1; done
{
  echo "== program A (c7azz_a.py)"; cat "$tmp/a"
  echo "== program B (c7azz_b.py)"; cat "$tmp/b"
  echo "== red controls, A: the claim 7.103205334136; 2047 for 2048 in the cubic; 2 log 2 for (5/2) log 2"
  for c in claim const log2; do grep -h "^row\|VERDICT" "$tmp/a.$c"; done
  echo "== red controls, B: the claim 7.103205334136; (y - 25)^2 for (y - 25)^3 in g; 2 log 2 for (5/2) log 2"
  for c in claim g log2; do grep -h "^row\|on the\|satisfy\|VERDICT" "$tmp/b.$c"; done
  echo "== agreement: 1 + sigma/tau (25 digits), A against B"
  ka=$(grep "^key" "$tmp/a"); kb=$(grep "^key" "$tmp/b")
  if [ -n "$ka" ] && [ "$ka" = "$kb" ]; then echo "agreement: True"; echo "$ka"; else echo "agreement: False ($ka / $kb)"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: mu(pi) <= 1 + sigma/tau = 7.1032053341370017275057734...; the row's 7.103205334137 is its truncation, 1.7e-15 below it; three red controls refused by each)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
