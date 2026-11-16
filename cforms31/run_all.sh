#!/bin/sh
# Closed forms on eight pages that had no replay folder: 6a, 11b, 35a (irrational closed forms) and 88b, 23b, 46a, 56a,
# 61a (rationals), each against the decimal the page prints, by two separately written programs (cforms_a.py: integer
# fixed point at 2^-200 with its own sqrt and log; cforms_b.py: decimal intervals and exact integer comparisons),
# then four red controls that both must refuse. Nothing is downloaded; no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+, standard library only; a few seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/cforms31.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/cforms_a.py" > "$tmp/m.a" 2>&1
python3 -I "$here/cforms_b.py" > "$tmp/m.b" 2>&1
for c in 6a-plus 11b-print 46a-15 61a-sign; do
  python3 -I "$here/cforms_a.py" --control $c > "$tmp/$c.a" 2>&1
  python3 -I "$here/cforms_b.py" --control $c > "$tmp/$c.b" 2>&1
done
vals() { grep -h "^row " "$tmp/m.$1" | sed 's/:.*= \([0-9]*\.[0-9]\{30\}\).*/ \1/; s/^\(row 61a [^:]*\):.*/\1/'; }
{
  echo "== program A (cforms_a.py)"; cat "$tmp/m.a"
  echo "== program B (cforms_b.py)"; cat "$tmp/m.b"
  echo "== red control 1: (3 + sqrt 5)/2 in place of (3 - sqrt 5)/2"; grep -h "CONTROL\|MISMATCH\|VERDICT" "$tmp/6a-plus.a" "$tmp/6a-plus.b"
  echo "== red control 2: log2(3/2) printed as 0.58497"; grep -h "CONTROL\|MISMATCH\|VERDICT" "$tmp/11b-print.a" "$tmp/11b-print.b"
  echo "== red control 3: 4 - 2/13 in place of 4 - 2/15"; grep -h "CONTROL\|MISMATCH\|VERDICT" "$tmp/46a-15.a" "$tmp/46a-15.b"
  echo "== red control 4: lambda_1 = 1/4 + theta^2 in place of 1/4 - theta^2"; grep -h "CONTROL\|MISMATCH\|VERDICT" "$tmp/61a-sign.a" "$tmp/61a-sign.b"
  echo "== agreement: every row and its first 30 digits, A against B"
  if [ -n "$(vals a)" ] && [ "$(vals a)" = "$(vals b)" ]; then echo "agreement: True ($(vals a | wc -l | tr -d ' ') rows)"; else echo "agreement: False"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED\|SELF-TEST FAILED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 8 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (18 rows on 8 pages, both programs, the same 30 digits; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
