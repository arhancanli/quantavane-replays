#!/bin/sh
# Erdos problem 838, upper-bound side of Rohan Arun's minima (@ c301008): two separately written exact programs count, for
# every published integer-coordinate witness, the subsets of size >= 3 in convex position and check that no three points
# are collinear. A = erdos838_a.py (no point inside a triangle of three others); B = erdos838_b.c (convex hull of every
# subset by Andrew's monotone chain, 128-bit orientation). Then three red controls each. No program of the source was
# read or run. The lower bounds (order-type enumeration) are not checked here.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+ and a C compiler; a few seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
c="$root/erdos838"
out="$here/results/erdos838.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
${CC:-cc} -O2 -o "$tmp/b" "$here/erdos838_b.c" || { echo "RESULT: FAIL (compile)"; exit 1; }
D14=574/574/574/575/575/575/576/578/578/580/581/582/586/588
set -- "$c/witness10.txt=189" "$c/witness11.txt=283" "$c/witness12.txt=408" \
  "$c/real13.jsonl:jsonl:0=570" "$c/real13.jsonl:jsonl:1=570" "$c/real13.jsonl:jsonl:2=570" \
  "$c/real13.jsonl:jsonl:3=570" "$c/real13.jsonl:jsonl:4=570" "$c/Witness13.lean:lean=570,286,241,41,2" \
  "$c/witness780-small.txt=780,364,350,64,2@$D14" "$c/witness780.txt=780,364,350,64,2@$D14"
python3 -I "$here/erdos838_a.py" "$@" > "$tmp/a" 2>&1
"$tmp/b" "$@" > "$tmp/bo" 2>&1
W12="$c/witness12.txt=408"; W14="$c/witness780-small.txt=780,364,350,64,2@$D14"
{
  echo "== program A (erdos838_a.py)"; cat "$tmp/a"
  echo "== program B (erdos838_b.c)"; cat "$tmp/bo"
  n=0
  for ctl in collinear move; do for w in "$W12" "$W14"; do
    n=$((n + 1))
    echo "== red control $n: $ctl on $(basename "${w%%=*}")"
    python3 -I "$here/erdos838_a.py" --control $ctl "$w" 2>&1 | grep "CONTROL\|FAIL\|VERDICT"
    "$tmp/b" --control $ctl "$w" 2>&1 | grep "CONTROL\|FAIL\|VERDICT" | grep -v "deletion totals differ at rank [1-9]"
  done; done
  echo "== red control 5: claimed total 407 for witness12"
  python3 -I "$here/erdos838_a.py" "$c/witness12.txt=407" 2>&1 | grep "FAIL\|VERDICT"
  "$tmp/b" "$c/witness12.txt=407" 2>&1 | grep "FAIL\|VERDICT"
} > "$out" 2>&1
sed -n 's/^A: //p' "$tmp/a" | sed 's/, sign string matches//' | grep -v VERDICT > "$tmp/a.n"
sed -n 's/^B: //p' "$tmp/bo" | grep -v VERDICT > "$tmp/b.n"
same=0; cmp -s "$tmp/a.n" "$tmp/b.n" && [ -s "$tmp/a.n" ] && same=1
signs=$(grep -c "sign string matches" "$tmp/a")
rm -rf "$tmp"
pass=$(grep -c "^[AB]: VERDICT PASS" "$out"); ref=$(grep -c "^[AB]: VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback" "$out")
if [ "$same" = 1 ] && [ "$signs" = 5 ] && [ "$pass" = 2 ] && [ "$ref" = 10 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (A and B agree witness by witness: no three collinear; 189, 283, 408, 570 (six 13-point witnesses) and 780 = 364 + 350 + 64 + 2 convex polygons; deleting any point of the 14-point witness leaves at least 574; every red control refused)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (A and B agree: $same, sign strings $signs, pass $pass, refused $ref)" >> "$out"; tail -1 "$out"; exit 1
