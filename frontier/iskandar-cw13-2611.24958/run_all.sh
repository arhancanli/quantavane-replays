#!/bin/sh
# The thirteen constant-weight codes of arXiv 2611.24958 (P. R. A. Iskandar): two separately written exact programs.
# A = cw_a.py (no (w - d/2 + 1)-subset in two words, by hashing; orbit rebuild as sets); B = cw_b.c (64-bit bitsets:
# pairwise intersections inside coordinate-pair buckets, plus all pairs directly for codes of at most 30,000 words;
# orbit rebuild by bit permutation, sorted). Both: length, weight, no repeated word, size, minimum distance >= d, and
# the full list equals the union of the orbits of the representatives. Then three red controls each.
# No program of the source was read or run. Usage: ./run_all.sh [CERTS_ROOT] (default ./certs). Python 3.9+ and a C
# compiler; under a minute, one process.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
anc="$root/anc"
out="$here/results/cw13.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
${CC:-cc} -O2 -o "$tmp/b" "$here/cw_b.c" || { echo "RESULT: FAIL (compile)"; exit 1; }
python3 -I "$here/cw_a.py" "$anc" > "$tmp/a" 2>&1
"$tmp/b" "$anc" > "$tmp/bo" 2>&1
{
  echo "== program A (cw_a.py)"; cat "$tmp/a"
  echo "== program B (cw_b.c)"; cat "$tmp/bo"
  n=0
  for c in distance repeat weight; do
    n=$((n + 1))
    case $c in
      distance) echo "== red control $n: in A(16,4,7) word 2 replaced by word 1 with one 1 moved (distance 2)";;
      repeat) echo "== red control $n: in A(16,4,7) word 2 replaced by a copy of word 1";;
      weight) echo "== red control $n: in A(16,4,7) word 2 replaced by word 1 with one extra 1 (weight 8)";;
    esac
    python3 -I "$here/cw_a.py" "$anc" --control $c 2>&1 | grep "CONTROL\|FAIL\|VERDICT"
    "$tmp/b" "$anc" $c 2>&1 | grep "CONTROL\|FAIL\|VERDICT"
  done
} > "$out" 2>&1
rm -rf "$tmp"
ok_a=$(grep -c "^A: A(.*clashing 0, .*equals list: yes" "$out"); ok_b=$(grep -c "^B: A(.*repeats 0, weight ok yes, .*equals list: yes" "$out")
tot_a=$(grep -c "^A: total words 408729, total representatives 11386" "$out"); tot_b=$(grep -c "^B: codes 13, total words 408729, total representatives 11386" "$out")
pass=$(grep -c "^[AB]: VERDICT PASS" "$out"); ref=$(grep -c "^[AB]: VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback" "$out")
if [ "$ok_a" = 13 ] && [ "$ok_b" = 13 ] && [ "$tot_a" = 1 ] && [ "$tot_b" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (A and B: all 13 codes valid with the claimed sizes, 408,729 words, minimum distance >= d, each equal to the union of the orbits of its 11,386 representatives in total; every red control refused)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (A ok $ok_a, B ok $ok_b, pass $pass, refused $ref)" >> "$out"; tail -1 "$out"; exit 1
