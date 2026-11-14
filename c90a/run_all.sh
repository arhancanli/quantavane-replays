#!/bin/sh
# 90a: the page's elementary covering system 0 mod 2, 0 mod 3, 1 mod 4, 1 mod 6, 11 mod 12 (lower bound C_90 >= 2),
# by two separately written programs (A: every residue of the period; B: inclusion-exclusion density in exact
# rationals), with four red controls that both must refuse. No download. Python 3.9+, standard library only; seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/elementary.txt"
mkdir -p "$here/results"
both() { echo "== $1"; python3 -I "$here/c90a_cover.py" "$2"; python3 -I "$here/c90a_density.py" "$2"; }
{
  both "the page's system" "0:2 0:3 1:4 1:6 11:12"
  both "red control: 11 mod 12 replaced by 10 mod 12" "0:2 0:3 1:4 1:6 10:12"
  both "red control: 1 mod 6 replaced by 5 mod 6" "0:2 0:3 1:4 5:6 11:12"
  both "red control: the class 1 mod 6 removed" "0:2 0:3 1:4 11:12"
  both "red control: a covering with a repeated modulus (0 mod 2, 1 mod 2)" "0:2 1:2"
} > "$out" 2>&1
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 8 ]; then echo "RESULT: ALL CHECKS PASS (the system accepted by both programs; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0; fi
echo "RESULT: FAIL (pass $pass, refused $ref)" >> "$out"; tail -1 "$out"; exit 1
