#!/bin/sh
# 33a: the finite inputs of [DM2013] A(2) >= 0.316999... (Duursma and Mak, arXiv 1102.4127v4, Proposition 5.5 and
# Theorem 1.3, criterion Theorem 1.1), by two separately written programs (c33a_a.py: point counts over F_(2^d), d <= 10, by a character
# sum and Moebius inversion, conductors from Laurent series in x - x0; c33a_b.py: closed points and Frobenius
# vectors, conductors from Laurent series in A, and the genus again from the zeta functions), then red controls that
# both must refuse. Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; under a minute, under 100 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c33a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c33a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c33a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m
pair x1 --claim 0.317
pair x2 --variant h-x3
pair x3 --variant A-x5x3
pair x4 --S "5:2,6:17,8:15,10:4"
{
  echo "== 33a [DM2013]: program A (c33a_a.py)"
  cat "$tmp/m.a"
  echo "== program B (c33a_b.py)"
  cat "$tmp/m.b"
  for x in x1 x2 x3 x4; do
    case $x in
      x1) echo "== red control 1: the bound claimed as A(2) >= 0.317";;
      x2) echo "== red control 2: y^2 replaced by x^3 in the span of h";;
      x3) echo "== red control 3: A = x^5 + x^3 + 1 (also irreducible) in place of x^5 + x^2 + 1";;
      x4) echo "== red control 4: S with 17 places of degree 6";;
    esac
    grep -h "CONTROL\|a_d(k)\|S asks\|S needs\|MISMATCH\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
  echo "== agreement"
  ka=$(grep "^a_d(k)" "$tmp/m.a" | sed 's/.*://'); kb=$(grep "^a_d(k)" "$tmp/m.b" | sed 's/.*://')
  ha=$(grep "^a_d(H)" "$tmp/m.a" | sed 's/.*://' | tr -d ' []'); hb=$(grep "^closed points of H by degree" "$tmp/m.b" | sed 's/.*://' | tr -d ' []')
  ba=$(grep -o "24576/77527 = 0.31699923897480877629" "$tmp/m.a" | head -1); bb=$(grep -o "24576/77527 = 0.31699923897480877629" "$tmp/m.b" | head -1)
  echo "a_d(k): A$ka  B$kb"
  echo "a_d(H): A $ha  B $hb"
  echo "bound:  A $ba  B $bb"
  if [ "$ka" = "$kb" ] && [ "$ha" = "$hb" ] && [ -n "$ba" ] && [ "$ba" = "$bb" ]; then echo "pair exit 0"; else echo "pair exit 1"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
agree=$(grep -c "^pair exit 0" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 8 ] && [ "$bad" = 0 ] && [ "$agree" = 1 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: a_d(H) = (6,0,1,1,6), a_d(k) = (192,0,0,0,2,16,0,16,0,64), g(k) = 343, inequality (1.1) value 57, A(2) >= 24576/77527 = 0.316999238974...; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
