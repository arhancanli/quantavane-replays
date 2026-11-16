#!/bin/sh
# 8a: the lower bound C_8 >= 2/log(gamma_1) ~ 0.755106. Two separately written programs prove that zeta has a zero
# 1/2 + i gamma with 14.13472514173469378 < gamma < 14.13472514173469380 (a sign change of Hardy's Z, each with Im of
# e^{i theta} zeta enclosing 0 as a self-test), and enclose 2/log gamma: c8a_a.py (Euler-Maclaurin, Stirling at shift
# 20, integer intervals at 2^-320), c8a_b.py (Cohen-Rodriguez Villegas-Zagier acceleration of eta with its error bound,
# Stirling at shift 30, decimal intervals). Then two red controls that both must refuse.
# Nothing is downloaded; no program of any source was read or run. Usage: ./run_all.sh   Python 3.9+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c8a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
for n in main theta-sign bracket; do
  arg=""; [ "$n" != main ] && arg="--control $n"
  python3 -I "$here/c8a_a.py" $arg > "$tmp/$n.a" 2>&1
  python3 -I "$here/c8a_b.py" $arg > "$tmp/$n.b" 2>&1
done
zs() { grep -h "^Z(" "$1" | grep -o "in \[[^]]*\]"; grep -h "^row 8a" "$1" | grep -o "= [0-9]\.[0-9]\{20\}"; }
{
  echo "== program A (c8a_a.py)"; cat "$tmp/main.a"
  echo "== program B (c8a_b.py)"; cat "$tmp/main.b"
  echo "== red control 1: theta with the wrong sign"; grep -h "CONTROL\|^Z(\|sign change\|VERDICT" "$tmp/theta-sign.a" "$tmp/theta-sign.b"
  echo "== red control 2: a bracket above gamma_1, (14.13472514173469381, 14.13472514173469383)"; grep -h "CONTROL\|^Z(\|sign change\|VERDICT" "$tmp/bracket.a" "$tmp/bracket.b"
  echo "== agreement: Z at both ends (7 digits) and 2/log gamma_1 (20 digits), A against B"
  if [ -n "$(zs "$tmp/main.a")" ] && [ "$(zs "$tmp/main.a")" = "$(zs "$tmp/main.b")" ]; then echo "agreement: True"; zs "$tmp/main.a"; else echo "agreement: False"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 4 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: a zero of zeta at 1/2 + i gamma, 14.13472514173469378 < gamma < 14.13472514173469380; 2/log gamma = 0.75510606143218443826... >= 0.755106; two red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
