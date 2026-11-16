#!/bin/sh
# 8a, pull request #373: the finite inputs of the deduction C_8 <= 2/log(gamma_1) ~ 0.755107*. Two separately written
# programs prove (1) a zero of zeta at 1/2 + i gamma_1 with 14.13472514173469378 < gamma_1 < 14.13472514173469380 (as in
# folder c8a), (2) N(89/5) = 1, so this zero is the lowest (Re zeta(sigma + 17.8 i) > 0 for 1/2 <= sigma <= 2, hence
# |S(17.8)| < 1/2, and theta(17.8)/pi + 1 = 0.99243...), and (3) 2/log gamma_1 <= 0.755107, gamma_1^4 < 4 * 10^4 and
# 8/log(3 * 10^12) < 0.279: c8a373_a.py (Euler-Maclaurin with sigma as an interval, integer intervals at 2^-320),
# c8a373_b.py (Cohen-Rodriguez Villegas-Zagier for eta with sigma as an interval, decimal intervals). Then three red
# controls that both must refuse. [OAI2026] and [PT2021] are not checked.
# Nothing is downloaded; no program of any source was read or run. Usage: ./run_all.sh   Python 3.9+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c8a373.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
for n in main height22 printed bracket; do
  arg=""; [ "$n" != main ] && arg="--control $n"
  python3 -I "$here/c8a373_a.py" $arg > "$tmp/$n.a" 2>&1
  python3 -I "$here/c8a373_b.py" $arg > "$tmp/$n.b" 2>&1
done
key() { grep -h "^(3) 2/log" "$1" | grep -o "\[0\.[0-9]\{22\}, 0\.[0-9]\{22\}\]"; grep -h "integers within 1/2" "$1" | grep -o "1/2: \[[0-9, ]*\]"; }
{
  echo "== program A (c8a373_a.py)"; cat "$tmp/main.a"
  echo "== program B (c8a373_b.py)"; cat "$tmp/main.b"
  echo "== red control 1: the same count at height 22 (above the second zero, 21.022...), claimed N = 1"; grep -h "CONTROL\|^(2)\|N(\|VERDICT" "$tmp/height22.a" "$tmp/height22.b"
  echo "== red control 2: 0.755106 (below 2/log gamma_1) claimed as the upper bound"; grep -h "CONTROL\|^(3)\|VERDICT" "$tmp/printed.a" "$tmp/printed.b"
  echo "== red control 3: a bracket above gamma_1, (14.13472514173469381, 14.13472514173469383)"; grep -h "CONTROL\|^(1)\|VERDICT" "$tmp/bracket.a" "$tmp/bracket.b"
  echo "== agreement: 2/log gamma_1 (22 digits) and the integer within 1/2 of theta(17.8)/pi + 1, A against B"
  if [ -n "$(key "$tmp/main.a")" ] && [ "$(key "$tmp/main.a")" = "$(key "$tmp/main.b")" ]; then echo "agreement: True"; key "$tmp/main.a"; else echo "agreement: False"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED\|NOT SHOWN" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: the zero 1/2 + i gamma_1 with 14.13472514173469378 < gamma_1 < 14.13472514173469380 is the only zero with 0 < gamma <= 17.8; 2/log gamma_1 = 0.75510606143218443826... <= 0.755107; gamma_1^4 < 40000; 8/log(3*10^12) < 0.279; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
