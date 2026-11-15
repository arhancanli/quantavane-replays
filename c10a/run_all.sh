#!/bin/sh
# 10a: the Davie-Reeds value K_DR ~ 1.67696 ([Dav1984], [Ree1991]) from its definition (Jones and Malavolta,
# arXiv 2603.30039: C* the root of 4 phi(C)^2 - 4 Phi(-C) + 1 in (0, 1), K_DR = 1/R(C*)), and the closed forms of
# [G1953] sinh(pi/2) and pi/2, [K1979] pi/(2 log(1 + sqrt 2)) and [SLXCKKM26] 6 pi/11 against the page's decimals,
# by two separately written programs (c10a_a.py: integer fixed point, K_DR = (pi/2) exp(C*^2) at the root;
# c10a_b.py: decimal intervals, K_DR = 1/R(C) by the full formula), then red controls that both must refuse.
# Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c10a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c10a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c10a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m
pair x1 --print-kdr 1.67695
pair x2 --variant phi-plus
pair x3 --print-sinh 2.30131
vals() { grep -h "^row \|^K_DR" "$tmp/m.$1" | grep -o "[0-9]\.[0-9]\{20\}" | sort; }
{
  echo "== program A (c10a_a.py)"; cat "$tmp/m.a"
  echo "== program B (c10a_b.py)"; cat "$tmp/m.b"
  echo "== red control 1: K_DR printed as 1.67695"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x1.a" "$tmp/x1.b"
  echo "== red control 2: Phi(C) in place of Phi(-C) in the equation for C*"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x2.a" "$tmp/x2.b"
  echo "== red control 3: sinh(pi/2) printed as 2.30131"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x3.a" "$tmp/x3.b"
  echo "== agreement: the first 20 digits of every value, A against B"
  va=$(vals a); vb=$(vals b)
  echo "$va" | tr "\n" " "; echo
  if [ -n "$va" ] && [ "$va" = "$vb" ]; then echo "agreement: True"; else echo "agreement: False"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: C* = 0.2557302131662100143618..., K_DR = 1.6769566742155762370778... ~ 1.67696; sinh(pi/2), pi/2, Krivine and 6 pi/11 round to the page's decimals; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
