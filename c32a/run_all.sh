#!/bin/sh
# 32a [Li25]: the closed form sum_{k>=1} 2^(-k-1) k log2 k - log2 log2 e ~ 0.76 against the page's decimal, by
# two separately written programs (c32a_a.py: integer fixed point with its own log; c32a_b.py: decimal ln with
# directed rounding), each with an explicit tail bound, then red controls that both must refuse.
# Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c32a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c32a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c32a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m
pair x1 --variant no-half
pair x2 --print 0.75
vals() { grep -h "row \[Li25\]\|^value in" "$tmp/m.$1" | grep -o "[0-9]\.[0-9]\{20\}" | head -1; }
{
  echo "== program A (c32a_a.py)"; cat "$tmp/m.a"
  echo "== program B (c32a_b.py)"; cat "$tmp/m.b"
  echo "== red control 1: 2^(-k) in place of 2^(-k-1)"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x1.a" "$tmp/x1.b"
  echo "== red control 2: the value printed as 0.75"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x2.a" "$tmp/x2.b"
  echo "== agreement: the first 20 digits of every value, A against B"
  va=$(vals a); vb=$(vals b)
  echo "$va" | tr "\n" " "; echo
  if [ -n "$va" ] && [ "$va" = "$vb" ]; then echo "agreement: True"; else echo "agreement: False"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 4 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: 0.75976490283449090729718935... ~ 0.76; two red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
