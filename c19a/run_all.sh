#!/bin/sh
# 19a [E1956]: Esseen's lower bound c_E = (sqrt 10 + 3)/(6 sqrt(2 pi)) ~ 0.4097321837 from its definition (the
# maximum over normalized two-point laws of Esseen's lattice limit divided by beta_3), by two separately written
# programs (c19a_a.py: the maximiser in closed form by exact algebra in Q(sqrt 10); c19a_b.py: branch-and-bound over p
# with mean-value bounds, no calculus), then three red controls that both must refuse.
# Nothing is downloaded; no program of any source was read or run. Usage: ./run_all.sh   Python 3.9+; a second.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c19a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
for n in main print lattice beta; do
  arg=""; [ "$n" != main ] && arg="--control $n"
  python3 -I "$here/c19a_a.py" $arg > "$tmp/$n.a" 2>&1
  python3 -I "$here/c19a_b.py" $arg > "$tmp/$n.b" 2>&1
done
val() { grep -h "^row 19a" "$1" | grep -o "c_E = [0-9]\.[0-9]\{30\}"; }
{
  echo "== program A (c19a_a.py)"; cat "$tmp/main.a"
  echo "== program B (c19a_b.py)"; cat "$tmp/main.b"
  echo "== red control 1: the value printed as 0.4097321838"; grep -h "CONTROL\|^row\|VERDICT" "$tmp/print.a" "$tmp/print.b"
  echo "== red control 2: the lattice term h/2 doubled"; grep -h "CONTROL\|^row\|VERDICT" "$tmp/lattice.a" "$tmp/lattice.b"
  echo "== red control 3: the factor p^2 + q^2 of beta_3 dropped"; grep -h "CONTROL\|^row\|VERDICT" "$tmp/beta.a" "$tmp/beta.b"
  echo "== agreement: the first 30 digits, A against B (main run and the two changed definitions)"
  for n in main lattice beta; do
    va=$(val "$tmp/$n.a"); vb=$(val "$tmp/$n.b")
    if [ -n "$va" ] && [ "$va" = "$vb" ]; then echo "$n: equal ($va)"; else echo "$n: DIFFERENT"; fi
  done
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED\|DIFFERENT" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: c_E = 0.409732183702396342990568823497... ~ 0.4097321837; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
