#!/bin/sh
# 57a [AG1937] and 57b [Rad1943]: the closed forms Gamma(1/3) Gamma(11/12)/(sqrt(1 + sqrt 3) Gamma(1/4)) ~ 0.4719
# and Gamma(1/3) Gamma(5/6)/Gamma(1/6) ~ 0.5433 against the page's decimals, by two separately written programs
# (c57ab_a.py: Stirling's series in integer fixed point; c57ab_b.py: the incomplete-gamma series in decimal
# intervals), then red controls that both must refuse. Closed forms named in the rows; the extremal-function
# analysis behind them is not checked.
# Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c57ab.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c57ab_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c57ab_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m
pair x1 --variant ag-no-factor
pair x2 --variant rad-swapped
pair x3 --print-ag 0.4718
vals() { grep -h "^row \|^Gamma(1/3)\|^Gamma(1/4)\|^Gamma(11/12)\|^Gamma(5/6)\|^Gamma(1/6)" "$tmp/m.$1" | grep -o "[0-9]\.[0-9]\{20\}" | sort; }
{
  echo "== program A (c57ab_a.py)"; cat "$tmp/m.a"
  echo "== program B (c57ab_b.py)"; cat "$tmp/m.b"
  echo "== red control 1: the Ahlfors-Grunsky form without 1/sqrt(1 + sqrt 3)"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x1.a" "$tmp/x1.b"
  echo "== red control 2: Rademacher's form with Gamma(1/6) and Gamma(5/6) swapped"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x2.a" "$tmp/x2.b"
  echo "== red control 3: the Ahlfors-Grunsky value printed as 0.4718"; grep -h "CONTROL\|row \|no sign change\|VERDICT" "$tmp/x3.a" "$tmp/x3.b"
  echo "== agreement: the first 20 digits of every value, A against B"
  va=$(vals a); vb=$(vals b)
  echo "$va" | tr "\n" " "; echo
  if [ -n "$va" ] && [ "$va" = "$vb" ]; then echo "agreement: True"; else echo "agreement: False"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED" "$out")
agree=$(grep -c "^agreement: True" "$out")
if [ "$agree" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: Gamma(1/3), Gamma(1/4), Gamma(11/12), Gamma(5/6), Gamma(1/6) and both rows agree to 20 digits; [AG1937] = 0.47186165345268178487... ~ 0.4719, [Rad1943] = 0.54325896534297670695... ~ 0.5433; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
