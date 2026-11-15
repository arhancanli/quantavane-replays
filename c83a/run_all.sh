#!/bin/sh
# 83a: rigorous bounds for the Wirsing constant C_83 = |lambda_2| of the Gauss-Kuzmin-Wirsing operator, the value
# behind the rows [Schweiger1980] <= 0.38197 and [MR1987] 0.30366299 <= C_83 <= 0.30366327. Both programs prove
# lo psi <= V psi <= hi psi on [0, 1] for Wirsing's positive operator V and the polynomial psi of psi.json (found by
# c83a_search.py, not trusted), hence lo <= |lambda_2| <= hi (Collatz-Wielandt; rho(V) = |lambda_2| is Wirsing's
# theorem, not re-proved). c83a_a.py: partial fractions and Hurwitz zeta by Euler-Maclaurin, exact rational
# intervals, Taylor cells; c83a_b.py: the summand as a power series in 1/(x+k+1/2), Riemann zeta at integers by
# alternating-series acceleration, one global Taylor polynomial, decimal intervals. Then red controls that both must
# refuse. A value from its definition. Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; about five minutes, under 100 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c83a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c83a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c83a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
python3 -I "$here/c83a_search.py" --out "$tmp/psi_again.json" > "$tmp/search.txt" 2>&1
pair m
pair x1 --claim-lower 0.303663002898733
pair x2 --claim-upper 0.303663002898732
pair x3 --perturb 9:1/1073741824
pair x4 --variant wrong-weight
{
  echo "== 83a: the search step (c83a_search.py, floating search, not part of the proof) run again"; cat "$tmp/search.txt"
  if cmp -s "$tmp/psi_again.json" "$here/psi.json"; then echo "psi.json reproduced byte for byte: True"; else echo "psi.json reproduced byte for byte: False"; fi
  echo "== program A (c83a_a.py)"; cat "$tmp/m.a"
  echo "== program B (c83a_b.py)"; cat "$tmp/m.b"
  for x in x1 x2 x3 x4; do
    case $x in
      x1) echo "== red control 1: the lower bound claimed as 0.303663002898733";;
      x2) echo "== red control 2: the upper bound claimed as 0.303663002898732";;
      x3) echo "== red control 3: coefficient 9 of psi changed by 2^-30";;
      x4) echo "== red control 4: the point term of V with (x+k)^2 in place of (x+k)^3";;
    esac
    grep -h "CONTROL\|>= 0 on\|FAILS\|not proved\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED" "$out")
same=$(grep -c "reproduced byte for byte: True" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 8 ] && [ "$same" = 1 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: psi > 0 and 0.303663002898732 psi <= V psi <= 0.303663002898733 psi on [0, 1], so 0.303663002898732 <= C_83 <= 0.303663002898733: [Schweiger1980] 0.38197 and [MR1987] 0.30366299, 0.30366327 hold; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad same=$same)" >> "$out"; tail -1 "$out"; exit 1
fi
