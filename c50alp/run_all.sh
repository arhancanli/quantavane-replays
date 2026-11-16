#!/bin/sh
# 50a, the lower bound [LP24] 0.595 (Lee and Parekh, arXiv 2401.03616, Theorem 2): min over x in [-1, 1/3) of
# rho_p(x) = [p (1 - F_3(x)) + (1 - p)(1 + (6/5) max(0, -(1 + 3x)))]/(1 - 3x) at the paper's p = 0.674, certified >= 0.595
# by two separately written programs (c50alp_a.py: hypergeometric series in exact rationals; c50alp_b.py: closed forms
# through arcsin in decimal intervals; both from our [GP19] programs in ../c50a). Red controls: the claim 0.5958 (both);
# p = 0.9 (A) and p = 0.5 (B). Nothing is downloaded; no program of any source read or run.
# Usage: ./run_all.sh   Python 3.9+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c50alp.txt"
mkdir -p "$here/results"
{
  echo "== program A (c50alp_a.py)"; python3 -I "$here/c50alp_a.py"
  echo "== program B (c50alp_b.py)"; python3 -I "$here/c50alp_b.py"
  echo "== red control 1: the claim 0.5958 (the minimum is 0.59572...), A then B"
  python3 -I "$here/c50alp_a.py" --target 0.5958 | grep "does not\|stopped\|RESULT\|VERDICT"
  python3 -I "$here/c50alp_b.py" --target 0.5958 | grep "does not\|cell of width\|RESULT\|VERDICT"
  echo "== red control 2: p = 0.9 (A) and p = 0.5 (B) in place of 0.674"
  python3 -I "$here/c50alp_a.py" --p 9/10 | grep "does not\|stopped\|RESULT\|VERDICT"
  python3 -I "$here/c50alp_b.py" --p 0.5 | grep "does not\|cell of width\|RESULT\|VERDICT"
} > "$out" 2>&1
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 4 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: min over x of rho_0.674 >= 0.595, so the max over p is at least 0.595; red controls refused)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
