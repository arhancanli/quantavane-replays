#!/bin/sh
# 52a: the first-moment upper bound [FP1983] 5.191 for the random 3-SAT threshold, from the page's model (m clauses
# uniformly with replacement from the 8 C(n,3) clauses on three distinct variables), by two separately written
# programs (c52a_a.py: exhaustive count over all formulas for small n and m, then exact integer powers;
# c52a_b.py: the satisfied fraction of the clauses for every assignment, n <= 9, then decimal intervals), then red
# controls that both must refuse. A value from its definition (the row's own argument), not a certificate.
# Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c52a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c52a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c52a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m
pair x1 --claim 5.190
pair x2 --variant two-literal
{
  echo "== 52a [FP1983]: program A (c52a_a.py)"; cat "$tmp/m.a"
  echo "== program B (c52a_b.py)"; cat "$tmp/m.b"
  echo "== red control 1: the bound claimed as 5.190"
  grep -h "< 1 (exact\|row \|VERDICT" "$tmp/x1.a" "$tmp/x1.b"
  echo "== red control 2: clauses on two variables (each satisfied with probability 3/4)"
  grep -h "CONTROL\|False\|VERDICT" "$tmp/x2.a" "$tmp/x2.b"
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 4 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: E[#satisfying assignments] = 2^n (7/8)^m in the page's model, r1 = log 2/log(8/7) = 5.1908930696... < 5.191; two red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
