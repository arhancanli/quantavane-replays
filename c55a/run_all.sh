#!/bin/sh
# 55a: the computational inputs of three upper bounds on the acyclic chromatic index coefficient, from the papers'
# text: [KLS2026] 3.142 (rho >= 0.6677 for even-outdegree >= 4 trees, and gamma = 1.142), [FLM2020] 3.569 and
# [GKPT2017] 3.74 (min phi(x)/x < 1), by two separately written programs (c55a_a.py: multiset table, exact interval
# arithmetic; c55a_b.py: the functional equation, decimal intervals), then red controls that both must refuse.
# Values from their definitions, not certificates. Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; about five minutes, under 200 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c55a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c55a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c55a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m
pair x1 --x0 0.66771
pair x2 --gamma 1.141
pair x3 --variant with-degree-2
{
  echo "== 55a: program A (c55a_a.py)"; cat "$tmp/m.a"
  echo "== program B (c55a_b.py)"; cat "$tmp/m.b"
  for x in x1 x2 x3; do
    case $x in
      x1) echo "== red control 1: the supersolution at x = 0.66771, above rho";;
      x2) echo "== red control 2: gamma = 1.141 in place of 1.142";;
      x3) echo "== red control 3: the tree class with outdegree 2 allowed (the MSET_2 term not removed)";;
    esac
    grep -h "CONTROL\|supersolution at\|gamma = \|no supersolution\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
  echo "== agreement"
  ta=$(grep "^t_n, n = 1..41" "$tmp/m.a"); tb=$(grep "^t_n, n = 1..41" "$tmp/m.b")
  sa=$(grep "^t_99 " "$tmp/m.a"); sb=$(grep "^t_99 " "$tmp/m.b")
  echo "A: $sa"; echo "B: $sb"
  if [ -n "$ta" ] && [ "$ta" = "$tb" ] && [ -n "$sa" ] && [ "$sa" = "$sb" ] \
     && grep -q "supersolution at x = 6677/10000: True" "$tmp/m.a" && grep -q "supersolution at x = 0.6677: True" "$tmp/m.b"; then
    echo "pair exit 0"; else echo "pair exit 1"; fi
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
agree=$(grep -c "^pair exit 0" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ] && [ "$agree" = 1 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: the same t_n; rho >= 0.6677 by a supersolution; (1 - e^(-1/1.142))^(3/4) < 0.6677, so 3.142; FLM2020 at gamma = 1.569 and GKPT2017 at gamma = 1.74 give min phi(x)/x < 1; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
