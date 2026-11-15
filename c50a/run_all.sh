#!/bin/sh
# 50a: [GP19] 0.498, the product-state approximation ratio for Quantum Max Cut (Gharibian and Parekh, arXiv 1909.08846v1,
# Theorem 4), recomputed from its definition by two separately written programs (c50a_a.py: the hypergeometric series
# in exact rationals; c50a_b.py: closed forms through arcsin in decimal intervals), then red controls that both must
# refuse, and the Goemans-Williamson anchor (r = 1) in program B. A value from its definition, not a certificate.
# Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; about five minutes, under 100 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c50a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c50a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c50a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m --target 0.4987
pair x1 --target 0.4988
pair x2 --target 0.4987 --g-unsquared
python3 -I "$here/c50a_b.py" --r 1 --target 0.8785 > "$tmp/gw.b" 2>&1 &
python3 -I "$here/c50a_b.py" --r 1 --target 0.8786 > "$tmp/gw2.b" 2>&1
wait
{
  echo "== 50a [GP19], r = 3: program A (c50a_a.py)"; cat "$tmp/m.a"
  echo "== program B (c50a_b.py)"; cat "$tmp/m.b"
  echo "== red control 1: the minimum claimed to be at least 0.4988"
  grep -h "cell \[\|RESULT\|VERDICT" "$tmp/x1.a" "$tmp/x1.b"
  echo "== red control 2: (Gamma(2)/Gamma(3/2)) not squared in F_3"
  grep -h "CONTROL\|least\|cell \[\|RESULT\|VERDICT" "$tmp/x2.a" "$tmp/x2.b"
  echo "== anchor (program B): r = 1, Goemans-Williamson, F_1 = (2/pi) arcsin t: at least 0.8785, and a claim of 0.8786 refused"
  grep -h "least\|RESULT\|VERDICT" "$tmp/gw.b" "$tmp/gw2.b"
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
if [ "$pass" = 3 ] && [ "$ref" = 5 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: min over [-1, 1/3) of (1 - F_3(t))/(1 - 3t) >= 0.4987, so the row's 0.498 holds; two red controls refused by both; Goemans-Williamson anchor >= 0.8785 and not >= 0.8786)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
