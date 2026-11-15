#!/bin/sh
# 80a: the Krauth-Mezard threshold alpha_star of the zero-margin Ising perceptron from its definition (Huang,
# arXiv 2404.18902, Sec. 3; Ding and Sun, arXiv 1809.07742, Proposition 1.3), the value behind the rows [H2024] and
# [DS2025] (alpha_star ~ 0.833), and the numerical step of [AT2024] 0.847 (Altschuler and Tikhomirov, arXiv
# 2401.15092), by two separately written programs (c80a_a.py: decimal intervals, fixed point in q; c80a_b.py:
# integer fixed-point intervals with its own exp/log, fixed point in psi, the tail integrals in t = cz), then red
# controls that both must refuse. Values from their definitions, not certificates; the papers' numerical hypotheses
# (monotonicity of G_star, the AT condition, the two-variable and one-variable maximisations) are not checked.
# Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; about five minutes, under 100 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c80a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c80a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c80a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m
pair at --at-alpha 0.847
pair i1 --alpha-lo 0.8330785997 --alpha-hi 0.8330785998
pair x1 --alpha-lo 0.833078600 --alpha-hi 0.833078601
pair x2 --alpha-lo 0.833078598 --alpha-hi 0.833078599
pair x3 --variant no-prefactor
pair x4 --variant sign-reversed
pair x5 --at-alpha 0.8465
{
  echo "== 80a, alpha_star: program A (c80a_a.py)"; cat "$tmp/m.a"
  echo "== program B (c80a_b.py)"; cat "$tmp/m.b"
  echo "== [AT2024] 0.847: program A"; cat "$tmp/at.a"
  echo "== [AT2024] 0.847: program B"; cat "$tmp/at.b"
  echo "== information: the narrower interval (0.8330785997, 0.8330785998)"
  grep -h "G_star in\|G_star(\|RESULT\|VERDICT" "$tmp/i1.a" "$tmp/i1.b"
  for x in x1 x2 x3 x4 x5; do
    case $x in
      x1) echo "== red control 1: alpha_star claimed in (0.833078600, 0.833078601)";;
      x2) echo "== red control 2: alpha_star claimed in (0.833078598, 0.833078599)";;
      x3) echo "== red control 3: F_q without its factor (1-q)^(-1/2) (R_alpha without 1/(1-q))";;
      x4) echo "== red control 4: G with the sign of alpha E log Phi(cZ) reversed";;
      x5) echo "== red control 5: [AT2024] with alpha = 0.8465";;
    esac
    grep -h "CONTROL\|bracketed\|inside\|G_star in\|G_star(\|GD(\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED" "$out")
if [ "$pass" = 6 ] && [ "$ref" = 10 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: G_star(0.833078599) > 0 > G_star(0.833078600), q_star inside (q_lb, q_ub), so alpha_star in (0.833078599, 0.833078600) and alpha_star ~ 0.833; for information also in (0.8330785997, 0.8330785998); [AT2024]: GD(0.847, 1/2) < -log 2; five red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
