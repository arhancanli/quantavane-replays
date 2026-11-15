#!/bin/sh
# 10a [P2026] K_G <= 1.7813* (Peng, arXiv 2609.20074v1): a replay of the computational lemma (Lemma 5.1: b_1 > 0.881850816,
# sum_{3<=n<=301} |b_n| < 3e-8, B_3 < 107) and of the closing inequality gamma + sum_{n>=3} |b_n| < b_1, from the paper's
# text and its exact parameters only, by two separately written programs on Arb (python-flint):
#   p26_a.py: the Hermite head by composite Gauss-Legendre with Bernstein-ellipse error bounds (all 22,952 pairs at
#             once); the integral I of the integrand displayed in Section 5.2 by an adaptive partition with ball upper bounds (the paper's F_ij
#             recurrence, Bell polynomials by Faa di Bruno, Table 1); its own exterior bound;
#   p26_b.py: the head by Arb's Hermite function, a different panel layout and the derivative form of the
#             Gauss-Legendre remainder with Cauchy estimates; b_n by Horner in rho; on A's partition (whose exact
#             tiling it checks) its own integrand (closed-form derivatives of the Gaussian, Bell by series
#             composition) and its own exterior bound.
# Both check the identity D^3 H(t) = T_(-t) Phi_3(t) at real t against their own heads. Then red controls.
# Nothing is downloaded; no program of the source was read or run. Analytic lemmas of the paper not re-proved.
# Usage: PYTHON=/path/to/python ./run_all.sh   needs python-flint 0.6.0 (pip install python-flint==0.6.0, PyPI);
# about 90 minutes on a loaded 18-core Mac, two processes, under 1 GB.
here=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
out="$here/results/c10ap26.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
v=$($PY -c "import flint; print(flint.__version__)" 2>/dev/null)
if [ "$v" != "0.6.0" ]; then echo "python-flint 0.6.0 is required (found: '$v'); see USAGE in run_all.sh"; exit 1; fi
A="$PY -I $here/p26_a.py"; B="$PY -I $here/p26_b.py"
cd "$tmp" || exit 1
$A head head_a.json > a_head.txt 2>&1 &
$B head head_b.json > b_head.txt 2>&1; wait
$A selftest head_a.json > a_self.txt 2>&1 &
$B selftest head_b.json > b_self.txt 2>&1; wait
$A b3 b3_a.json boxes.txt.gz 150 1000000 1e9 > a_b3.txt 2>&1
$B b3 boxes.txt.gz b3_b0.json 0/2 > b_b3.txt 2>&1 &
$B b3 boxes.txt.gz b3_b1.json 1/2 > b_b3_1.txt 2>&1; wait
$B b3sum b3_b.json b3_b0.json b3_b1.json >> b_b3.txt 2>&1
$A final head_a.json b3_a.json > a_final.txt 2>&1 &
$B final head_b.json b3_b.json > b_final.txt 2>&1; wait
# red controls
$A final head_a.json b3_a.json 17812/10000 > x1.a 2>&1 &
$B final head_b.json b3_b.json 17812/10000 > x1.b 2>&1; wait
$A head head_ax.json --c1-shift -1/10000 > x2h.a 2>&1 &
$B head head_bx.json --c1-shift -1/10000 > x2h.b 2>&1; wait
$A final head_ax.json b3_a.json > x2.a 2>&1 &
$B final head_bx.json b3_b.json > x2.b 2>&1; wait
$A selftest head_a.json --kappa-sign > x3.a 2>&1 &
$B selftest head_b.json --kappa-sign > x3.b 2>&1; wait
{
  for f in a_head b_head a_self b_self a_b3 b_b3 b_b3_1 a_final b_final; do echo "== $f"; cat $f.txt; done
  echo "== red control 1: the claim K_G <= 1.7812"; cat x1.a x1.b
  echo "== red control 2: c_1 lowered by 1/10000 (heads recomputed)"; grep -h "CONTROL\|b_1 in\|VERDICT\|margin" x2h.a x2h.b x2.a x2.b
  echo "== red control 3: kappa_(1,1) with its sign flipped (the integrand self-tests)"; grep -h "CONTROL\|self-test\|VERDICT" x3.a x3.b
} > "$out" 2>&1
cd "$here"; rm -rf "$tmp"
pass=$(grep -c "^VERDICT PASS" "$out"); self=$(grep -c "^VERDICT SELFTEST PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out")
bad=$(grep -c -i "traceback\|error" "$out"); cov=$(grep -c "^coverage: True" "$out")
if [ "$pass" = 2 ] && [ "$self" = 2 ] && [ "$ref" = 6 ] && [ "$cov" = 1 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: b_1, the head sum and B_3 as in Lemma 5.1, and gamma + head + tail < b_1, so K_G <= 1.7813 given the paper's analytic lemmas; the partition tiles the domain; both integrands pass the D^3 H identity; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass selftest=$self refused=$ref coverage=$cov errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
