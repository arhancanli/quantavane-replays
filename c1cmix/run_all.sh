#!/bin/sh
# 1c, the remark of PR #354 on test-measure upper bounds: replay the mixture certificate (two step functions with
# K = 4096 steps and rational weights lam_1 + lam_2 = 1) with two separately written exact programs (A: direct sums in
# Python integers and Fractions; B: 16-bit limbs with exact 64-bit correlations, integers only), then the smaller
# K = 1024 certificate, then three red controls. The repository's two verifiers were neither read nor run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+; program B needs numpy. About ten seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/mix.txt"
mkdir -p "$here/results"
V="5279542067372831214509056145021140361894623145945986164923459 / 12855504354070095327289033191535703182400661039125234952783872"
A="python3 -I $here/c1cmix_direct.py"; B="python3 $here/c1cmix_limbs.py"
{
  echo "== K = 4096 (the statement's certificate)"
  $A "$root/mix/mix_K4096_P2.json"; $B "$root/mix/mix_K4096_P2.json"
  echo "stated V: $V"
  echo "== K = 1024 (the smaller certificate; its README gives 0.4105305028)"
  $A "$root/mix/mix_K1024_P2.json"; $B "$root/mix/mix_K1024_P2.json"
  echo "== red control 1 (K = 4096): the weight at index 780 of the first component plus 1"
  $A "$root/mix/mix_K4096_P2.json" --bump 0 780 | grep "V = [0-9]* /"; $B "$root/mix/mix_K4096_P2.json" --bump 0 780 | grep "V = [0-9]* /"
  echo "== red control 2 (K = 4096): the two values of lam exchanged"
  $A "$root/mix/mix_K4096_P2.json" --swap-lam | grep "V = "; $B "$root/mix/mix_K4096_P2.json" --swap-lam | grep "V = "
  echo "== red control 3 (K = 4096): the first weight of the second component plus 1"
  $A "$root/mix/mix_K4096_P2.json" --bump 1 0 | grep "V = [0-9]* /"; $B "$root/mix/mix_K4096_P2.json" --bump 1 0 | grep "V = [0-9]* /"
} > "$out" 2>&1
main=$(grep -c "V = $V\$" "$out"); pass=$(grep -c "VERDICT PASS" "$out"); k780=$(grep -c "minimum at k = \[780\]" "$out")
if [ "$main" = 2 ] && [ "$pass" = 4 ] && [ "$(grep -c -i traceback "$out")" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs give the stated V; the fraction appears $main times in all, so every red control differs; k = 780 lines: $k780)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (stated V found $main times, pass $pass)" >> "$out"; tail -1 "$out"; exit 1
