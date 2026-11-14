#!/bin/sh
# Replay the 13a upper bound [K2026] (PR #314), C_13a <= area(K) = 7170601298323360123534337/29109471743520000000000000
# = 0.2463322372..., including the covering certificate (24,898-node tree, 12,663 leaves), with two separately written
# exact programs (A: c13acov_a.py, B: c13acov_b.py), then red controls that both must refuse.
# No program of the source repository was read or run. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs)
# Python 3.9+, standard library only; about three minutes (A about 100 s, B about 15 s, at most two processes; each red
# control stops at its first refusal).
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/k2026cov.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c13acov_a.py" "$root" > "$tmp/a" 2>&1 &
python3 -I "$here/c13acov_b.py" "$root" > "$tmp/b" 2>&1
wait
for x in "a h2" "a swapRL" "a overlap" "a mult" "a esc" "b h3" "b dir" "b gap" "b mult" "b pminus"; do
  set -- $x
  python3 -I "$here/c13acov_$1.py" "$root" --control $2 > "$tmp/$1.$2" 2>&1
done
{
  echo "== program A (c13acov_a.py)"; grep -v "run time" "$tmp/a"
  echo "== program B (c13acov_b.py)"; grep -v "run time" "$tmp/b"
  for x in "a h2" "a swapRL" "a overlap" "a mult" "a esc" "b h3" "b dir" "b gap" "b mult" "b pminus"; do
    set -- $x
    case $2 in
      h2) echo "== red control A1: h_2 = 130782/250000 in the escape rows (one unit lower in the last place)";;
      swapRL) echo "== red control A2: the first two right contact directions swapped";;
      overlap) echo "== red control A3: the root split's first cell widened to p <= 23, overlapping the second";;
      mult) [ $1 = a ] && echo "== red control A4: leaf 0 (B = 1 exactly), its escape multiplier one unit lower" || echo "== red control B4: leaf 0, every multiplier doubled (B = 2, sum still 0)";;
      pminus) echo "== red control B5: a direction (0, -1) from a spec given to P- instead of P+ (against the stated rule)";;
      esc) echo "== red control A5: leaf 0 given an escape row of spec 1462, which is not on its path";;
      h3) echo "== red control B1: h_3 = 446198/1000000 in the escape rows (one unit lower in the last place)";;
      dir) echo "== red control B2: the contact direction u^R_6 reflected in the x-axis";;
      gap) echo "== red control B3: the root split's first cell narrowed by one row of p (a gap)";;
    esac
    grep "CONTROL\|FAIL\|VERDICT" "$tmp/$1.$2"
  done
} > "$out" 2>&1
same=0
wa=$(grep "A: largest" "$tmp/a" | awk '{print $NF}'); wb=$(grep "B: largest" "$tmp/b" | awk '{print $NF}')
[ -n "$wa" ] && [ "$wa" = "$wb" ] && same=1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 10 ] && [ "$bad" = 0 ] && [ "$same" = 1 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs accept every leaf, with the same largest squared norm $wa and least B = 1; the area is the paper's fraction; ten red controls refused)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (pass $pass, refused $ref, same largest norm $same)" >> "$out"; tail -1 "$out"; exit 1
