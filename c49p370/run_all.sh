#!/bin/sh
# Replay the 49a block of PR #370 (252,303 sets on 28 points, mu_3 >= 1.5678643*) with two separately written exact
# programs: A = our 4 November checker ../c49/sunflower_check.py, unchanged, through c49p370_a.py (box method on two
# point partitions, its red controls); B = c49p370_b.py (decision diagram and memoised pair and triple joins; its
# own self-test against the definitions; its own red controls). No program of any source was read or run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about two minutes, one process.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/pin370.txt"
F="$root/pr370/best_1.5678643845728069793620.json"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c49p370_a.py" "$root" --controls > "$tmp/a" 2>&1
python3 -I "$here/c49p370_b.py" --selftest > "$tmp/bself" 2>&1
python3 -I "$here/c49p370_b.py" "$F" > "$tmp/b" 2>&1
for c in sunflower2 subset drop printed; do
  python3 -I "$here/c49p370_b.py" "$F" --control $c > "$tmp/b.$c" 2>&1
done
{
  echo "== program A (../c49/sunflower_check.py, unchanged, via c49p370_a.py), with its red controls"
  grep -v "run time" "$tmp/a"
  echo "== program B (c49p370_b.py): self-test, then the certificate"
  cat "$tmp/bself"; grep -v "run time" "$tmp/b" | sed -E 's/, [0-9.]+ s\)/)/'
  for c in sunflower2 subset drop printed; do
    case $c in
      sunflower2) echo "== B red control 1: one set planted that makes a 3-sunflower with two members and is comparable with none";;
      subset) echo "== B red control 2: a member with one point removed (a proper subset of a member)";;
      drop) echo "== B red control 3: the added set 64768 removed (N00 = 33979)";;
      printed) echo "== B red control 4: the printed bound one unit higher, 1.5678644";;
    esac
    grep "CONTROL\|planted\|FAIL\|VERDICT" "$tmp/b.$c"
  done
} > "$out" 2>&1
rm -rf "$tmp"
a_ok=$(grep -c "RESULT: ALL CHECKS PASS" "$out"); pass=$(grep -c "B: VERDICT PASS" "$out"); ref=$(grep -c "B: VERDICT REFUSED" "$out")
bad=$(grep -c -i "traceback" "$out")
same=$(grep -c "1.5678643845728069793620" "$out")
if [ "$a_ok" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 4 ] && [ "$bad" = 0 ] && [ "$same" -ge 3 ]; then
  echo "RESULT: ALL CHECKS PASS (A and B: Sperner and 3-sunflower-free, port counts 33980, 64728, 64728, 88867, rho^(1/26) = 1.5678643845728069793620...; every red control refused)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (A $a_ok, B pass $pass, B refused $ref)" >> "$out"; tail -1 "$out"; exit 1
