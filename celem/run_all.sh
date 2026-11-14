#!/bin/sh
# Replay five elementary rows on four pages that had no replay folder (20b: the cube and the simplices; 25a: the cube;
# 39a: the cube; 11a: n = 1 and f(x) = x) with two separately written exact programs, then five red controls that both
# must refuse. Usage: ./run_all.sh   (Python 3.9+, standard library only; about ten seconds).
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/celem.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/celem_a.py" > "$tmp/main.a" 2>&1 &
python3 -I "$here/celem_b.py" > "$tmp/main.b" 2>&1
wait
for c in cube20b simplex polar vertex poincare; do
  python3 -I "$here/celem_a.py" --control $c > "$tmp/$c.a" 2>&1 &
  python3 -I "$here/celem_b.py" --control $c > "$tmp/$c.b" 2>&1
  wait
done
{
  echo "== the five rows"; cat "$tmp/main.a" "$tmp/main.b"
  for c in cube20b simplex polar vertex poincare; do
    case $c in
      cube20b) echo "== red control 1: 20b, the cube of side 2 with its volume left at 1";;
      simplex) echo "== red control 2: 20b, one second moment of the 5-simplex raised by 1/1000";;
      polar) echo "== red control 3: 25a, the cube taken as its own polar";;
      vertex) echo "== red control 4: 39a, the 16 vertices of a 4-cube projected to three coordinates (repeated points)";;
      poincare) echo "== red control 5: 11a, D_j f taken without its factor 1/2";;
    esac
    grep -h "CONTROL\|FAIL\|VERDICT" "$tmp/$c.a" "$tmp/$c.b"
  done
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 10 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs accept the five rows; five red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (pass $pass, refused $ref)" >> "$out"; tail -1 "$out"; exit 1
