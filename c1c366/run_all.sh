#!/bin/sh
# Replay PR #366's C_1c certificate (K = 98,304 cells, 98,368 weights) with two separately written exact programs and
# one red control each (the weight at the deciding lag plus 1). The repository's verify.py was neither read nor run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Program A: Python 3.9+, standard library only (python3 -I).
# Program B needs numpy. About 3 minutes.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
f="$root/pr366/K98304_0.4103735.json"
out="$here/results/pr366.txt"
mkdir -p "$here/results"
{
  python3 -I "$here/c1c366_kron.py" "$f"
  python3 "$here/c1c366_limbs.py" "$f"
  echo "red control (weight 15880 plus 1):"
  python3 -I "$here/c1c366_kron.py" "$f" --bump 15880 | tail -1
  python3 "$here/c1c366_limbs.py" "$f" --bump 15880 | tail -1
} > "$out" 2>&1
a=$(grep -c "equals stated fraction: True" "$out"); b=$(grep -c "equals stated fraction: False" "$out")
if [ "$a" = 2 ] && [ "$b" = 2 ]; then echo "RESULT: ALL CHECKS PASS (both programs give the stated fraction; both red controls differ)" >> "$out"; tail -1 "$out"; exit 0; fi
echo "RESULT: FAIL" >> "$out"; tail -1 "$out"; exit 1
