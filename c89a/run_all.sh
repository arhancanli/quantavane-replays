#!/bin/sh
# Replay the four C89 (Ramsey multiplicity of K4) upper bounds of the 89a page and PR #270 exactly: the self-tests
# (brute force over every 4-tuple on random block-structured graphons), then each row with its red controls; one
# output file per row in results/.  The first round of #270 (vid1) is replayed too.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+ (standard library only) and a C compiler;
# about 45 minutes on a loaded 18-core Mac (the 7680-class row with its controls about 27).
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/bin" "$here/results"
cc -O2 -o "$here/bin/k4_block" "$here/k4_block.c" && cc -O2 -o "$here/bin/k4_dense" "$here/k4_dense.c" || exit 1
status=0
python3 -I "$here/k4_check.py" selftest > "$here/results/selftest.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/selftest.txt" | tr '\n' ' '; echo
for row in vsj sk kos vid; do
  python3 -I "$here/k4_check.py" $row "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
python3 -I "$here/k4_check.py" vid1 "$root" > "$here/results/vid1.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/vid1.txt" | tr '\n' ' '; echo
exit $status
