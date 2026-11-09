#!/bin/sh
# Replay the 3a rows [Syl2026] k = 2 (C_3a > 1.19281555*) and [M2026] (1.19230254*, two certificates, with [K2026b]'s
# 1.19102809 as a control) exactly from their theorem statements, and the finite inputs of [O2026] (C_3a > 1.21418,
# PR #148, o2026_check.py), the [GGSWT2025] sets and the [G2026] construction; red controls included.
# Usage: ./run_all.sh [CERTS_ROOT] [--workers=N]   (default: ./certs). Python 3.9+, standard library only, and cc.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
shift 2>/dev/null
mkdir -p "$here/results"
status=0
python3 -I "$here/c3a_check.py" syl2 "$root" --controls "$@" > "$here/results/syl2.txt" 2>&1 || status=1
for row in m2026 m2026b k2026b; do
  python3 -I "$here/c3a_check.py" $row "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
done
python3 -I "$here/o2026_check.py" --controls > "$here/results/o2026.txt" 2>&1 || status=1
# the two [GGSWT2025] (AlphaEvolve) rows, 1.1479 and 1.1584, from the notebook's sets (C helper bin/ae_count)
mkdir -p "$here/bin" && cc -O2 -o "$here/bin/ae_count" "$here/ae_count.c" || status=1
python3 -I "$here/ae_check.py" "$root/ae/sums_differences_problems.ipynb" "$here/bin" --controls > "$here/results/ggswt2025.txt" 2>&1 || status=1
# the [G2026] row, 1.1740744 (PR #71), from the stated construction (exact counts of U+U and U-U; no download)
python3 -I "$here/g2026_check.py" --controls > "$here/results/g2026.txt" 2>&1 || status=1
# the [MI2026] row, 1.1835129324 (PR #95): the same programs at base 33, 420 digits, cap 1392, with the C helper
# bin/diff2d for |U-U| (about 3 hours, under 0.2 GB; the counts file is compared only after our computation)
cc -O2 -o "$here/bin/diff2d" "$here/diff2d.c" || status=1
python3 -I "$here/mi2026_check.py" "$root" "$here/bin" --controls > "$here/results/mi2026.txt" 2>&1 || status=1
for f in syl2 m2026 m2026b k2026b o2026 ggswt2025 g2026 mi2026; do grep -H "RESULT\|run time" "$here/results/$f.txt" | tr '\n' ' '; echo; done
exit $status
