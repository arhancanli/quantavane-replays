#!/bin/sh
# Replay the 3a rows [Syl2026] k = 2 (C_3a > 1.19281555*) and [M2026] (1.19230254*, two certificates, with [K2026b]'s
# 1.19102809 as a control) exactly from their theorem statements, and the finite inputs of [O2026] (C_3a > 1.21418,
# PR #148, o2026_check.py); red controls included.
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
for f in syl2 m2026 m2026b k2026b o2026; do grep -H "RESULT\|run time" "$here/results/$f.txt" | tr '\n' ' '; echo; done
exit $status
