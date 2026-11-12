#!/bin/sh
# Replay of two 3a lower-bound rows through the Gyarmati-Hennecart-Ruzsa lemma: [Ru96] U = {0, 1, 3} (exact counts) and
# [Z2025] Zheng's limit bound theta >= 1.173077* (his table reproduced; the bound certified at B = 5, r0 = 0.805053,
# a0 = 0.466823 in 60-digit interval arithmetic); red controls. Needs ./fetch.sh first (Zheng's LaTeX source).
# Also [GHR2007] U = {0,1,3,6,13,17,21} (1.1078...) and [K2026c] (the transfer identity, the decimals and the README
# cell of PR #364), which need no download (ghr_k2026c_check.py).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 15 seconds.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results
python3 -I z2025_check.py "$root" --controls > results/z2025_ru96.txt 2>&1; status=$?
python3 -I ghr_k2026c_check.py --controls > results/ghr2007_k2026c.txt 2>&1 || status=1
grep -H "RESULT\|controls all" results/z2025_ru96.txt results/ghr2007_k2026c.txt
exit $status
