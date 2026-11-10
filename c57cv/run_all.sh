#!/bin/sh
# Replay both halves of the 57c row [vS2026] U <= 0.65628589* from the files fetch.sh downloads, with red controls:
# the inradius (vs_inradius_check.py, exact, two arithmetics) and the conformal radius (vs_conformal_check.py with
# bin/vs_nodes, built here from vs_nodes.c). Output in results/. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs)
# Python 3.9+, standard library only, and a C compiler; about four minutes.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results bin
cc -O2 -ffp-contract=off -o bin/vs_nodes vs_nodes.c -lm || exit 1
status=0
python3 -I vs_inradius_check.py "$root/inradius_leaves_G.npz" --controls > results/vs_inradius.txt 2>&1 || status=1
python3 -I vs_conformal_check.py "$root/density_G_n64.json" --tex="$root/bloch_slit_discs.tex" --controls > results/vs_conformal.txt 2>&1 || status=1
grep -H "RESULT\|rho(D) >=\|U <=" results/*.txt
exit $status
