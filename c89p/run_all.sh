#!/bin/sh
# Replays of three 89a rows (all added in PR #255): [PPSS2024], the 768-vertex graph of Parczyk, Pokutta, Spiegel and
# Szabo, its K_4 homomorphism density (with the looped complement) by two routes in k4hom.c, exactly
# 4551721 * 2^-24 * 3^-2 (ppss_check.py); [Tho1997] and [EZL2015], Thomason's product and the nested blow-up of
# Even-Zohar and Linial, exact 4-point profiles and k4hom on the explicit graphs (ezl_check.py); red controls.
# Needs ./fetch.sh first. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+ and a C compiler; about
# two minutes (three full counts of the 768-vertex graph and a random graph take most of it). Peak memory under 50 MB.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results bin
cc -O2 -Wall -o bin/k4hom k4hom.c || exit 1
status=0
python3 -I ppss_check.py "$root" bin --controls > results/ppss2024.txt 2>&1 || status=1
grep -H "RESULT\|controls all" results/ppss2024.txt
python3 -I ezl_check.py "$root" bin --controls > results/tho1997_ezl2015.txt 2>&1 || status=1
grep -H "RESULT\|controls all" results/tho1997_ezl2015.txt
exit $status
