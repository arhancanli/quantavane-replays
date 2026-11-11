#!/bin/sh
# Replay of PR #283's statement (37a): no Boolean function of degree 4 has sensitivity 10, from the definition and the
# certificate files alone (fetched by fetch.sh; read as data), with red controls; one output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; a few seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c37_check.py" "$root/kosh" --controls > "$here/results/kos2026.txt" 2>&1; status=$?
grep -H "RESULT\|run time" "$here/results/kos2026.txt"
# the two explicit lower-bound rows of the page ([HKP2011] Kushilevitz's h, [BdW2002]/[T2013] not-all-equal g), from
# the polynomials printed on the page (./fetch_page.sh): exhaustive on the cube, with controls; results/hkp2011_bdw2002.txt
python3 -I "$here/kushilevitz_check.py" "$root" --controls > "$here/results/hkp2011_bdw2002.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/hkp2011_bdw2002.txt"
exit $status
