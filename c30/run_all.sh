#!/bin/sh
# Replay of the 30a row [N2026], gr(Av(1324)) >= 10.412263 (PR #313): the evaluation of Norton's Table 1 from the
# paper's LaTeX source (strip functions differentiated exactly and compared with (H2), (H3); profile; interval
# evaluation of the bracket, its slopes and the tangent bound; Phi > 0), with controls. One output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 10 seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/norton_check.py" "$root" --controls > "$here/results/n2026.txt" 2>&1; status=$?
grep -H "RESULT\|controls all" "$here/results/n2026.txt"
exit $status
