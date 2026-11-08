#!/bin/sh
# Replays of two 30a rows: [N2026] gr(Av(1324)) >= 10.412263 (PR #313), the evaluation of Norton's Table 1 from the
# paper's LaTeX source (strip functions differentiated exactly and compared with (H2), (H3); profile; interval
# evaluation of the bracket, its slopes and the tangent bound; Phi > 0), and the candidate 521/50 of PR #317 (the exact
# numerical point of its Section 6, with its envelope profile), with controls. One output file per row in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 10 seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/norton_check.py" "$root" --controls > "$here/results/n2026.txt" 2>&1 || status=1
python3 -I "$here/norton_check.py" "$root" --controls --row=pin2026 > "$here/results/pin2026.txt" 2>&1 || status=1
grep -H "RESULT\|controls all" "$here/results/n2026.txt" "$here/results/pin2026.txt"
exit $status
