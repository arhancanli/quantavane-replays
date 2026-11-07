#!/bin/sh
# Replay [J2026]'s C_91 >= 0.216276 (PR #273): the self-tests (canonical form against the edit distance, block
# products against enumeration of all letters, the automaton against brute-force expected counts), then the
# certificate with its red controls; one output file each in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
python3 -I "$here/ed_check.py" selftest > "$here/results/selftest.txt" 2>&1 || status=1
python3 -I "$here/ed_check.py" cert "$root" --controls > "$here/results/j2026.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/selftest.txt" "$here/results/j2026.txt"
exit $status
