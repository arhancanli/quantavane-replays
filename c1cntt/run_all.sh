#!/bin/sh
# Replay of PR #354's K = 2^22 and K = 2^24 certificates for C_1c with our own NTT code: exact autocorrelation
# minimum, stated fraction and lag, printed bounds, and controls. One output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). A C compiler and Python 3.9+, standard library only;
# about 3 minutes and under 1 GB of memory.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results" "$here/bin"
cc -O2 -o "$here/bin/c1c_ntt" "$here/c1c_ntt.c" || exit 1
python3 -I "$here/c1cntt_check.py" "$root" "$here/bin/c1c_ntt" --controls > "$here/results/sag2026_k22_k24.txt" 2>&1
status=$?
grep -H "RESULT\|FINDING\|controls all" "$here/results/sag2026_k22_k24.txt"
exit $status
