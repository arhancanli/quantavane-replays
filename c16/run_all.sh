#!/bin/sh
# Replay both sides of C16 in PR #270 (4 Nov): the upper bound from the certificate fetched by fetch.sh, and the lower
# bound of the untruncated conditional extremal at the PR's mu = 0x1.f42a1d5c1aabdp-4; red controls included.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c16_check.py" "$root/vid/certificate.json" --lower 0x1.f42a1d5c1aabdp-4 --controls > "$here/results/vid.txt" 2>&1
status=$?
grep -H "RESULT\|run time" "$here/results/vid.txt"
exit $status
