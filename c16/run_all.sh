#!/bin/sh
# Replay both sides of C16 in PR #270 (4 Nov): the upper bound from the certificate fetched by fetch.sh, and the lower
# bound of the untruncated conditional extremal at the PR's mu = 0x1.f42a1d5c1aabdp-4; red controls included. Then the
# lower bound of PR #319 (Julsgaard), the same family at mu = 128042/2^20.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; under a minute.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c16_check.py" "$root/vid/certificate.json" --lower 0x1.f42a1d5c1aabdp-4 --controls > "$here/results/vid.txt" 2>&1
status=$?
grep -H "RESULT\|run time" "$here/results/vid.txt"
# PR #319 (Julsgaard): the same trial function at his mu = 128042/2^20 (no certificate file needed)
python3 -I "$here/c16_check.py" --lower-only 128042/1048576 > "$here/results/julsgaard_319.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/julsgaard_319.txt"
exit $status
