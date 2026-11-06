#!/bin/sh
# Replay both sides of C16 in PR #270 (4 Nov): the upper bound from the certificate fetched by fetch.sh, and the lower
# bound of the untruncated conditional extremal at the PR's mu = 0x1.f42a1d5c1aabdp-4; red controls included. Then the
# lower bound of PR #319 (Julsgaard), the same family at mu = 128042/2^20, and the night update of #270 (c16_pinned_v2,
# 18 digits).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 7 minutes.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
python3 -I "$here/c16_check.py" "$root/vid/certificate.json" --lower 0x1.f42a1d5c1aabdp-4 --controls > "$here/results/vid.txt" 2>&1
status=$?
grep -H "RESULT\|run time" "$here/results/vid.txt"
# PR #270 update of 5 Nov, night (c16_pinned_v2): both sides again, G and H enclosed to about 1e-21 near the peak
python3 -I "$here/c16_check.py" "$root/vid2/certificate.json" --lower 0x1.f42a1bd5af466p-4 --controls \
    --upper-width 1e-21 --max-rows 262144 --lower-factor 128 > "$here/results/vid2.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/vid2.txt"
# PR #319 (Julsgaard): the same trial function at his mu = 128042/2^20 (no certificate file needed)
python3 -I "$here/c16_check.py" --lower-only 128042/1048576 > "$here/results/julsgaard_319.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/julsgaard_319.txt"
# PR #326 (16a row, maximusrufus): the trial function truncated at R = 4000, mu = 1221/10000; red controls included
python3 -I "$here/c16_trunc_check.py" "$root/maxr/lower_certificate.json" --controls > "$here/results/maxr_trunc.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/maxr_trunc.txt"
exit $status
