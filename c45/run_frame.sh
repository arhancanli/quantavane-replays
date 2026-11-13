#!/bin/sh
# Replay of the C45 row [S2026] of PR #290 (29 base primes in the frame T0 = 12252240, 25 cluster primes):
# C45 < 0.489481743486002. Our own code (c45_frame.c, c45_frame_check.py, with the group laws of c45_check.py); the
# claimant's histogram, fetched by fetch.sh, is read as data only and compared with ours bin by bin.
# Usage: ./run_frame.sh [CERTS_ROOT] [STATE_DIR]   (defaults: ./certs, ./work/state_s2026). A C compiler and
# Python 3.9+. The base histogram enumerates about 7.4e11 pairs: about three hours on two threads and under 2 GB;
# a stopped run resumes from STATE_DIR.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
state=${2:-$here/work/state_s2026}
mkdir -p "$here/results" "$here/bin" "$state"
cc -O2 -o "$here/bin/c45_frame" "$here/c45_frame.c" -lpthread || exit 1
status=0
python3 -I "$here/c45_frame_check.py" selftest "$here/bin" > "$here/results/frame_selftest.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/frame_selftest.txt"
(cd "$here" && python3 -I c45_frame_check.py s2026 "$here/bin" "$state" --threads 2 --passes 4 --certs "$root" --controls \
   > results/s2026.txt 2>&1) || status=1
grep -H "RESULT" "$here/results/s2026.txt"
exit $status
