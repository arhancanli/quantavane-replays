#!/bin/sh
# Replay of three C45 upper bounds (Romanoff's constant) from the 45a page's obstruction inequality and the prime sets
# alone: [G2026] (36 primes, PR #82), [Y2026] (38 primes, PR #125, the page's record) and the 21-base-prime row [P2026]
# of PR #281 (49 primes). Our own histogram of nu0 is checked against brute force on small sets (selftest), against
# two exact moments, and against the packet's histograms (fetched by fetch.sh, read as data only).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). A C compiler and Python 3.9+; the 21-prime row takes about
# 13 minutes and 0.6 GB, the others seconds.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results" "$here/bin"
cc -O2 -o "$here/bin/c45_comb" "$here/c45_comb.c" || exit 1
status=0
python3 -I "$here/c45_check.py" selftest "$here/bin" > "$here/results/selftest.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/selftest.txt"
for row in griego yoo pin21; do
  python3 -I "$here/c45_check.py" $row "$here/bin" "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT" "$here/results/$row.txt"
done
exit $status
