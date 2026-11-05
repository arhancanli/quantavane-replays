#!/bin/sh
# Replay the asterisked 82a upper bounds ([Pin26], [Epi26]; PR #242) and the open PRs #280 and #301 from the files
# fetched by fetch.sh; one output file per row in results/, red controls included.
# Usage: ./run_all.sh [CERTS_ROOT] [--workers=N]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
shift 2>/dev/null
mkdir -p "$here/results"
status=0
python3 -I "$here/c82_check.py" selftest "$@" > "$here/results/selftest.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/selftest.txt" | tr '\n' ' '; echo
for row in pin epi doche pr280 pr280b pr280c pr280d pr301 pr301a pr330 pr330b pr330c; do   # pin: every red control; the others: R1-R4
  python3 -I "$here/c82_check.py" $row "$root" --controls "$@" > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
exit $status
