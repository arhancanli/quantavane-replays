#!/bin/sh
# One command per replay: fetch the certificate (when the folder downloads one) and run the folder's own run_all.sh.
#   ./check.sh --list          every replay folder, with the first line of its RESULT.md or of its run_all.sh
#   ./check.sh FOLDER [...]    for each folder: ./fetch.sh into FOLDER/certs if that is still empty, then ./run_all.sh;
#                              prints PASS or FAIL per folder (the exit status of run_all.sh) and ends with ALL PASS
# Most folders need only Python 3.9+; some also a C compiler. Run times are in each folder's RESULT.md or run_all.sh
# (from seconds to a few hours). Certificates are downloaded from the pinned commits and checked against SHA256SUMS.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
if [ "$1" = "--list" ] || [ -z "$1" ]; then
  for d in c*/ frontier/*/; do
    d=${d%/}; [ -f "$d/run_all.sh" ] || continue
    if [ -f "$d/RESULT.md" ]; then t=$(head -1 "$d/RESULT.md" | sed 's/^# *//'); else t=$(sed -n '2p' "$d/run_all.sh" | sed 's/^# *//'); fi
    printf '%-10s %s\n' "$d" "$t"
  done
  [ -z "$1" ] && echo "usage: ./check.sh FOLDER [...]   (./check.sh --list for the folders)"
  exit 0
fi
status=0
for d in "$@"; do
  d=${d%/}
  if [ ! -f "$d/run_all.sh" ]; then echo "FAIL $d: no such replay folder"; status=1; continue; fi
  if [ -f "$d/fetch.sh" ] && { [ ! -d "$d/certs" ] || [ -z "$(ls -A "$d/certs" 2>/dev/null)" ]; }; then
    (cd "$d" && sh ./fetch.sh > /dev/null 2>&1) || { echo "FAIL $d: fetch.sh failed (network, or a SHA-256 differs)"; status=1; continue; }
  fi
  if (cd "$d" && sh ./run_all.sh > /dev/null 2>&1); then echo "PASS $d"; else echo "FAIL $d (see $d/results/)"; status=1; fi
done
if [ $status = 0 ]; then echo "ALL PASS"; else echo "FAIL"; fi
exit $status
