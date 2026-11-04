#!/bin/sh
# Replay every asterisked 49a row from the files fetched by fetch.sh, plus the method self-test, the gluing-lemma
# instance tests and the DEGKM control; one output file per row in results/, red controls included.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results"
status=0
for row in selftest lemma irsi_v10 irsi_v11 irsi_v12 irsi_degkm cla_F cla_FH syl_unit syl_Bprime kos; do
  python3 -I "$here/sunflower_check.py" "$row" "$root" --controls > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
exit $status
