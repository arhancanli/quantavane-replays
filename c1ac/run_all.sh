#!/bin/sh
# Replay the six certificates fetched by fetch.sh with ../autocorr_check.py (exact integers, red controls), and check
# each result against what its pull request states: the reduced fraction, the deciding lag and the printed bound.
# Red controls in every run: a corrupted coefficient and a changed weight are refused by the polynomial identity,
# scaling by 3 leaves the bound unchanged, and the printed bound moved one unit in its last digit past the exact value
# (0.4103768837, 0.4103604996, 0.4103608340, 1.5027427355, 1.4546733457, 1.5027427340) is refused; that last control makes the checker exit 1,
# so each run must exit 1 with exactly that one "NOT a valid" line. One output file per certificate in results/.
# Exits 0 only when every expected line is printed and nothing else fails.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about five minutes.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1   # paths in the output stay relative to this folder
root=${1:-certs}
mkdir -p results
C="python3 -I ../autocorr_check.py"
status=0
run() {  # kind file printed fraction lag-line moved-printed
  k=$1; f=$2; pv=$3; frac=$4; lag=$5; bad=$6
  out="results/${f%.json}.txt"
  $C $k "$root/$f" --printed=$pv --printed=$bad --controls > "$out" 2>&1
  [ $? -eq 1 ] || { echo "FAILED: $out did not exit 1 (the moved printed value must be refused)" >&2; status=1; }
  [ "$(grep -c 'NOT a valid' "$out")" -eq 1 ] || { echo "FAILED: $out: expected exactly one refused printed value" >&2; status=1; }
  [ "$(grep -c 'control: .*: True$' "$out")" -eq 3 ] || { echo "FAILED: $out: a red control did not hold" >&2; status=1; }
  for want in "bound = $frac" "$lag" "printed $pv: " "printed $bad: " "NOT a valid"; do
    grep -q -F "$want" "$out" || { echo "FAILED: $out lacks '$want'" >&2; status=1; }
  done
  grep -F "printed $pv: " "$out" | grep -q -v "NOT a valid" || { echo "FAILED: $out: the stated bound $pv is refused" >&2; status=1; }
  grep -F "printed $bad: " "$out" | grep -q "NOT a valid" || { echo "FAILED: $out: the moved value $bad was not refused" >&2; status=1; }
}
# #354, 1c >= 0.4103768: the fraction, minimum at k = 46008 (S = 1509900967, min A = 3568945316176), as the PR states
run 1c sean_c1c_K262144.json 0.4103768 935577600963641344/2279800930147535089 "argmin_k: 46008" 0.4103768837
# #270, 9 Nov 13:22 UTC update (c1c_v4): 1c >= 0.4103604995, minimum only at k = 98288
run 1c vidigoat_c1c_K98304_v4.json 0.4103604995 5282641762424933775753216/12873173145647069582596129 "argmin_k: 98288" 0.4103604996
# #270, 9 Nov 16:56 UTC update (c1c_v5): 1c >= 0.4103608339, minimum at k = 45081
run 1c vidigoat_c1c_K98304_v5.json 0.4103608339 2556315663995378551725293568/6229433836731244234885554241 "argmin_k: 45081" 0.4103608340
# #270, 9 Nov 16:56 UTC update (c1a_v2): 1a <= 1.5027427356, peak at index 96751
run 1a vidigoat_c1a_n262144_v2.json 1.5027427356 29134419116089766846043850668467284344832/19387496227028661445940960928888412813689 "argmax_m: 96751" 1.5027427355
# #270, 10 Nov update (eb710a8, c1d): 1d <= 1.4546733458, peak at index 165229
run 1d vidigoat_c1d_n163840.json 1.4546733458 646233848324886328286774558720/444246710235775351931155445761 "argmax_m: 165229" 1.4546733457
# #270, 10 Nov update (eb710a8, c1a_v3): 1a <= 1.5027427341, peak at index 262143
run 1a vidigoat_c1a_n262144_v3.json 1.5027427341 116537598371962553887369519059399403896832/77549933019818026777245548252121251545225 "argmax_m: 262143" 1.5027427340
grep -H "bound =\|printed\|control" results/*.txt
exit $status
