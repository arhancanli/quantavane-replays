#!/bin/sh
# Replay two step-function rows with ../autocorr_check.py (exact integers), each in both readings of the printed
# decimals (as binary64 and, with --decimal, as printed), with red controls; one output file per run in results/:
#   1a [MV2009] C_1a <= 1.50992: Matolcsi and Vinuesa's 208-step function (arXiv 0907.1379v2, appendix);
#   1d [GGSWT2025] C_1d <= 1.4688: AlphaEvolve's signed 150-step function (alphaevolve_results, section B.3, cell 65).
# Then signed_note.py: Matolcsi and Vinuesa's signed 150-step example, for the note on the 1d [MV2010] row.
# Every run must exit 0. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs; run ./fetch.sh first). Python 3.9+,
# standard library only; a few seconds.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1   # paths in the output stay relative to this folder
root=${1:-certs}
mkdir -p results
C="python3 -I ../autocorr_check.py"
status=0
run() {  # kind file printed
  out="results/${2%.json}.txt"
  $C $1 "$root/$2" --printed=$3 --controls > "$out" 2>&1 || { echo "UNEXPECTED: $out"; status=1; }
  $C $1 "$root/$2" --printed=$3 --decimal > "${out%.txt}.decimal.txt" 2>&1 || { echo "UNEXPECTED: ${out%.txt}.decimal.txt"; status=1; }
  grep -H "bound = \|printed" "$out" "${out%.txt}.decimal.txt"
}
run 1a mv2009_n208.json 1.50992
run 1d ae_c3prime.json 1.4688
python3 -I signed_note.py "$root/mv2009_n150s.json" > results/mv2009_n150s_note.txt 2>&1 || { echo "UNEXPECTED: note"; status=1; }
grep -H "RESULT" results/mv2009_n150s_note.txt
exit $status
