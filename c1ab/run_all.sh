#!/bin/sh
# Replay eleven 1a and 1b rows from the heights fetched by fetch.sh with ../autocorr_check.py (exact integers), with
# red controls; one output file per row in results/. Both readings of decimal text are run: as binary64 (what a
# float parser gives) and as the printed decimals (--decimal). --rescale scales a 1b function whose heights do not
# sum to n/2 exactly down or up to n/2 (and is refused if a height then leaves [0, 1]).
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about a minute.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1   # paths in the output stay relative to this folder
root=${1:-certs}
mkdir -p results
C="python3 -I ../autocorr_check.py"
status=0
run() {  # kind file printed extra-flags...
  k=$1; f=$2; pv=$3; shift 3
  out="results/${f%.json}.txt"
  $C $k "$root/$f" --printed=$pv --controls "$@" > "$out" 2>&1 || status=1
  $C $k "$root/$f" --printed=$pv --decimal "$@" > "${out%.txt}.decimal.txt" 2>&1 || status=1
}
run_decimal() {  # a paper's printed decimals: the decimal reading only (as binary64 the heights miss n/2)
  k=$1; f=$2; pv=$3; shift 3
  $C $k "$root/$f" --printed=$pv --decimal --controls "$@" > "results/${f%.json}.decimal.txt" 2>&1 || status=1
}
run 1a c1a_alphaevolve_2025.json 1.5053
run 1a c1a_alphaevolve_dec2025.json 1.503164
run 1a tog_c1a_ours_2026.json 1.502862
run 1a ttt_ac1_sequence.json 1.502870
run 1a c1a_ea2026_2603.json 1.5027435320
run 1b c1b_alphaevolve_2025.json 0.380924 --rescale
run 1b ttt_erdos_sequence.json 0.380876 --rescale
run 1b c1b_together_ai_2026.json 0.380871 --rescale
run 1b c1b_ea2026_2507.json 0.3808585749 --rescale
run 1b atris_c1b_construction.json 0.3808585504
run_decimal 1b c1b_haugland_2016.json 0.380926
grep -H "printed\|normalisation\|Error" results/*.txt
exit $status
