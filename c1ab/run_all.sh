#!/bin/sh
# Replay fourteen 1a and 1b rows from the heights fetched by fetch.sh with ../autocorr_check.py (exact integers), with
# red controls; one output file per row in results/. Both readings of decimal text are run: as binary64 (what a
# float parser gives) and as the printed decimals (--decimal). --rescale scales a 1b function whose heights do not
# sum to n/2 exactly down or up to n/2 (and is refused if a height then leaves [0, 1]).
# Expected outcomes: every run exits 0, except the three whose printed value the row's own heights do not support,
# which must exit 1 and print "NOT a valid ... bound as printed": tog_c1a_ours_2026 (both readings) and
# c1b_haugland_2016 (the decimal reading). The script exits 1 only on an unexpected outcome and names the row and why.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about a minute.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1   # paths in the output stay relative to this folder
root=${1:-certs}
mkdir -p results
C="python3 -I ../autocorr_check.py"
status=0
expected() {  # output file -> the exit status its run must have
  case "$1" in
    results/tog_c1a_ours_2026.txt|results/tog_c1a_ours_2026.decimal.txt|results/c1b_haugland_2016.decimal.txt) echo 1 ;;
    *) echo 0 ;;
  esac
}
one() {  # output-file command...: run it, then compare its exit status with the table above
  out=$1; shift
  "$@" > "$out" 2>&1
  rc=$?
  want=$(expected "$out")
  if [ "$rc" -ne "$want" ]; then
    echo "UNEXPECTED: $out exited $rc, expected $want" >&2
    status=1
  elif [ "$want" -eq 1 ] && ! grep -q "NOT a valid .* bound as printed" "$out"; then
    echo "UNEXPECTED: $out exited 1 without a 'NOT a valid ... bound as printed' line" >&2
    status=1
  fi
}
run() {  # kind file printed extra-flags...
  k=$1; f=$2; pv=$3; shift 3
  out="results/${f%.json}.txt"
  one "$out" $C $k "$root/$f" --printed=$pv --controls "$@"
  one "${out%.txt}.decimal.txt" $C $k "$root/$f" --printed=$pv --decimal "$@"
}
run_decimal() {  # a paper's printed decimals: the decimal reading only (as binary64 the heights miss n/2)
  k=$1; f=$2; pv=$3; shift 3
  one "results/${f%.json}.decimal.txt" $C $k "$root/$f" --printed=$pv --decimal --controls "$@"
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
# added 10 Nov 2026: SimpleTES (1a and 1b) and ThetaEvolve (1a)
run 1a c1a_simpletes_2026.json 1.503871
run 1b c1b_simpletes_2026.json 0.380868 --rescale
run 1a c1a_thetaevolve_2025.json 1.503133
grep -H "printed\|normalisation\|Error" results/*.txt
# the 1a lower bound [PBV2026] 1.292: the five anchors and the strict-failure inequality (about three minutes)
one results/pbv2026.txt python3 -I pbv_check.py "$root/pbv2026_multiscale_arcsine_1292.json" --controls
grep -H "RESULT" results/pbv2026.txt
exit $status
