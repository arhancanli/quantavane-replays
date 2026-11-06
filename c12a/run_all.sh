#!/bin/sh
# Fetch the data if needed (fetch.sh), build c12a_check and replay the computation. For each spec file every stored
# number is recomputed and compared (results/SPEC.txt, progress in results/SPEC.log); the cell sums of tiny to k6 are
# recomputed again by all orders in every cell, with none of the identities (results/SPEC_orders.txt), and seven cell
# sums each of k7, k8 and k9 by Held-Karp in every cell (results/SPEC_hk.txt); then the red controls written by
# make_controls.sh (results/controls.txt). k9 takes hours (1 to 2 h on 12 threads of a loaded Apple M5 Pro):
# ./run_all.sh quick leaves it out. Threads: J (default 12). Needs cc (C11 with unsigned __int128), awk, sed, bc and
# python3 (for one relative path).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
[ -d certs/tsp ] || ./fetch.sh || exit 2
mkdir -p bin results
cc -std=c11 -O3 -pthread c12a_check.c -o bin/c12a_check -lm || exit 2
J=${J:-12}
sp=certs/tsp/data/specs
status=0
run() {  # run NAME SPEC [OPTION...]
  name=$1; spec=$2; shift 2
  ./bin/c12a_check "$spec" -j "$J" "$@" > "results/$name.txt" 2> "results/$name.log"
  st=$?
  printf '%-15s status %s  %s\n' "$name" "$st" "$(grep RESULT "results/$name.txt")"
  [ "$st" = 0 ] || status=1
}
specs="tiny rounding bhh k4 k5 k6 k7 k8"
[ "$1" = quick ] || specs="$specs k9"
for s in $specs; do run "$s" "$sp/$s.txt"; done
for s in tiny rounding bhh k4 k5 k6; do run "${s}_orders" "$sp/$s.txt" --orders; done
run k7_hk "$sp/k7.txt" --hk 0,1,4,4096,7000,12345,16383
run k8_hk "$sp/k8.txt" --hk 0,1,4,16384,30000,54321,65535
[ "$1" = quick ] || run k9_hk "$sp/k9.txt" --hk 0,1,4,65536,100000,199999,262143
./make_controls.sh || exit 2
ctl() {  # ctl SPEC EXPECTED_STATUS DESCRIPTION
  echo "== $3: $1, expected status $2"
  ./bin/c12a_check "$1" -j "$J" 2> /dev/null
  st=$?
  if [ "$st" = "$2" ]; then echo "status $st, as expected"; else echo "status $st, NOT AS EXPECTED"; status=1; fi
  echo
}
{
  ctl controls/k5_sum_plus_one.txt 1 "(a) k5 values file with the cell sum S_100 increased by 1"
  ctl controls/k7_threshold_below.txt 1 "(b) k7 spec with the threshold 10^-12 below the printed (rounded-up) bound"
  ctl controls/k7_threshold_at.txt 0 "(b, green) k7 spec with the threshold equal to the printed bound"
  ctl controls/k5_mesh127.txt 1 "(c) k5 spec and values file with mesh 127 in place of 128, nothing else changed"
  ctl controls/k5_mesh127_orig.txt 1 "(c') k5 spec with mesh 127 and the unaltered values file"
} > results/controls.txt
grep "^== \|^status" results/controls.txt
exit $status
