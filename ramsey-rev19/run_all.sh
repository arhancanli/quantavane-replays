#!/bin/sh
# Check every witness fetched by fetch.sh with two exact programs, under nice 10, each run capped in CPU seconds.
# Usage: ./run_all.sh [CERTS_ROOT] [all] [only=NAME]   (default: ./certs; Python 3.9+, a C compiler)
#   Without "all", only the witnesses marked "cheap" run (under an hour of CPU on an 18-core Mac); with "all" the
#   others run too, uncapped (many core-days, see README). One log per run in results/, red controls included.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
tier=cheap; only=
for a in "$@"; do case $a in all) tier=all ;; only=*) only=${a#only=} ;; esac; done
bin=$root/build; mat=$root/matrices; ctl=$root/controls; res=$here/results
mkdir -p "$bin" "$mat/nbhd" "$mat/relabelled" "$ctl" "$res"
cc -O2 -o "$bin/ramsey_check" "$here/checkers/ramsey_check.c" || exit 2
cc -O2 -o "$bin/ramsey_indep" "$here/checkers/ramsey_indep.c" || exit 2
# the same source with room for 2048 vertices, for the two witnesses above 1024 (check_one.sh picks it)
sed -e 's/^#define MAXN 1024$/#define MAXN 2048/' -e 's/^#define W 16$/#define W 32/' "$here/checkers/ramsey_indep.c" \
    > "$bin/ramsey_indep_2048.c" && cc -O2 -o "$bin/ramsey_indep_2048" "$bin/ramsey_indep_2048.c" || exit 2
py() { python3 -I "$@"; }
run() { "$here/check_one.sh" "$bin" "$res" "$@"; }
# name | format | file | n | complement? | s | t | extra | tier | CPU cap (s) for the full check
while IFS='|' read -r name fmt file n comp s t extra cls cap; do
  case $name in ''|'#'*) continue ;; esac
  [ -n "$only" ] && [ "$only" != "$name" ] && continue
  case $cls in cheap*) ;; *) [ "$tier" = cheap ] && [ -z "$only" ] && continue ;; esac
  [ "$tier" = all ] && cap=100000000
  label=; [ "${fmt#circulant:}" != "$fmt" ] && label="--label ${fmt#circulant:}" && fmt=circulant
  flag=; [ "$comp" = yes ] && flag=--complement
  py "$here/to_matrix.py" $fmt "$root/$file" "$mat/$name.txt" --n $n $flag $label > "$res/$name.convert.txt" || continue
  m=$mat/$name.txt
  case $extra in
    nbhd)  # circulant R(3,t): both programs from vertex 0 (exact by rotation); the plain full check is too slow
      py "$here/nbhd.py" "$m" 1 "$mat/nbhd/$name.N1.txt" >> "$res/$name.convert.txt"
      py "$here/nbhd.py" "$m" 0 "$mat/nbhd/$name.N0.txt" >> "$res/$name.convert.txt"
      run $name.nbhd1 300 check $(( s - 1 )) $(( n + 1 )) "$mat/nbhd/$name.N1.txt"
      run $name.nbhd0 $cap check $(( n + 1 )) $(( t - 1 )) "$mat/nbhd/$name.N0.txt"
      run $name $cap indep "$m" ;;
    relabel)  # circulant: ramsey_indep from vertex 0, then again on a renumbered copy with no symmetry
      run $name $cap check $s $t "$m"
      run $name $cap indep "$m"
      py "$here/relabel.py" "$m" "$mat/relabelled/$name.perm.txt" >> "$res/$name.convert.txt"
      run $name.perm $cap indep "$mat/relabelled/$name.perm.txt" ;;
    *)
      run $name $cap check $s $t "$m"
      run $name $cap indep "$m" ;;
  esac
  if [ "$name" = yu_R4_20_n251 ]; then  # Yu's second graph, x -> 116x, as section 4 of the paper builds it
    py "$here/yu_gprime.py" "$m" "$mat/${name}_Gprime.txt" > "$res/${name}_Gprime.convert.txt"
    run ${name}_Gprime $cap check $s $t "$mat/${name}_Gprime.txt"
    run ${name}_Gprime $cap indep "$mat/${name}_Gprime.txt"
  fi
  case $extra$cls in *ctl*)
    py "$here/controls.py" "$m" $s $t "$ctl/$name.ctrl_red.txt" "$ctl/$name.ctrl_blue.txt" > "$res/$name.controls.txt"
    for c in red blue; do
      run $name.ctrl_$c 900 check $s $t "$ctl/$name.ctrl_$c.txt"
      run $name.ctrl_$c 900 indep "$ctl/$name.ctrl_$c.txt"
    done ;;
  esac
done <<'EOF'
# the bounds since revision 18; "ctl" adds the two one-edge red controls
nart_R3_17_n92|numpy|nart/R3_17_93.txt|92|no|3|17|-|cheap,ctl|600
nart_R4_15_n162|numpy|nart/R4_15_163.txt|162|no|4|15|-|cheap,ctl|600
nart_R4_16_n173|numpy|nart/R4_16_174.txt|173|no|4|16|relabel|cheap,ctl|600
nart_R4_18_n208|numpy|nart/R4_18_209.txt|208|no|4|18|-|cheap|600
nart_R4_19_n218|numpy|nart/R4_19_219.txt|218|no|4|19|relabel|cheap|600
nart_R4_20_n236|numpy|nart/R4_20_237.txt|236|no|4|20|relabel|cheap|600
nart_R4_21_n243|numpy|nart/R4_21_244.txt|243|no|4|21|-|cheap|600
yu_R4_20_n251|circulant:S02|yu/CERTIFICATE_R4_20.txt|251|no|4|20|relabel|cheap,ctl|900
cremer_R6_8_n134|rows|cremer/n134-gpt-6.1-sol.txt|134|no|6|8|-|cheap,ctl|600
cremer_R8_10_n344_1|rows|cremer/n344-gpt-6.1-sol-1.txt|344|no|8|10|-|cheap,ctl|600
cremer_R8_10_n344_2|rows|cremer/n344-gpt-6.1-sol-2.txt|344|no|8|10|-|cheap|600
qv_R5_15_n278|rows|quantavane/R5_15_K278.txt|278|no|5|15|relabel|cheap,ctl|900
coniglio_R3_24_n150|dimacs|coniglio/R_3_24_graphSize_150.clq|150|yes|3|24|nbhd|cheap|1800
coniglio_R3_25_n159|dimacs|coniglio/R_3_25_graphSize_159.clq|159|yes|3|25|nbhd|cheap|1800
coniglio_R3_26_n162|dimacs|coniglio/R_3_26_graphSize_162.clq|162|yes|3|26|nbhd|cheap|1800
coniglio_R3_28_n180|dimacs|coniglio/R_3_28_graphSize_180.clq|180|yes|3|28|nbhd|cheap|1800
coniglio_R3_29_n193|dimacs|coniglio/R_3_29_graphSize_193.clq|193|yes|3|29|nbhd|cheap|1800
coniglio_R3_30_n198|dimacs|coniglio/R_3_30_graphSize_198.clq|198|yes|3|30|nbhd|cheap|1800
coniglio_R3_31_n209|dimacs|coniglio/R_3_31_graphSize_209.clq|209|yes|3|31|nbhd|cheap|1800
coniglio_R3_32_n218|dimacs|coniglio/R_3_32_graphSize_218.clq|218|yes|3|32|nbhd|cheap|1800
coniglio_R3_33_n229|dimacs|coniglio/R_3_33_graphSize_229.clq|229|yes|3|33|nbhd|all|1800
coniglio_R3_34_n238|dimacs|coniglio/R_3_34_graphSize_238.clq|238|yes|3|34|nbhd|all|1800
coniglio_R3_35_n252|dimacs|coniglio/R_3_35_graphSize_252.clq|252|yes|3|35|nbhd|all|1800
coniglio_R3_36_n260|dimacs|coniglio/R_3_36_graphSize_260.clq|260|yes|3|36|nbhd|all|1800
coniglio_R3_37_n271|dimacs|coniglio/R_3_37_graphSize_271.clq|271|yes|3|37|nbhd|all|1800
coniglio_R3_38_n280|dimacs|coniglio/R_3_38_graphSize_280.clq|280|yes|3|38|nbhd|all|1800
coniglio_R3_39_n293|dimacs|coniglio/R_3_39_graphSize_293.clq|293|yes|3|39|nbhd|all|1800
coniglio_R3_40_n304|dimacs|coniglio/R_3_40_graphSize_304.clq|304|yes|3|40|nbhd|all|1800
coniglio_R3_41_n316|dimacs|coniglio/R_3_41_graphSize_316.clq|316|yes|3|41|nbhd|all|1800
coniglio_R3_42_n326|dimacs|coniglio/R_3_42_graphSize_326.clq|326|yes|3|42|nbhd|all|1800
coniglio_R3_43_n339|dimacs|coniglio/R_3_43_graphSize_339.clq|339|yes|3|43|nbhd|all|1800
coniglio_R3_44_n347|dimacs|coniglio/R_3_44_graphSize_347.clq|347|yes|3|44|nbhd|all|1800
coniglio_R3_45_n363|dimacs|coniglio/R_3_45_graphSize_363.clq|363|yes|3|45|nbhd|all|1800
coniglio_R3_46_n370|dimacs|coniglio/R_3_46_graphSize_370.clq|370|yes|3|46|nbhd|all|1800
coniglio_R3_47_n384|dimacs|coniglio/R_3_47_graphSize_384.clq|384|yes|3|47|nbhd|all|1800
coniglio_R3_48_n390|dimacs|coniglio/R_3_48_graphSize_390.clq|390|yes|3|48|nbhd|all|1800
coniglio_R3_49_n410|dimacs|coniglio/R_3_49_graphSize_410.clq|410|yes|3|49|nbhd|all|1800
nart_R7_13_n512|numpy|nart/R7_13_513.txt|512|no|7|13|-|all|120
nart_R8_12_n510|numpy|nart/R8_12_511.txt|510|no|8|12|-|all|120
nart_R8_18_n1051|numpy|nart/R8_18_1052.txt|1051|no|8|18|-|all|120
nart_R9_11_n616|numpy|nart/R9_11_617.txt|616|no|9|11|-|all|120
nart_R10_12_n855|numpy|nart/R10_12_856.txt|855|no|10|12|-|all|120
tatarevic_R12_12_n1640|rows|tatarevic/r1212_1640.graph|1640|no|12|12|-|all|120
qv_R8_12_n518|rows|quantavane/R8_12_K518.txt|518|no|8|12|-|all|100000
EOF
grep -h "^# summary" "$res"/*.log | sort
