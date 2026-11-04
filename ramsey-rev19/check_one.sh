#!/bin/sh
# One exact check of one matrix, under nice 10, with its CPU time measured by /usr/bin/time and a CPU-seconds cap.
# Usage: check_one.sh BIN_DIR LOG_DIR NAME CAP_SECONDS check S T MATRIX   (engine/native/ramsey_check)
#        check_one.sh BIN_DIR LOG_DIR NAME CAP_SECONDS indep MATRIX      (independent/ramsey_indep)
# Writes LOG_DIR/NAME.<check|indep>.log: the command, the checker's output, exit status, /usr/bin/time -l report.
bin=$1; logs=$2; name=$3; cap=$4; which=$5; shift 5
case $which in
  check) cmd="$bin/ramsey_check $1 $2 $3"; m=$3 ;;
  indep) cmd="$bin/ramsey_indep $1"; m=$1
         n=$(awk 'END {print NR}' "$m"); [ "$n" -gt 1024 ] && cmd="$bin/ramsey_indep_2048 $1" ;;
  *) echo "usage"; exit 2 ;;
esac
log="$logs/$name.$which.log"
{
  echo "# $(date '+%Y-%m-%d %H:%M:%S %z') host $(hostname -s) load $(sysctl -n vm.loadavg)"
  echo "# matrix $m sha256 $(shasum -a 256 "$m" | cut -d' ' -f1)"
  echo "# command: nice -n 10 $cmd   (CPU cap ${cap} s)"
} > "$log"
( ulimit -t "$cap"; /usr/bin/time -l nice -n 10 $cmd ) >> "$log" 2> "$log.time"
status=$?
echo "# exit status $status" >> "$log"
cat "$log.time" >> "$log" && rm -f "$log.time"
user=$(awk '/ user /{for(i=1;i<=NF;i++) if($i=="user") print $(i-1)}' "$log")
sys=$(awk '/ sys/{for(i=1;i<=NF;i++) if($i=="sys") print $(i-1)}' "$log")
real=$(awk '/ real /{for(i=1;i<=NF;i++) if($i=="real") print $(i-1)}' "$log")
echo "# summary $name $which exit $status cpu_user $user cpu_sys $sys wall $real" >> "$log"
tail -1 "$log"
