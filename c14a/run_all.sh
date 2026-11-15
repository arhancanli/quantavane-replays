#!/bin/sh
# 14a: the determinations behind the rows C_14 >= 4 [LS1965] (BB(3)) and C_14 >= 5 [B1983] (BB(4)), by two separately
# written programs that enumerate every 2-, 3- and 4-state machine in tree normal form and decide each one
# (c14a_a.c: re-simulation, then closed sets of local views only; c14a_b.c: resumed runs, then cyclers, translated
# cyclers and its own closed sets of local views), then red controls that both must refuse. Nothing is downloaded.
# Usage: ./run_all.sh   a C compiler; about 15 minutes on two cores, under 600 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c14a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
cc -O2 -o "$tmp/a" "$here/c14a_a.c" || exit 1
cc -O2 -o "$tmp/b" "$here/c14a_b.c" || exit 1
pair() {  # pair NAME ARGS-A -- ARGS-B
  name=$1; shift
  aargs=""; while [ "$1" != "--" ]; do aargs="$aargs $1"; shift; done; shift
  "$tmp/a" $aargs > "$tmp/$name.a" 2>&1 &
  "$tmp/b" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair n2 2 2000 -- 2 2000
pair n3 3 2000 -- 3 2000
pair n4 4 2000 -- 4 2000
pair x1 3 20 -- 3 20
pair x2 4 100 -- 4 100
pair x3 3 2000 --maxrad 1 --hdepth 0 -- 3 2000 --no-translated --hmax 0 --rmax 1
pair x4 3 2000 --selftest -- 3 2000 --selftest
{
  for k in n2 n3 n4; do
    echo "== $k: program A (c14a_a.c)"; grep -v "^HOLDOUT" "$tmp/$k.a"
    echo "== $k: program B (c14a_b.c)"; grep -v "^UNDECIDED" "$tmp/$k.b"
  done
  echo "== n4: the machines neither program decides (A, then B; sorted)"
  grep "^HOLDOUT" "$tmp/n4.a" | awk '{print $2}' | sort > "$tmp/ha"; grep "^UNDECIDED" "$tmp/n4.b" | awk '{print $2}' | sort > "$tmp/hb"
  cat "$tmp/ha"
  if cmp -s "$tmp/ha" "$tmp/hb"; then echo "the same $(wc -l < "$tmp/ha" | tr -d ' ') machines in both programs"; else echo "DIFFERENT undecided sets"; fi
  echo "== red control 1: n = 3 with a step limit of 20, so the 21-step machines must be left undecided"
  grep -h "HOLDOUT\|UNDECIDED\|RESULT" "$tmp/x1.a" "$tmp/x1.b"
  echo "== red control 2: n = 4 with a step limit of 100; the 107-step champion 1RB1LB_1LA0LC_---1LD_1RD0RA must be left undecided"
  grep -h "1RB1LB_1LA0LC_---1LD_1RD0RA\|RESULT" "$tmp/x2.a" "$tmp/x2.b"
  echo "== red control 3: n = 3 with the tests cut down (radius 1, no history; B also without translated cyclers)"
  grep -h "RESULT" "$tmp/x3.a" "$tmp/x3.b"
  echo "== red control 4: every closed-set test (and in B the cycler tests) applied to every halting 3-state machine"
  grep -h "SELFTEST" "$tmp/x4.a" "$tmp/x4.b"
} > "$out" 2>&1
rm -rf "$tmp"
ok=1
[ "$(grep -c 'RESULT n = 2: S = 6, Sigma = 4 (every machine decided)' "$out")" = 2 ] || ok=0
[ "$(grep -c 'RESULT n = 3: S = 21, Sigma = 6 (every machine decided)' "$out")" = 2 ] || ok=0
[ "$(grep -c 'RESULT n = 4: undecided (12 ' "$out")" = 2 ] || ok=0
[ "$(grep -c 'S >= 107, Sigma >= 13' "$out")" = 2 ] || ok=0
grep -q "the same 12 machines in both programs" "$out" || ok=0
[ "$(grep -c '^HOLDOUT 1RB---_1LB0RC_1LC1LA\|^UNDECIDED 1RB---_1LB0RC_1LC1LA' "$out")" = 2 ] || ok=0
[ "$(grep -c 'RESULT n = 3: undecided (1 ' "$out")" = 2 ] || ok=0
[ "$(grep -c '^HOLDOUT 1RB1LB_1LA0LC_---1LD_1RD0RA\|^UNDECIDED 1RB1LB_1LA0LC_---1LD_1RD0RA' "$out")" = 2 ] || ok=0
[ "$(grep -c 'RESULT n = 3: undecided' "$out")" = 4 ] || ok=0
[ "$(grep -c 'halting machines: 0 false proofs' "$out")" = 2 ] || ok=0
if [ "$ok" = 1 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: BB(2) S = 6, Sigma = 4; BB(3) S = 21, Sigma = 6 with every 3-state machine decided; BB(4): S >= 107, Sigma >= 13 attained, all but the same 12 of 858,909 machines decided; four red controls as required)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED" >> "$out"; tail -1 "$out"; exit 1
fi
