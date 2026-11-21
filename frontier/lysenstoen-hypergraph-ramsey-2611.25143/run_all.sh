#!/bin/sh
# Lysenstoen, arXiv 2611.25143: six hypergraph Ramsey lower bounds, checked by program A (C) and program B (JavaScript),
# then the red controls (each must be refused). Needs a C compiler, Python 3 and node >= 18. About two minutes.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results work
cc -O2 -o work/lys_a lys_a.c || exit 1
python3 -I extract_a.py certs work/a_in > /dev/null || exit 1
status=0
run() { name=$1; shift; "$@" > "results/$name.txt" 2>&1; rc=$?; echo "$name exit=$rc"; return $rc; }
run A_t37 work/lys_a tri 37 4 5 work/a_in/c37.txt || status=1
run A_t70 work/lys_a tri 70 4 6 work/a_in/c70.txt || status=1
run A_t98 work/lys_a tri 98 5 5 work/a_in/c98.txt || status=1
run A_f83 work/lys_a five work/a_in/s83.txt 6 || status=1
run A_p200 work/lys_a proj 199 work/a_in/t199.txt 5 6 psl || status=1
run A_p600 work/lys_a proj 599 work/a_in/t599.txt 6 6 pgl || status=1
for c in t37 t70 t98 f83 p200 p600; do run B_$c node lys_b.mjs certs $c || status=1; done
# red controls: a nonzero exit is the expected refusal
python3 -I plant_a.py work/a_in work/ctl > /dev/null || exit 1
refused() { name=$1; shift; if "$@" > "results/control_$name.txt" 2>&1; then echo "control $name NOT refused"; status=1; else echo "control $name refused"; fi; }
refused A_plant37 work/lys_a tri 37 4 5 work/ctl/c37_plant.txt
refused A_plant70 work/lys_a tri 70 4 6 work/ctl/c70_plant.txt
refused A_plant98 work/lys_a tri 98 5 5 work/ctl/c98_plant.txt
refused A_plant83 work/lys_a five work/ctl/s83_plant.txt 6
refused A_flip600 work/lys_a proj 599 work/ctl/t599_flip.txt 6 6 pgl
refused A_size4_200 work/lys_a proj 199 work/a_in/t199.txt 4 6 psl
for c in "plant37 t37" "plant70 t70" "plant98 t98" "plant83 f83" "flip600 p600" "flip200 p200"; do
  set -- $c; refused B_$1 env LYS_CONTROL=$1 node lys_b.mjs certs $2
done
[ $status = 0 ] && echo "ALL PASS (six colourings confirmed by both programs; every red control refused)" || echo FAIL
exit $status
