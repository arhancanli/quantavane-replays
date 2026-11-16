#!/bin/sh
# 5a [St26] C_5a <= 0.98173* and [CHO25] C_5a <= 0.98183: the bound b of Carter, Hunter and O'Bryant's Lemmas 3.1
# and 3.2 (arXiv 2310.20032, Section 3), maximised exactly over every w, with each row's own parameters (St26: the
# six-decimal parameters of its certified record; CHO25: the proof of Theorem 3.3), by two separately written exact
# programs (c5a_st26_a.py: primal LPs on cells of pair-difference levels; c5a_st26_b.py: dual LPs on the paper's
# interlacing cells), then red controls that both must refuse. C_5a <= b/2 (Proposition 1.3).
# Nothing is downloaded; no program or LP certificate of either source was read or run.
# Usage: ./run_all.sh   Python 3.9+, standard library only; one process at a time, about 25 minutes.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c5ast26.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
run() {  # run NAME ARGS...: A then B, one at a time
  name=$1; shift
  python3 -I "$here/c5a_st26_a.py" "$@" > "$tmp/$name.a" 2>&1
  python3 -I "$here/c5a_st26_b.py" "$@" > "$tmp/$name.b" 2>&1
}
run st26 st26
run cho25 cho25
run x1 st26 --claim 0.98172
run x2 st26 --factor1
run x3 st26 --drop 0,1
frac() { grep -h "^b\* = [0-9]*/[0-9]*$" "$1" | head -1; }
{
  echo "== row 5a [St26] 0.98173*, program A"; cat "$tmp/st26.a"
  echo "== row 5a [St26] 0.98173*, program B"; cat "$tmp/st26.b"
  echo "== row 5a [CHO25] 0.98183, program A"; cat "$tmp/cho25.a"
  echo "== row 5a [CHO25] 0.98183, program B"; cat "$tmp/cho25.b"
  echo "== red control 1: the St26 claim lowered to 0.98172"; grep -h "program\|^b\* = [0-9]\.\|row \|VERDICT" "$tmp/x1.a" "$tmp/x1.b"
  echo "== red control 2: Lemma 3.2 with tau/c^2 in place of 2 tau/c^2"; grep -h "program\|^b\* = [0-9]\.\|row \|VERDICT" "$tmp/x2.a" "$tmp/x2.b"
  echo "== red control 3: the two near-equal windows c_1 = 0.664209 and c_2 = 0.664134 removed (only c_3 kept)"; grep -h "program\|^b\* = [0-9]\.\|row \|VERDICT" "$tmp/x3.a" "$tmp/x3.b"
  echo "== agreement: the exact maxima b*, A against B"
  for n in st26 cho25 x2 x3; do
    fa=$(frac "$tmp/$n.a"); fb=$(frac "$tmp/$n.b")
    if [ -n "$fa" ] && [ "$fa" = "$fb" ]; then echo "$n: equal"; else echo "$n: DIFFERENT ($fa | $fb)"; fi
  done
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error\|VERDICT FAILED\|NOT REFUSED\|SELF-TEST FAILED" "$out")
agree=$(grep -c ": equal$" "$out")
if [ "$agree" = 4 ] && [ "$pass" = 4 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs, the same exact maxima: St26 b*/2 <= 0.98173 and CHO25 b*/2 <= 0.98183; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
