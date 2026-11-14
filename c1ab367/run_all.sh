#!/bin/sh
# Replay the step-function certificates of PR #367 (1a: the row's file and the file its description names; the
# [Cum2026] file it records) and PR #368 (1b) with two separately written exact programs (A: Kronecker substitution in
# Python integers; B: plain double loops in C, its own reader), then red controls that both must refuse.
# No program of any source was read or run. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs)
# Python 3.9+ and a C compiler; about five minutes on two cores, under 300 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/c1ab367.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
F367="$root/pr367/C1_262144_1.5026648332.txt"
F367D="$root/pr367d/C1_262144_1.5026648336.txt"
FCUM="$root/cum2026/49c83eb6f9f82e244578fab15d056757886679700bd8c8c653f0b100ae882f2d"
F368="$root/pr368/C1b_15360_0.3808564712.json"
C367=735013322210029338350039415555883008/489139897323017668625465543540606041
C367D=11484584234108926456768515550896128/7642811608550585358028269349510225
C368=3707849734734953936493216938911489/9735556609752801803494680617287680
pair() {  # pair NAME KIND FILE [ARGS...]: both programs at once, outputs kept apart
  name=$1; shift
  python3 -I "$here/c1ab367_kron.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c1ab367_direct.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair r367 1a "$F367" --claim $C367 --argmax 344084 --printed 1.5026648332
pair r367d 1a "$F367D" --claim $C367D --argmax 316694
pair rcum 1a "$FCUM" --prefix 1.50266484077609912906120 --printed 1.50266484078
pair r368 1b "$F368" --claim $C368 --argmax 6384 --printed 0.3808564712
pair i367d 1a "$F367D" --printed 1.5026648336
pair x1 1a "$F367" --bump 1000 1 --claim $C367
pair x2 1a "$F367" --printed 1.5026648331
pair x3 1a "$F367" --bump 5 -2000000000000000
pair x4 1a "$FCUM" --bump 1000 1000000000000000 --prefix 1.50266484077609912906120
pair x5 1a "$FCUM" --printed 1.50266484077
pair x6 1b "$F368" --bump 2490 1
pair x7 1b "$F368" --move 5753 11656 1000 --claim $C368
pair x8 1b "$F368" --printed 0.3808564711
sel() { grep -h -v "run time" "$@"; }
{
  echo "== PR #367 (1a), the row in its files: C1_262144_1.5026648332.txt (intervalplace/frontier @ 9edc297)"; sel "$tmp/r367.a" "$tmp/r367.b"
  echo "== PR #367 (1a), the certificate its description and title name: C1_262144_1.5026648336.txt (@ 8c9fa72)"; sel "$tmp/r367d.a" "$tmp/r367d.b"
  echo "== [Cum2026] (1a), recorded by PR #367: the ancillary file 49c83eb6...2f2d"; sel "$tmp/rcum.a" "$tmp/rcum.b"
  echo "== PR #368 (1b): C1b_15360_0.3808564712.json (intervalplace/frontier @ 6d6c4ba)"; sel "$tmp/r368.a" "$tmp/r368.b"
  echo "== information: the description's value 1.5026648336 compared with the exact value of its file (not a red control)"
  grep "printed upper bound\|VERDICT" "$tmp/i367d.a" "$tmp/i367d.b" | sed 's/^.*\/i367d\.[ab]://'
  for x in x1 x2 x3 x4 x5 x6 x7 x8; do
    case $x in
      x1) echo "== red control 1: PR #367 row, height 1000 plus one unit";;
      x2) echo "== red control 2: PR #367 row, claimed as an upper bound 1.5026648331 (below the exact value)";;
      x3) echo "== red control 3: PR #367 row, height 5 made negative";;
      x4) echo "== red control 4: [Cum2026], height 1000 plus 10^15 units, against the stated 23 digits";;
      x5) echo "== red control 5: [Cum2026], claimed as an upper bound 1.50266484077 (below the exact value)";;
      x6) echo "== red control 6: PR #368, height 2490 (equal to Q) plus one unit";;
      x7) echo "== red control 7: PR #368, 1000/Q moved from height 11656 to height 5753 (mass unchanged)";;
      x8) echo "== red control 8: PR #368, claimed as an upper bound 0.3808564711 (below the exact value)";;
    esac
    grep "CONTROL\|FAIL\|exact value\|VERDICT" "$tmp/$x.a" "$tmp/$x.b" | sed "s/^.*\/$x\.[ab]://"
  done
} > "$out" 2>&1
fp_ok=1
for r in r367 r367d rcum r368; do
  fa=$(grep "fingerprint" "$tmp/$r.a" | awk '{print $NF}'); fb=$(grep "fingerprint" "$tmp/$r.b" | awk '{print $NF}')
  va=$(grep "exact value" "$tmp/$r.a" | awk '{print $NF}'); vb=$(grep "exact value" "$tmp/$r.b" | awk '{print $NF}')
  if [ -z "$fa" ] || [ "$fa" != "$fb" ] || [ -z "$va" ] || [ "$va" != "$vb" ]; then fp_ok=0; fi
done
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback" "$out")
# expected: 8 PASS (four certificates, two programs); the information run refused by both (2); 16 refusals in controls
if [ "$fp_ok" = 1 ] && [ "$pass" = 8 ] && [ "$ref" = 18 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs give the same exact value and the same fingerprint of every sequence for all four certificates; eight red controls refused by both; the description's 1.5026648336 is below the exact value of its file, as recorded above)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (pass $pass, refused $ref, fingerprints agree $fp_ok)" >> "$out"; tail -1 "$out"; exit 1
