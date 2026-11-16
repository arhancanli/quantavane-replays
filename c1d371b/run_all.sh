#!/bin/sh
# Replay the 1d step function of PR #371 (its eighth certificate, 15 Nov 19:46 UTC: 1.4546676818) with two separately written exact programs (A: one Kronecker squaring with
# balanced digits in Python integers; B: its own reader and plain loops over every knot in C, signed 128-bit
# integers), then red controls that both must refuse. No program of any source was read or run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs)
# Python 3.9+ and a C compiler; a few minutes on two cores, under 200 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/c1d371b.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
F="$root/pr371b/C1d_327680_1.4546676818.txt"
C=40122795831148524928230686720/27582104377369449232322657763
pair() {  # pair NAME [ARGS...]: both programs at once, outputs kept apart
  name=$1; shift
  python3 -I "$here/c1d371b_a.py" "$F" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c1d371b_b.py" "$F" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair r371b --claim $C --argmax 116123 --printed 1.4546676818
pair x1 --bump 50000 1 --claim $C
pair x2 --printed 1.4546676817
pair x3 --bump 0 -1942324222640 --claim $C
pair x4 --bump 327679 1000000 --claim $C
{
  echo "== PR #371 (1d), eighth certificate: C1d_327680_1.4546676818.txt (intervalplace/frontier @ 1ab20c7), program A then program B"
  grep -h -v "run time" "$tmp/r371b.a" "$tmp/r371b.b"
  for x in x1 x2 x3 x4; do
    case $x in
      x1) echo "== red control 1: height 50000 plus one unit, against the stated fraction";;
      x2) echo "== red control 2: claimed as an upper bound 1.4546676817 (below the exact value)";;
      x3) echo "== red control 3: height 0 (971162111320) negated, against the stated fraction";;
      x4) echo "== red control 4: the last height plus 10^6 units (S changes), against the stated fraction";;
    esac
    grep -h "CONTROL\|exact value\|stated fraction\|printed upper\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
} > "$out" 2>&1
fa=$(grep "fingerprint" "$tmp/r371b.a" | awk '{print $NF}'); fb=$(grep "fingerprint" "$tmp/r371b.b" | awk '{print $NF}')
va=$(grep "exact value" "$tmp/r371b.a" | awk '{print $NF}'); vb=$(grep "exact value" "$tmp/r371b.b" | awk '{print $NF}')
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|core error" "$out")
if [ -n "$fa" ] && [ "$fa" = "$fb" ] && [ -n "$va" ] && [ "$va" = "$vb" ] && [ "$pass" = 2 ] && [ "$ref" = 8 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs give the stated fraction, the single maximizing knot 116123 and the same fingerprint of all 655,359 knot values; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad; fingerprints $fa / $fb)" >> "$out"; tail -1 "$out"; exit 1
fi
