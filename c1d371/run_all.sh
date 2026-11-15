#!/bin/sh
# Replay the 1d step function of PR #371 with two separately written exact programs (A: one Kronecker squaring with
# balanced digits in Python integers; B: its own reader and plain loops over every knot in C, signed 128-bit
# integers), then red controls that both must refuse. No program of any source was read or run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs)
# Python 3.9+ and a C compiler; about four minutes on two cores, under 200 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/c1d371.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
F="$root/pr371/C1d_327680_1.4546733294.txt"
C=8076601902879958587898310656/5552175694520934184034268245
pair() {  # pair NAME [ARGS...]: both programs at once, outputs kept apart
  name=$1; shift
  python3 -I "$here/c1d371_a.py" "$F" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c1d371_b.py" "$F" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair r371 --claim $C --argmax 99803 --printed 1.4546733294
pair x1 --bump 50000 1 --claim $C
pair x2 --printed 1.4546733293
pair x3 --bump 0 -1811939703810 --claim $C
pair x4 --bump 327679 1000000 --claim $C
{
  echo "== PR #371 (1d): C1d_327680_1.4546733294.txt (intervalplace/frontier @ 316b76c), program A then program B"
  grep -h -v "run time" "$tmp/r371.a" "$tmp/r371.b"
  for x in x1 x2 x3 x4; do
    case $x in
      x1) echo "== red control 1: height 50000 plus one unit, against the stated fraction";;
      x2) echo "== red control 2: claimed as an upper bound 1.4546733293 (below the exact value)";;
      x3) echo "== red control 3: height 0 (905969851905) negated, against the stated fraction";;
      x4) echo "== red control 4: the last height plus 10^6 units (S changes), against the stated fraction";;
    esac
    grep -h "CONTROL\|exact value\|stated fraction\|printed upper\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
} > "$out" 2>&1
fa=$(grep "fingerprint" "$tmp/r371.a" | awk '{print $NF}'); fb=$(grep "fingerprint" "$tmp/r371.b" | awk '{print $NF}')
va=$(grep "exact value" "$tmp/r371.a" | awk '{print $NF}'); vb=$(grep "exact value" "$tmp/r371.b" | awk '{print $NF}')
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|core error" "$out")
if [ -n "$fa" ] && [ "$fa" = "$fb" ] && [ -n "$va" ] && [ "$va" = "$vb" ] && [ "$pass" = 2 ] && [ "$ref" = 8 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs give the stated fraction, the single maximizing knot 99803 and the same fingerprint of all 655,359 knot values; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad; fingerprints $fa / $fb)" >> "$out"; tail -1 "$out"; exit 1
fi
