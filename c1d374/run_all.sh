#!/bin/sh
# Replay the 1d step function of PR #374 with two separately written exact programs (A: Python's float read as an
# exact integer, one Kronecker squaring with balanced digits in Python integers; B: its own JSON scanner with the C
# library's strtod, plain loops over every knot in C, signed 128-bit integers, after a Cauchy-Schwarz overflow guard),
# then red controls that both must refuse. No program of any source was read or run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs)
# Python 3.9+ and a C compiler; under half a minute, under 300 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/c1d374.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
F="$root/pr374/c1d_n102400_1.4513064007.json"
C=513655524515252183593642627316297728000/353926313748615382700541734286037793169
pair() {  # pair NAME [ARGS...]: both programs at once, outputs kept apart
  name=$1; shift
  python3 -I "$here/c1d374_a.py" "$F" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c1d374_b.py" "$F" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair r374 --claim $C --argmax 30289 --printed 1.4513064006879361
pair x1 --bump 50000 1 --claim $C
pair x2 --printed 1.4513064006879360
pair x3 --bump 0 -83573763203730464 --claim $C
pair x4 --bump 102399 1000000 --claim $C
{
  echo "== PR #374 (1d): c1d_n102400_1.4513064007.json (cweill/c1d-upper-bound @ 7889e93), program A then program B"
  echo "   (checked against the stated fraction, the stated knot k = 30289 and the printed 1.4513064006879361)"
  grep -h -v "run time" "$tmp/r374.a" "$tmp/r374.b"
  for x in x1 x2 x3 x4; do
    case $x in
      x1) echo "== red control 1: height 50000 plus one unit, against the stated fraction";;
      x2) echo "== red control 2: claimed as an upper bound 1.4513064006879360";;
      x3) echo "== red control 3: height 0 (41786881601865232) negated, against the stated fraction";;
      x4) echo "== red control 4: the last height plus 10^6 units (S changes), against the stated fraction";;
    esac
    grep -h "CONTROL\|control\|exact value\|stated fraction\|printed upper\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
} > "$out" 2>&1
fa=$(grep "fingerprint" "$tmp/r374.a" | awk '{print $NF}'); fb=$(grep "fingerprint" "$tmp/r374.b" | awk '{print $NF}')
va=$(grep "exact value" "$tmp/r374.a" | awk '{print $NF}'); vb=$(grep "exact value" "$tmp/r374.b" | awk '{print $NF}')
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|core error" "$out")
if [ -n "$fa" ] && [ "$fa" = "$fb" ] && [ -n "$va" ] && [ "$va" = "$vb" ] && [ "$pass" = 2 ] && [ "$ref" = 8 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs give the stated fraction, the single maximizing knot 30289 and the same fingerprint of all 204,799 knot values; the printed 1.4513064006879361 holds, 1.2e-18 above the exact value; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad; fingerprints $fa / $fb)" >> "$out"; tail -1 "$out"; exit 1
fi
