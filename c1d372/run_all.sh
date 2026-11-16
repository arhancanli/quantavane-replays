#!/bin/sh
# Replay the 1d step function of PR #372 with two separately written exact programs (A: one Kronecker squaring with
# balanced digits in Python integers; B: its own JSON scanner and plain loops over every knot in C, signed 128-bit
# integers, after a Cauchy-Schwarz overflow guard), then red controls that both must refuse. No program of any
# source was read or run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs)
# Python 3.9+ and a C compiler; under ten seconds, under 200 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/c1d372.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
F="$root/pr372/c006_31bbf628_int56.json"
C=21267647941510998263770007610303168512/14626277862557668614837367632277156521
pair() {  # pair NAME [ARGS...]: both programs at once, outputs kept apart
  name=$1; shift
  python3 -I "$here/c1d372_a.py" "$F" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c1d372_b.py" "$F" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair r372 --claim $C --argmax 12511 --printed 1.4540711000681049
pair p372 --printed 1.4540711000681048
pair x1 --bump 25000 1 --claim $C
pair x2 --printed 1.4540711000681047
pair x3 --bump 0 -35225177806011560 --claim $C
pair x4 --bump 51199 1000000 --claim $C
{
  echo "== PR #372 (1d): c006_31bbf628_int56.json (ArthurAnache/Hackathon-AI-Reasearch @ c0f8d18), program A then program B"
  echo "   (checked against the stated fraction, the stated knot k = 12511, and 1.4540711000681049, the printed figure rounded up)"
  grep -h -v "run time" "$tmp/r372.a" "$tmp/r372.b"
  echo "== the printed figure 1.4540711000681048 itself, as an upper bound (reader's note: it is below the exact value)"
  grep -h "exact value\|printed upper\|VERDICT" "$tmp/p372.a" "$tmp/p372.b"
  for x in x1 x2 x3 x4; do
    case $x in
      x1) echo "== red control 1: height 25000 plus one unit, against the stated fraction";;
      x2) echo "== red control 2: claimed as an upper bound 1.4540711000681047";;
      x3) echo "== red control 3: height 0 (17612588903005780) negated, against the stated fraction";;
      x4) echo "== red control 4: the last height plus 10^6 units (S changes), against the stated fraction";;
    esac
    grep -h "CONTROL\|exact value\|stated fraction\|printed upper\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
} > "$out" 2>&1
fa=$(grep "fingerprint" "$tmp/r372.a" | awk '{print $NF}'); fb=$(grep "fingerprint" "$tmp/r372.b" | awk '{print $NF}')
va=$(grep "exact value" "$tmp/r372.a" | awk '{print $NF}'); vb=$(grep "exact value" "$tmp/r372.b" | awk '{print $NF}')
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|core error" "$out")
if [ -n "$fa" ] && [ "$fa" = "$fb" ] && [ -n "$va" ] && [ "$va" = "$vb" ] && [ "$pass" = 2 ] && [ "$ref" = 10 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs give the stated fraction, the single maximizing knot 12511 and the same fingerprint of all 102,399 knot values; the printed 1.4540711000681048 is 3.5e-17 below the exact value, 1.4540711000681049 holds; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad; fingerprints $fa / $fb)" >> "$out"; tail -1 "$out"; exit 1
fi
