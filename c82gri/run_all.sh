#!/bin/sh
# Replay the 82a row [Gri26] C_82 <= 0.2536331090204145 with two separately written programs: program 1, our existing
# ../c82/c82_check.py, through c82gri_first.py, which only registers the row; and program 2, ../c82v6/c82v6_second.py
# (Newton inclusion with exact evaluation, Horner values at panel centres, Taylor polynomials in the panel variable),
# its integral in two parts, then its Jensen self-test and the combination with the red controls.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs; run ./fetch.sh first). Python 3.9+, standard library only;
# single-threaded (load permitting, the two parts of program 2 may run at once).
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
cert="$root/gri26/certificate.json"
sha=ccbda5c371802b765e7a01f5a9a5a7b7b84c907719a28ff79a7ef44ded0f3aee
P2="$here/../c82v6/c82v6_second.py"
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c82gri_first.py" "$root" --controls --workers=1 > "$here/results/gri26_first.txt" 2>&1
python3 -I "$P2" selftest > "$tmp/selftest.txt" 2>&1
python3 -I "$P2" roots "$cert" "$tmp/approx.json" > "$tmp/roots.log" 2>&1
A="--approx $tmp/approx.json"
python3 -I "$P2" part "$cert" 0 2 "$tmp/pos0.json" pos $A > "$tmp/pos0.log" 2>&1
python3 -I "$P2" part "$cert" 1 2 "$tmp/pos1.json" pos $A > "$tmp/pos1.log" 2>&1
python3 -I "$P2" part "$cert" 0 1 "$tmp/jen.json" jensen $A > "$tmp/jen.log" 2>&1
{
  cat "$tmp/selftest.txt"
  python3 -I "$P2" combine "$cert" 0.2536331090204145 $sha "$tmp/pos0.json" "$tmp/pos1.json" "$tmp/jen.json" $A
  echo "  parts (each certifies the cached approximations again):"
  cat "$tmp/roots.log" "$tmp/pos0.log" "$tmp/pos1.log" "$tmp/jen.log" | sed 's/^/    /'
} > "$here/results/gri26_second.txt" 2>&1
rm -rf "$tmp"
p1=$(grep -c "RESULT: ALL CHECKS PASS" "$here/results/gri26_first.txt")
p2=$(grep -c "RESULT: ALL CHECKS PASS" "$here/results/gri26_second.txt")
b1=$(grep "^  B:" "$here/results/gri26_first.txt" | sed 's/.*\[\([0-9.]*\), \([0-9.]*\)\].*/\1 \2/')
b2=$(grep "^  B:" "$here/results/gri26_second.txt" | sed 's/.*\[\([0-9.]*\), \([0-9.]*\)\].*/\1 \2/')
over=$(python3 -I -c "import sys; a, b, c, d = map(__import__('fractions').Fraction, sys.argv[1:]); print(int(a <= d and c <= b))" $b1 $b2 2>/dev/null)
echo "program 1: B in [$b1]; program 2: B in [$b2]" > "$here/results/gri26.txt"
if [ "$p1" = 1 ] && [ "$p2" = 2 ] && [ "$over" = 1 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs accept the certificate, B < 0.2536331090204145 in both, and their enclosures of B overlap)" >> "$here/results/gri26.txt"
  tail -1 "$here/results/gri26.txt"; exit 0
fi
echo "RESULT: FAIL (program 1 pass $p1, program 2 pass $p2, overlap $over)" >> "$here/results/gri26.txt"; tail -1 "$here/results/gri26.txt"; exit 1
