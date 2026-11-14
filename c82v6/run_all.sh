#!/bin/sh
# Replay PR #280's update of 14 November 2026 (C82 v6) with two separately written programs:
#  program 1, our existing ../c82/c82_check.py (Gershgorin discs on the Weierstrass matrix, per-root series on the
#  panels), through c82v6_first.py, which only registers the row; and program 2, c82v6_second.py (Newton inclusion with
#  exact evaluation, Horner values at the panel centres, Taylor polynomials in the panel variable), whose integral runs
#  in two parts at once, then its Jensen self-test, then the combination with the red controls.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only; about 15 minutes on two
# cores, under 100 MB; every process stays below 20 minutes of CPU.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
cert="$root/pr280f/certificate.json"
sha=29dcbf21a10ddb600a9866fe8ee7e7b43e4357926d74f374a8b6fbf41ada71b5
tmp=$(mktemp -d)
mkdir -p "$here/results"
python3 -I "$here/c82v6_first.py" "$root" --controls --workers=2 > "$here/results/pr280f_first.txt" 2>&1
python3 -I "$here/c82v6_second.py" selftest > "$tmp/selftest.txt" 2>&1
python3 -I "$here/c82v6_second.py" roots "$cert" "$tmp/approx.json" > "$tmp/roots.log" 2>&1
A="--approx $tmp/approx.json"
python3 -I "$here/c82v6_second.py" part "$cert" 0 2 "$tmp/pos0.json" pos $A > "$tmp/pos0.log" 2>&1 &
python3 -I "$here/c82v6_second.py" part "$cert" 1 2 "$tmp/pos1.json" pos $A > "$tmp/pos1.log" 2>&1
wait
python3 -I "$here/c82v6_second.py" part "$cert" 0 1 "$tmp/jen.json" jensen $A > "$tmp/jen.log" 2>&1
{
  cat "$tmp/selftest.txt"
  python3 -I "$here/c82v6_second.py" combine "$cert" 0.2519694883 $sha "$tmp/pos0.json" "$tmp/pos1.json" "$tmp/jen.json" $A
  echo "  parts (each certifies the cached approximations again):"
  cat "$tmp/roots.log" "$tmp/pos0.log" "$tmp/pos1.log" "$tmp/jen.log" | sed 's/^/    /'
} > "$here/results/pr280f_second.txt" 2>&1
rm -rf "$tmp"
p1=$(grep -c "RESULT: ALL CHECKS PASS" "$here/results/pr280f_first.txt")
p2=$(grep -c "RESULT: ALL CHECKS PASS" "$here/results/pr280f_second.txt")
b1=$(grep "^  B:" "$here/results/pr280f_first.txt" | sed 's/.*\[\([0-9.]*\), \([0-9.]*\)\].*/\1 \2/')
b2=$(grep "^  B:" "$here/results/pr280f_second.txt" | sed 's/.*\[\([0-9.]*\), \([0-9.]*\)\].*/\1 \2/')
over=$(python3 -I -c "import sys; a, b, c, d = map(__import__('fractions').Fraction, sys.argv[1:]); print(int(a <= d and c <= b))" $b1 $b2 2>/dev/null)
echo "program 1: B in [$b1]; program 2: B in [$b2]" > "$here/results/pr280f.txt"
if [ "$p1" = 1 ] && [ "$p2" = 2 ] && [ "$over" = 1 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs accept the certificate, B < 0.2519694883 in both, and their enclosures of B overlap)" >> "$here/results/pr280f.txt"
  tail -1 "$here/results/pr280f.txt"; exit 0
fi
echo "RESULT: FAIL (program 1 pass $p1, program 2 pass $p2, overlap $over)" >> "$here/results/pr280f.txt"; tail -1 "$here/results/pr280f.txt"; exit 1
