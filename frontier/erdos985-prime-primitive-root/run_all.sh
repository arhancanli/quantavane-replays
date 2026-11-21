#!/bin/sh
# Erdos problem 985, finite part (Howe, Huffman, Kohlmeyer, Kruer @ 99fd0a7): two separately written exact programs.
# A = erdos985_a.py (exponent test for primitive roots, depth-first enumeration of the certificate's candidates, trial
# division); B = erdos985_b.c (order by walking powers for p <= 2791; segmented sieve for omega(n-1) and primality of the
# candidates; bit sieve for the extension). Both: (a) the 405 odd primes p <= 2791, (b) every prime modulus of the
# companion certificate, (c) every odd prime p <= 10^7; then two red controls each. No program of the source was read or run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+ and a C compiler; about one minute, one process.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
T="$root/erdos985/small-prime-witnesses.tsv"; W="$root/erdos985/witnesses"
out="$here/results/erdos985.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
${CC:-cc} -O2 -o "$tmp/b" "$here/erdos985_b.c" || { echo "RESULT: FAIL (compile)"; exit 1; }
python3 -I "$here/erdos985_a.py" "$T" "$W" --bound 10000000 > "$tmp/a" 2>&1
"$tmp/b" "$T" "$W" 10000000 > "$tmp/bo" 2>&1
for c in witness strip; do
  python3 -I "$here/erdos985_a.py" "$T" "$W" --bound 1000 --control $c > "$tmp/a.$c" 2>&1
  "$tmp/b" "$T" "$W" 1000 $c > "$tmp/b.$c" 2>&1
done
{
  echo "== program A (erdos985_a.py)"; cat "$tmp/a"
  echo "== program B (erdos985_b.c)"; cat "$tmp/bo"
  echo "== red control 1: table row p = 7 given q = 2 (order 3 mod 7)"
  grep -h "CONTROL\|FAIL\|VERDICT" "$tmp/a.witness" "$tmp/b.witness"
  echo "== red control 2: first prime modulus of the rank-3 strip given a prime that is not a primitive root"
  grep -h "CONTROL\|FAIL\|VERDICT" "$tmp/a.strip" "$tmp/b.strip"
} > "$out" 2>&1
sed -n 's/^A (/(/p' "$tmp/a" > "$tmp/a.n"; sed -n 's/^B (/(/p' "$tmp/bo" > "$tmp/b.n"
same=0; cmp -s "$tmp/a.n" "$tmp/b.n" && [ -s "$tmp/a.n" ] && same=1
rm -rf "$tmp"
pass=$(grep -c "^[AB]: VERDICT PASS" "$out"); ref=$(grep -c "^[AB]: VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback" "$out")
if [ "$same" = 1 ] && [ "$pass" = 2 ] && [ "$ref" = 4 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (A and B agree line by line: 405 odd primes p <= 2791 and all 18,802 prime moduli of the certificate have a prime primitive root < p, as does every odd prime p <= 10^7; every red control refused)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (A and B agree: $same, pass $pass, refused $ref)" >> "$out"; tail -1 "$out"; exit 1
