#!/bin/sh
# 48a: C_48 = c_0^2 from the system of equations of Davis and Power (arXiv 2604.03170v1, (3)-(6)), by two
# separately written programs (c48a_a.py: balls on exact rationals, Q by its alternating series, pi by Machin;
# c48a_b.py: decimal intervals, Q by its positive series, pi by Bailey-Borwein-Plouffe), then red controls that both
# must refuse. Nothing is downloaded and no program of any source was read or run.
# Usage: ./run_all.sh   Python 3.9+; about a minute, under 100 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c48a.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
pair() {  # pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c48a_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c48a_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair m --claim-round 5.33386
pair x1 --claim-round 5.33387
pair x2 --claim-round 5.33386 --variant "t0=sqrt(ln2)"
pair x3 --claim-round 5.33386 --variant "p0-without-2"
pair x4 --claim-round 5.33386 --variant "H-with-Q-once"
{
  echo "== 48a [DP2026]: program A (c48a_a.py), then program B (c48a_b.py)"
  cat "$tmp/m.a" "$tmp/m.b"
  for x in x1 x2 x3 x4; do
    case $x in
      x1) echo "== red control 1: c_0^2 claimed to round to 5.33387";;
      x2) echo "== red control 2: t_0 = sqrt(ln 2) in place of sqrt(2 ln 2)";;
      x3) echo "== red control 3: p_0 = e^(-a^2/2), the factor 2 dropped";;
      x4) echo "== red control 4: H(x) = 2x e^(-x^2/2) + Q(x), the factor 2 on Q dropped";;
    esac
    grep -h "CONTROL\|C_48\|rounds to 5\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
  echo "== agreement: every quantity printed by both programs (enclosures must overlap)"
  python3 -I "$here/compare.py" "$tmp/m.a" "$tmp/m.b"; echo "pair exit $?"
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
agree=$(grep -c "^pair exit 0" "$out")
if [ "$pass" = 2 ] && [ "$ref" = 8 ] && [ "$bad" = 0 ] && [ "$agree" = 1 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs give C_48 = c_0^2 = 5.33386354161609721600557861141..., overlapping enclosures for every quantity; four red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agree=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
