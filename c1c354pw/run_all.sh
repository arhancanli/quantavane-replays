#!/bin/sh
# Replay the piecewise-polynomial C_1c certificate of PR #354 (13 November 2026) with two separately written exact
# programs (A: symbolic integrals of pairs of pieces, Bernstein coefficients; B: the jump formula for g, centred
# Taylor forms), then three red controls that both must refuse, and one informative run of program A (the bound T in
# place of T - margin). The repository's verifiers and construction code were neither read nor run.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+, standard library only (python3 -I).
# About 6 minutes on two cores (each run about 45 s), under 200 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
f="$root/pr354pw/pwcert_p20_dps40_S120.json"
out="$here/results/pr354pw.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
T=13866842223167558038152245155365550010625/5444517870735015415413993718908291383296
pair() {  # pair NAME [ARGS...]: both programs at once, outputs kept apart
  name=$1; shift
  python3 -I "$here/c1c354pw_a.py" "$f" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c1c354pw_b.py" "$f" "$@" > "$tmp/$name.b" 2>&1
  wait
}
pair main
python3 -I "$here/c1c354pw_a.py" "$f" --claim "$T" > "$tmp/claimT.a" 2>&1
pair above --claim "13866842223167558038152245155365550010626/5444517870735015415413993718908291383296"
pair bump --bump 4 0 -1/100000000000000000000
pair late --bump 9 0 -1/1000000000000000000000000000000
{
  echo "== the certificate"; grep -v "run time" "$tmp/main.a"; grep -v "run time" "$tmp/main.b"
  echo "== information (program A only; centred forms cannot certify a bound that is attained): T in place of T - margin"
  grep "CONTROL\|every cell\|certified\|VERDICT" "$tmp/claimT.a"
  echo "== red control 1: the bound on g raised to T + 2^-132 (one unit in the numerator of T)"
  grep "CONTROL\|every cell\|VERDICT" "$tmp/above.a" "$tmp/above.b" | sed 's/^.*\/above\.[ab]://'
  echo "== red control 2: the constant coefficient of piece 4 lowered by 10^-20"
  grep "CONTROL\|int f = 1\|every cell\|digest\|VERDICT" "$tmp/bump.a" "$tmp/bump.b" | sed 's/^.*\/bump\.[ab]://'
  echo "== red control 3: h on [1 - t1, 1) lowered from T by 10^-30"
  grep "CONTROL\|last two\|h = T\|int f = 1\|every cell\|digest\|VERDICT" "$tmp/late.a" "$tmp/late.b" | sed 's/^.*\/late\.[ab]://'
} > "$out" 2>&1
rm -rf "$tmp"
da=$(grep "^A: digest" "$out" | head -1 | awk '{print $NF}'); db=$(grep "^B: digest" "$out" | head -1 | awk '{print $NF}')
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback" "$out")
if [ -n "$da" ] && [ "$da" = "$db" ] && [ "$pass" = 3 ] && [ "$ref" = 6 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs accept the certificate and give the same exact g, digest $da; three red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (pass $pass, refused $ref, digests $da $db)" >> "$out"; tail -1 "$out"; exit 1
