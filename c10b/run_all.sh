#!/bin/sh
# 10b, the complex Grothendieck constant. No data is downloaded: every input is a definition or an exact decimal typed
# from the cited papers (see RESULT.md). No program of any source was read or run.
#  (1) [H1987] Haagerup's upper bound 8/(pi(K_0 + 1)) from its defining equation, by two separately written programs
#      (c10b_series.py: hypergeometric series in exact rationals; c10b_agm.py: elliptic integrals by the AGM in
#      outward-rounded decimal intervals); also, as information, Haagerup's conjectured value and Pisier's e^(1-gamma).
#  (2) [GFL26] (PR #196) K_G^C > 1.35584631827168: the finite conditions of the paper's Theorem 3.1 for its exact
#      parameters, by two separately written programs (c10b_gfl_a.py: balls on exact rationals, partial fractions,
#      E1 by its power series; c10b_gfl_b.py: decimal intervals, Laguerre forms, a recurrence over the cubic, E1 by
#      its continued fraction).
# Then red controls that both programs of each pair must refuse.
# Usage: ./run_all.sh   Python 3.9+; about three minutes on two cores, under 200 MB.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c10b.txt"
tmp=$(mktemp -d)
mkdir -p "$here/results"
h_pair() {  # h_pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c10b_series.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c10b_agm.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
g_pair() {  # g_pair NAME [ARGS...]
  name=$1; shift
  python3 -I "$here/c10b_gfl_a.py" "$@" > "$tmp/$name.a" 2>&1 &
  python3 -I "$here/c10b_gfl_b.py" "$@" > "$tmp/$name.b" 2>&1
  wait
}
h_pair h --claim-h 1.40491 --claim-p 1.52621
h_pair hx1 --claim-h 1.4049091327357955
h_pair hx2 --claim-h 1.40490
h_pair hx3 --claim-p 1.52620
g_pair g --claim-ratio 1.35584631827168 --claim-mu 0.716627735168304
g_pair gx1 --set rho 0.012
g_pair gx2 --set beta10 0.0148
g_pair gx3 --set zeta2 1.3819
g_pair gx4 --claim-ratio 1.35584631827169
g_pair gx5 --claim-mu 0.716627735168303
{
  echo "== (1) [H1987]: program 1 (c10b_series.py), then program 2 (c10b_agm.py)"
  cat "$tmp/h.a" "$tmp/h.b"
  echo "== (2) [GFL26], PR #196: program A (c10b_gfl_a.py), then program B (c10b_gfl_b.py)"
  grep -h -v "^ *[0-9.]* real" "$tmp/g.a" "$tmp/g.b"
  for x in hx1 hx2 hx3 gx1 gx2 gx3 gx4 gx5; do
    case $x in
      hx1) echo "== red control H1: [H1987] claimed as 1.4049091327357955 (the value cut after 16 decimals, below it)";;
      hx2) echo "== red control H2: [H1987] claimed as 1.40490";;
      hx3) echo "== red control H3: [P1978] claimed as 1.52620";;
      gx1) echo "== red control G1: rho = 0.012 (lowered by 3.9e-5): the Gram blocks still hold, S > 1";;
      gx2) echo "== red control G2: beta_10 = 0.0148 (raised by 0.7%): a Gram block fails, S < 1 still";;
      gx3) echo "== red control G3: zeta_2 = 1.3819 (the weight lowered): both conditions fail";;
      gx4) echo "== red control G4: the ratio claimed above 1.35584631827169";;
      gx5) echo "== red control G5: mu_W claimed below 0.716627735168303";;
    esac
    grep -h "CONTROL\|Gram blocks\|1 - S\|ratio (\|mu_W in\|printed upper\|claimed\|conditions\|VERDICT" "$tmp/$x.a" "$tmp/$x.b"
  done
  echo "== agreement: every quantity printed by both programs of a pair (enclosures must overlap)"
  python3 -I "$here/compare.py" "$tmp/h.a" "$tmp/h.b"; echo "pair H exit $?"
  python3 -I "$here/compare.py" "$tmp/g.a" "$tmp/g.b"; echo "pair G exit $?"
} > "$out" 2>&1
rm -rf "$tmp"
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
# the two programs of each pair must print the same enclosures to the shown digits (checked by hand in RESULT.md);
# here: 4 PASS (two pairs), 16 refusals (eight controls, two programs each), no error
agree=$(grep -c "^pair [HG] exit 0" "$out")
if [ "$pass" = 4 ] && [ "$ref" = 16 ] && [ "$bad" = 0 ] && [ "$agree" = 2 ]; then
  echo "RESULT: ALL CHECKS PASS ([H1987] and [GFL26] by two programs each; eight red controls refused by both programs of their pair)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad agreeing pairs=$agree)" >> "$out"; tail -1 "$out"; exit 1
fi
