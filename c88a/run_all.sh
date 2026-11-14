#!/bin/sh
# 88a: check the admissible tuple behind every upper-bound row with two separately written programs (A: occupied
# residue classes; B: a Chinese-remainder witness checked by gcds, and for [Zha14] a direct strike of every prime
# p <= k0), then the red controls, which both programs must refuse. The sieve theorems DHL[k, 2] are not checked.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Python 3.9+; program B's zhang mode needs numpy.
# About two minutes; under 1 GB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
out="$here/results/tuples.txt"
mkdir -p "$here/results"
A="python3 -I $here/c88a_tuples.py"; B="python3 -I $here/c88a_crt.py"; Bz="python3 $here/c88a_crt.py"
row() {  # row LABEL KIND FILE K D
  echo "== $1"; $A "$2" "$root/$3" "$4" "$5"; $B "$2" "$root/$3" "$4" "$5"
}
ctl() {  # ctl LABEL KIND FILE K D OLD NEW
  echo "== red control, $1: entry $6 replaced by $7"
  $A "$2" "$root/$3" "$4" "$5" --swap "$6" "$7" | tail -3; $B "$2" "$root/$3" "$4" "$5" --swap "$6" "$7" | tail -3
}
{
  echo "== [Zha14] 70,000,000: the k0 = 3,500,000 primes after k0"
  $A zhang 3500000 70000000; $Bz zhang 3500000 70000000
  row "[Pol14a] 4680: k0 = 632 (Sutherland's tables)" list sutherland/admissible_632_4680.txt 632 4680
  row "[May15] 600: k = 105 (Sutherland's tables)" list sutherland/admissible_105_600.txt 105 600
  row "[Pol14b] 246: k = 50 (Sutherland's tables)" list sutherland/admissible_50_246.txt 50 246
  row "[Zha26] 240: k = 49 (the tuple in the author's certificate.json)" json primebound/certificate.json 49 240
  row "[Sta26] 240: k = 49 (Sutherland's tables)" list sutherland/admissible_49_240.txt 49 240
  row "[OAI26] 186: k = 40 (the tuple of the Lean statement)" lean primegaps186/Challenge.lean 40 186
  row "a second 40-tuple of diameter 186 (Sutherland's tables; information)" list sutherland/admissible_40_186.txt 40 186
  echo "== red control, [Zha14]: the FIRST 3,500,000 primes instead"
  $A zhang 3500000 70000000 --from-first | tail -3; $Bz zhang 3500000 70000000 --from-first | tail -3
  ctl "[Pol14a] (odd entry)" list sutherland/admissible_632_4680.txt 632 4680 12 13
  ctl "[Pol14a] (even entry, fills p = 11)" list sutherland/admissible_632_4680.txt 632 4680 12 2
  ctl "[May15] (odd entry)" list sutherland/admissible_105_600.txt 105 600 10 11
  ctl "[May15] (even entry, fills p = 3)" list sutherland/admissible_105_600.txt 105 600 10 2
  ctl "[Pol14b] (odd entry)" list sutherland/admissible_50_246.txt 50 246 4 5
  ctl "[Pol14b] (even entry, fills p = 3)" list sutherland/admissible_50_246.txt 50 246 4 2
  ctl "[Zha26] (odd entry)" json primebound/certificate.json 49 240 2 3
  ctl "[Zha26] (even entry, fills p = 3)" json primebound/certificate.json 49 240 2 4
  ctl "[Sta26] (odd entry)" list sutherland/admissible_49_240.txt 49 240 4 5
  ctl "[Sta26] (even entry, fills p = 3)" list sutherland/admissible_49_240.txt 49 240 4 2
  ctl "[OAI26] (odd entry)" lean primegaps186/Challenge.lean 40 186 2 3
  ctl "[OAI26] (even entry, fills p = 3)" lean primegaps186/Challenge.lean 40 186 2 4
  echo "== red control, [Pol14b]: the same tuple offered as a 51-tuple, and with diameter 244"
  $A list "$root/sutherland/admissible_50_246.txt" 51 246 | tail -1; $B list "$root/sutherland/admissible_50_246.txt" 51 246 | tail -1
  $A list "$root/sutherland/admissible_50_246.txt" 50 244 | tail -1; $B list "$root/sutherland/admissible_50_246.txt" 50 244 | tail -1
} > "$out" 2>&1
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); err=$(grep -c -i "traceback\|error" "$out")
if [ "$pass" = 16 ] && [ "$ref" = 30 ] && [ "$err" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (8 tuples accepted by both programs; 15 red controls refused by both)" >> "$out"; tail -1 "$out"; exit 0
fi
echo "RESULT: FAIL (pass $pass, refused $ref, errors $err)" >> "$out"; tail -1 "$out"; exit 1
