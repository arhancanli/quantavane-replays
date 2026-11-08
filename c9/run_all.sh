#!/bin/sh
# Replay the 9a rows 134753^(1/10) ~ 3.258020 ([IRCR2026], added in PR #174: the set is independent in the strong
# power C_7^10, two routes, strong_check.c, and the printed 3.258020 is at most 134753^(1/10), exact integers) and
# M_40^(1/200) = 3.258789153908... ([Gao2026]: the base gadget, the k = 2 product code checked by strong_check.c, the
# recursion to M_40, gao_check.py); red controls. One output file per row in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). Needs a C compiler and Python 3.9+; about two minutes.
# Peak memory about 40 MB.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
root=${1:-certs}
mkdir -p results bin
cc -O2 -Wall -o bin/strong_check strong_check.c || exit 1
F="$root/ircr2026/CC_10_7_134753.txt"
out=results/ircr2026.txt
status=0
{
  echo "== the set"
  ./bin/strong_check 7 10 "$F" || status=1
  echo "== the printed value"
  python3 -I -c "
n, lo = 134753, 3258020
print('3.258020^10 <= 134753 exactly:', lo ** 10 <= n * 10 ** 60, '; 3.258021^10 <= 134753:', (lo + 1) ** 10 <= n * 10 ** 60)
" || status=1
  echo "== red controls (each must be refused)"
  ./bin/strong_check 7 10 "$F" --inject-neighbour && status=1
  ./bin/strong_check 7 10 "$F" --inject-duplicate && status=1
  ./bin/strong_check 7 10 "$F" --wrong-graph && status=1
} > $out 2>&1
grep -H "RESULT\|exactly" $out
out=results/gao2026.txt
{
  python3 -I gao_check.py "$root" "$root/gao2026/code_k2.txt" --controls || status=1
  echo "== the k = 2 product code (Lemma 4 on two copies of the base gadget)"
  ./bin/strong_check 7 10 "$root/gao2026/code_k2.txt" || status=1
} > $out 2>&1
grep -H "RESULT" $out
exit $status
