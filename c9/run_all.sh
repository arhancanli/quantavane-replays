#!/bin/sh
# Replay the 9a row 134753^(1/10) ~ 3.258020 ([IRCR2026], added in PR #174): the set is independent in the strong
# power C_7^10 (two routes, strong_check.c), and the printed 3.258020 is at most 134753^(1/10) (exact integers);
# red controls. One output file in results/. Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs).
# Needs a C compiler and Python 3.9+; about a minute (34 s for the main check). Peak memory about 40 MB.
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
exit $status
