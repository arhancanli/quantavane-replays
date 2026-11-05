#!/bin/sh
# Replay the C89c (Ramsey multiplicity of K_6) upper bounds of PRs #279 and #331 and the C89b (K_5) bounds of PRs #307
# and #277 exactly, and the PPSS 2024 value as a control; one output file per row in results/, red controls included,
# after the self-tests (brute force, f = 1/2, Franek 2002).
# Usage: ./run_all.sh [CERTS_ROOT] [--workers=N]   (default: ./certs). Python 3.9+, standard library only, and cc.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
shift 2>/dev/null
mkdir -p "$here/results"
status=0
python3 -I "$here/c89_check.py" selftest "$@" > "$here/results/selftest.txt" 2>&1 || status=1
python3 -I "$here/c89b_check.py" selftest "$@" > "$here/results/selftest_b.txt" 2>&1 || status=1
grep -H "RESULT" "$here/results/selftest.txt" "$here/results/selftest_b.txt"
for row in sk2026 sciel2026; do
  python3 -I "$here/c89_check.py" $row "$root" --controls "$@" > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
for row in pr307 pr277; do
  python3 -I "$here/c89b_check.py" $row "$root" --controls "$@" > "$here/results/$row.txt" 2>&1 || status=1
  grep -H "RESULT\|run time" "$here/results/$row.txt" | tr '\n' ' '; echo
done
python3 -I "$here/c89b_check.py" ppss "$root" "$@" > "$here/results/ppss.txt" 2>&1 || status=1
grep -H "RESULT\|run time" "$here/results/ppss.txt" | tr '\n' ' '; echo
exit $status
