#!/bin/sh
# Replay Martinet's 87a row (Example 4.2 of Invent. Math. 44 (1978), p. 69; PR #332's 4 sqrt(70035)) with red
# controls; one output file in results/. No download: the inputs are the five primes of the example.
# Usage: ./run_all.sh. Python 3.9+, standard library only; under a second.
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/results"
python3 -I "$here/martinet_check.py" --controls > "$here/results/mar1978.txt" 2>&1; status=$?
grep -H "RESULT\|r = 13\|printed\|unramified at every" "$here/results/mar1978.txt" | grep -v "refused"
exit $status
