#!/bin/sh
# Round taoreplays22 (22 Nov 2026): 22 closed-form bounds from constant pages of teorth/optimizationproblems (main 87b4ebb),
# each evaluated by two separately written programs (decimal / integer fixed point with AGM logarithms) inside tao22_check.py,
# compared with the printed digits, with a red control per row (the printed value moved by two units must fail).
# Usage: ./run_all.sh   Python 3.9+; under a second.
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/results"
python3 -I "$here/tao22_check.py" | tee "$here/results/tao22.txt"
