#!/bin/sh
# Replay of the 89c candidate lower bound of PR #308 from its certificate and graph list (fetched by fetch.sh): the
# certificate's structure, the completeness of the graph list, and the coefficient inequality for every 8-vertex graph,
# exactly, with red controls and a pure-Python recomputation for three graphs. One output file in results/.
# Usage: ./run_all.sh [CERTS_ROOT]   (default: ./certs). A C compiler and Python 3.9+; about a minute, under 400 MB.
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$here/results" "$here/bin"
cc -O2 -o "$here/bin/flag_count" "$here/flag_count.c" || exit 1
python3 -I "$here/flag_check.py" "$root" "$here/bin" --controls > "$here/results/pin2026.txt" 2>&1; status=$?
grep -H "RESULT\|controls all" "$here/results/pin2026.txt"
exit $status
