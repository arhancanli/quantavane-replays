#!/bin/sh
# Replay of the 73a lower-bound row [CS2019] (Flt(3) >= 2 + sqrt 2, Codenotti and Santos's hollow tetrahedron):
# hollowness and the lattice width exactly in Q(sqrt 2), with controls; no download.
# Usage: ./run_all.sh   (Python 3.9+, standard library only; about a second).
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
python3 -I cs2019_check.py --controls > results/cs2019.txt 2>&1; status=$?
grep -H "RESULT\|controls all" results/cs2019.txt
exit $status
