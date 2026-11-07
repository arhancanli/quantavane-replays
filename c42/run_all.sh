#!/bin/sh
# Fetch the certificates (data only; files already present are kept and checked again), then run the checks with
# their red controls; output in results/. thmA_check.py writes its own files (thmA_*.txt): the three Theorem A rows,
# the self-tests (with the convergence table) and the red controls; it takes about 25 minutes on a loaded machine.
set -e
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
./fetch.sh
mkdir -p results
python3 -I c42_check.py certs --controls --finite | tee results/c42_check_controls.txt
python3 -I kos_check.py certs --controls | tee results/c42_kos2026.txt
python3 -I thmA_check.py certs --selftest --controls --write
