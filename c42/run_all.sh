#!/bin/sh
# Fetch the certificates (data only; files already present are kept and checked again), then run both checks with
# their red controls; output in results/.
set -e
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
./fetch.sh
mkdir -p results
python3 -I c42_check.py certs --controls --finite | tee results/c42_check_controls.txt
python3 -I kos_check.py certs --controls | tee results/c42_kos2026.txt
