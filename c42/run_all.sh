#!/bin/sh
# Fetch the two certificates (data only), then run the check with its red controls and the finite-n look;
# output in results/.
set -e
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
[ -d certs/griego ] || ./fetch.sh
mkdir -p results
python3 -I c42_check.py certs --controls --finite | tee results/c42_check_controls.txt
