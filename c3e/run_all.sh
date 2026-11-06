#!/bin/sh
# Fetch the two certificate files (data only), then run the check with its red controls; output in results/.
set -e
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
[ -d certs/mizopl ] || ./fetch.sh
python3 -I c3e_check.py certs --controls | tee results/c3e_check_controls.txt
