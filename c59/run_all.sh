#!/bin/sh
# Fetch the certificate (data only), then run the check with its red controls; output in results/.
set -e
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
[ -d certs/gomyway ] || ./fetch.sh
python3 -I c59_check.py certs/gomyway/c59_cert.json --controls | tee results/c59_check_controls.txt
