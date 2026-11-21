#!/bin/sh
# Download the codes of P. R. A. Iskandar, arXiv 2611.24958v1 (ancillary data files only, not the three .py verifiers)
# into a new directory and check the SHA-256: CLAIMS.csv, ORBIT_MANIFEST.csv, codes/A{n}_{d}_{w}_{M}.txt (13 full word
# lists) and reps/..._reps.txt (13 orbit-representative lists). One request every 15 s (about 7 minutes).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
base=${ANC_BASE:-https://arxiv.org/src/2611.24958v1/anc}
mkdir -p "$root" && cd "$root"
if [ -d anc ] && [ -n "$(ls -A anc)" ]; then echo "refusing: $root/anc is not empty"; exit 1; fi
mkdir -p anc/codes anc/reps
first=1
for f in CLAIMS.csv ORBIT_MANIFEST.csv \
  A16_4_7_841 A20_4_5_876 A25_4_6_7800 A29_4_5_4491 A29_6_5_333 A36_4_5_11207 A38_4_5_13737 A39_4_5_15618 \
  A44_4_5_25585 A54_4_5_60632 A55_4_5_66144 A60_4_5_94105 A62_4_5_107360; do
  case $f in *.csv) list=$f;; *) list="codes/$f.txt reps/${f}_reps.txt";; esac
  for g in $list; do
    [ $first = 1 ] || sleep 15; first=0
    curl -s -f -L -A "quantavane-replays fetch.sh" -o "anc/$g" "$base/$g"
  done
done
shasum -a 256 -c "$here/SHA256SUMS"
