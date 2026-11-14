#!/bin/sh
# Download the piecewise-polynomial C_1c certificate of PR #354 (data only, no code) into a new directory and check
# its SHA-256 (the one stated in the pull request and in the repository's SHA256SUMS).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d pr354pw ] && [ -n "$(ls -A pr354pw)" ]; then echo "refusing: $root/pr354pw is not empty"; exit 1; fi
mkdir -p pr354pw
# SeanSagait/c1c-autocorrelation-certificate @ 308005a (13 November 2026)
curl -s -f -L -o pr354pw/pwcert_p20_dps40_S120.json \
  https://raw.githubusercontent.com/SeanSagait/c1c-autocorrelation-certificate/308005a40a6df640b4f2124ef7d812522b792f42/pwcert_p20_dps40_S120.json
(cd pr354pw && shasum -a 256 -c "$here/SHA256SUMS")
