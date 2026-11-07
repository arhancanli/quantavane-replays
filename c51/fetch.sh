#!/bin/sh
# Download the one data file replayed here (no code) into a new directory and check its SHA-256: the exact parameters
# of [R2026b], certificates/maximum_term.json of AndreRoehrig/optimization-constants-certificates at the commit PR #184
# and the 51a row pin (ed6387399cf6; v1.0.2). Its scripts (certify_maximum_term.py, audit.py) are neither downloaded
# nor run. A file already present is kept and checked again.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ ! -f r2026b/maximum_term.json ]; then
  if [ -d r2026b ] && [ -n "$(ls -A r2026b)" ]; then echo "refusing: $root/r2026b is not empty"; exit 1; fi
  mkdir -p r2026b
  curl -s -f -L -o r2026b/maximum_term.json \
    "https://raw.githubusercontent.com/AndreRoehrig/optimization-constants-certificates/ed6387399cf6de02eaa4cf10c43787633496619a/certificates/maximum_term.json"
fi
shasum -a 256 -c "$here/SHA256SUMS"
