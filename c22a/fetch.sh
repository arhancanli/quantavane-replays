#!/bin/sh
# Download the paper of [V2026] (#356) at the commit its row cites, into a new directory, and check its SHA-256.
# c22a_check.py does not read it: the replay is written from the paper's statements (Sections 4 to 6); the file is
# pinned here so that the text we worked from can be compared. The repository's verify.py is not downloaded.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
curl -s -f -L -o v2026_paper.pdf https://raw.githubusercontent.com/ArseniyVarlamov/ropelength-crossing-number/6d11254ed91fee29f29dae8fa891840b46c3b4df/paper.pdf
shasum -a 256 -c "$here/SHA256SUMS"
