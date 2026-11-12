#!/bin/sh
# Download the 22b sources at the commit of [Claude2026] (PR #267), each into a new directory, check their SHA-256,
# and take the seven q coefficients from code/capaxis.py with take_q.py (ast only; nothing executed). main.tex and
# numbers.tex are read as text, design.json as data; c22b_check.py compares design.json and the q strings with its
# own copies (--sources). No other file of the repository is downloaded.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
base=https://raw.githubusercontent.com/YouhuaLi/tight-alternating-knot-constant/debf7c92f93bae64486e6302e565a2d3b1290f9c
mkdir paper code
curl -s -f -L -o paper/main.tex "$base/paper/main.tex"
curl -s -f -L -o paper/numbers.tex "$base/paper/numbers.tex"
curl -s -f -L -o design.json "$base/design.json"
curl -s -f -L -o code/capaxis.py "$base/code/capaxis.py"
shasum -a 256 -c "$here/SHA256SUMS"
python3 -I "$here/take_q.py" code/capaxis.py q_from_capaxis.json > /dev/null
echo "q strings written to $root/q_from_capaxis.json"
