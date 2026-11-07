#!/bin/sh
# Download the C89 (Ramsey multiplicity of K4) step graphons replayed here (data only, no code), each into a new
# directory, at the commits the 89a rows and PR #270 pin; then check their SHA-256 against SHA256SUMS.
#   vsj   [VSJ2026] lukevs/k4-ramsey @ 2cf216b, the final 3840-class witness
#   sk    [SK2026]  JeetuSK0808/tao-constants-certificates @ 3d9431a
#   kos   [Kos2026] DennisKoshta/c89-k4-ramsey-multiplicity-certificate @ 4fdbe0c: its amplitudes only; the 7680-class
#                   witness is rebuilt from them and VSJ's table by k4_check.py, and must match the stated dense hash
#   vid   #270's [Pat2026] row, vidigoat/optimization-records @ fb80010 (the commit the row links), second round
#   vid1  the same commit, the first round (the value in the row's parenthesis)
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # get DIR FILE REPO COMMIT PATH: one file into a directory that must be new or empty
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "https://raw.githubusercontent.com/$3/$4/$5"
}
get vsj final3840.json.gz lukevs/k4-ramsey 2cf216ba885fbb1b2d1b2d4920efeaa6ba438842 data/constructions/final3840.json.gz
get sk c89_stepgraphon_D20.json.gz JeetuSK0808/tao-constants-certificates 3d9431a00abf03be3752f8f9f126dec81ac34987 \
    c89-k4-ramsey-multiplicity/c89_stepgraphon_D20.json.gz
get kos c89_sign3_crude_0.03013888570.delta.json.gz DennisKoshta/c89-k4-ramsey-multiplicity-certificate \
    4fdbe0cecbb2519a17b90468a592c789a818f744 delta/c89_sign3_crude_0.03013888570.delta.json.gz
get vid c89_stepgraphon_D30_round2.json.gz vidigoat/optimization-records fb800103581d5c63fa096805f73eff48fce32c98 \
    tao/certificates/c89/c89_stepgraphon_D30_round2.json.gz
get vid1 c89_stepgraphon_D30.json.gz vidigoat/optimization-records fb800103581d5c63fa096805f73eff48fce32c98 \
    tao/certificates/c89/c89_stepgraphon_D30.json.gz
cd "$root" && shasum -a 256 -c "$here/SHA256SUMS"
