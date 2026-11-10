#!/bin/sh
# Download the four step-function certificates of 9 November 2026 replayed here, at the commits their pull requests
# pin, into a new directory, and check their SHA-256 against SHA256SUMS:
#   #354 (SeanSagait): K262144_0.4103768.txt, the 1c lower bound 0.4103768 (K = 262144; the first number of the file is
#        K, then 262400 non-negative integer weights); copied into JSON by a few lines of Python below (integers only,
#        each token checked to be a canonical decimal integer; nothing downloaded is run);
#   #270 (vidigoat, updates of 9 Nov): c1c_K98304_v4.json (0.4103604995), c1c_K98304_v5.json (0.4103608339) and
#        c1a_n262144_v2.json (the 1a upper bound 1.5027427356); update of 10 Nov: c1d_n163840.json (the 1d upper bound
#        1.4546733458) and c1a_n262144_v3.json (1.5027427341).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
get() { curl -s -f -L -o "$1" "$2"; }
get sean_K262144_0.4103768.txt https://raw.githubusercontent.com/SeanSagait/c1c-autocorrelation-certificate/9e7f3d19b93a463ba6073f69d34b65d915c49f25/K262144_0.4103768.txt
V=https://raw.githubusercontent.com/vidigoat/optimization-records
get vidigoat_c1c_K98304_v4.json $V/358b200b9810ba0b4041497dd1bbd5a31d0952e5/tao/certificates/c1c_v4/c1c_K98304_v4.json
get vidigoat_c1c_K98304_v5.json $V/1b39a28d3cfe34ac73cf30a08b5afae9c2cdf4b2/tao/certificates/c1c_v5/c1c_K98304_v5.json
get vidigoat_c1a_n262144_v2.json $V/1b39a28d3cfe34ac73cf30a08b5afae9c2cdf4b2/tao/certificates/c1a_v2/c1a_n262144_v2.json
# #270, update of 10 Nov (eb710a8): the new 1d row 1.4546733458 and the 1a row updated in place to 1.5027427341
get vidigoat_c1d_n163840.json $V/a559680e416be0eae03955c58372efff4a551fbb/tao/certificates/c1d/c1d_n163840.json
get vidigoat_c1a_n262144_v3.json $V/a559680e416be0eae03955c58372efff4a551fbb/tao/certificates/c1a_v3/c1a_n262144_v3.json
python3 -I -c "
import json, re, sys
tok = open('sean_K262144_0.4103768.txt').read().split()
if not all(re.fullmatch(r'0|[1-9][0-9]*', t) for t in tok): sys.exit('a token is not a canonical non-negative integer')
K, w = int(tok[0]), [int(t) for t in tok[1:]]
if (K, len(w)) != (262144, 262400): sys.exit('unexpected K or weight count')
json.dump({'K': K, 'weights_int': w}, open('sean_c1c_K262144.json', 'w'))
"
shasum -a 256 -c "$here/SHA256SUMS"
