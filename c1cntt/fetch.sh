#!/bin/sh
# Download the two latest certificates of PR #354 (SeanSagait/c1c-autocorrelation-certificate: K = 2^22 at commit
# 867e0e6, K = 2^24 at commit f2b97a3, the second xz-compressed and unpacked here with Python's lzma module), each
# into its own new directory, and check their SHA-256. Only these data files are fetched; the repository's programs
# are not.  Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in k22 k24; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
base=https://github.com/SeanSagait/c1c-autocorrelation-certificate/raw
curl -s -f -L -o k22/K4194304_0.4103781.txt "$base/867e0e6f2858c0f59385cc663ac70c01bce8fd85/K4194304_0.4103781.txt"
curl -s -f -L -o k24/K16777216_0.4103782.txt.xz "$base/f2b97a3caa1fc2da7bd655024ce67658ff736bb4/K16777216_0.4103782.txt.xz"
python3 -I -c "import lzma, shutil, sys
with lzma.open(sys.argv[1]) as f, open(sys.argv[2], 'wb') as g: shutil.copyfileobj(f, g, 1 << 20)" \
  k24/K16777216_0.4103782.txt.xz k24/K16777216_0.4103782.txt
shasum -a 256 -c "$here/SHA256SUMS"
