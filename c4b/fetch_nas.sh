#!/bin/sh
# Download the LaTeX sources of Naslund's two papers (enaslund/sarkozy-lower-bound-0.758: commit e5d6937 for the
# 0.75806746 row [Nas2026a], commit 5a83f41 for the 0.7580758 row [Nas2026b]), each into its own new directory, decode
# the archive embedded in each (data; its SHA-256 must be the one the paper prints) and keep only the four JSON witness
# files, each checked against the SHA-256 of that paper's appendix. No program of either archive is extracted or read.
# Usage: ./fetch_nas.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
for d in nas2026a nas2026b; do
  if [ -d $d ] && [ -n "$(ls -A $d)" ]; then echo "refusing: $root/$d is not empty"; exit 1; fi
  mkdir -p $d
done
base=https://raw.githubusercontent.com/enaslund/sarkozy-lower-bound-0.758
curl -s -f -L -o nas2026a/paper.tex "$base/e5d6937/papers/square-difference-free-sets-of-exponent-0.75806746.tex"
curl -s -f -L -o nas2026b/paper.tex "$base/5a83f41/papers/square-difference-free-sets-of-exponent-0.7580758318.tex"
extract() {  # extract TEX OUTDIR ARCHIVE_SHA256 MEMBER...
  python3 -I -c '
import base64, hashlib, io, os, sys, zipfile
s = open(sys.argv[1], encoding="utf-8").read()
i = s.index("PAYLOAD = \"\"\"") + len("PAYLOAD = \"\"\"")
raw = base64.b64decode("".join(s[i:s.index("\"\"\"", i)].split()))
assert hashlib.sha256(raw).hexdigest() == sys.argv[3], "the embedded archive is not the one the paper names"
z = zipfile.ZipFile(io.BytesIO(raw))
for name in sys.argv[4:]:
    open(os.path.join(sys.argv[2], os.path.basename(name)), "wb").write(z.read(name))
' "$@"
}
extract nas2026a/paper.tex nas2026a 5dcaffc0c826e178ad9b7ec3b42007f92eb4f7f44a03fe9c24012247f87b1e6b \
  research-notes/capacity-bound/bound.json research-notes/odd-search/q215-depth3-product.json \
  research-notes/odd-push/q437-depth3-final.json research-notes/binary-push/certificate-alpha-0758067460.json
extract nas2026b/paper.tex nas2026b 088f7876b76fd83745932beec88548df796ac02a2c3b66907608c79b721a8295 \
  research-notes/capacity-bound/bound.json research-notes/odd-search/q215-depth3-product.json \
  research-notes/next-push/q437-unequal32s.json research-notes/next-push/binary-cert-0.7580758318008816.json
shasum -a 256 -c "$here/SHA256SUMS_nas"
