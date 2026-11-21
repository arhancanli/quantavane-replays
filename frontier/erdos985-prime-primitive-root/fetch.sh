#!/bin/sh
# Download the finite data of Howe, Huffman, Kohlmeyer and Kruer's Erdos problem 985 repository (data only, no code)
# into a new directory and check the SHA-256: H0W4RD-dev/erdos-985-prime-primitive-root @ 99fd0a7,
# small-prime-witnesses.tsv (Lemma 2.1, the 405 odd primes p <= 2791) and x-equals-p-witnesses.zip (the companion
# certificate, 7 strips of 16-bit words), unpacked with Python's zipfile module.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d erdos985 ] && [ -n "$(ls -A erdos985)" ]; then echo "refusing: $root/erdos985 is not empty"; exit 1; fi
mkdir -p erdos985/witnesses
base=https://raw.githubusercontent.com/H0W4RD-dev/erdos-985-prime-primitive-root/99fd0a70c602c8d501e35035e9a1a2bfbc8854aa
curl -s -f -L -o erdos985/small-prime-witnesses.tsv "$base/small-prime-witnesses.tsv"
curl -s -f -L -o erdos985/x-equals-p-witnesses.zip "$base/x-equals-p-witnesses.zip"
python3 -I -c "
import zipfile, sys, os
z = zipfile.ZipFile(sys.argv[1])
for name in z.namelist():
    if '/' in name or not name.startswith('strip-') or not name.endswith('.u16'):
        raise SystemExit('unexpected member ' + name)
    open(os.path.join(sys.argv[2], name), 'wb').write(z.read(name))
" erdos985/x-equals-p-witnesses.zip erdos985/witnesses
shasum -a 256 -c "$here/SHA256SUMS"
