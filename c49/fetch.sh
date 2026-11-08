#!/bin/sh
# Download the data files replayed here (no code), each into its own new directory, at the commits the 49a page pins,
# and check their SHA-256 against SHA256SUMS. Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # get DIR FILE REPO COMMIT PATH: one file into a directory that must be new or empty
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "https://raw.githubusercontent.com/$3/$4/$5"
}
IRSI=Aspect5/sunflower-free-capacity
get irsi_v10 cert_c49_chain20.json $IRSI 4d23c9ad391c9a087cad03a3c3346b5f79badae0 cert_c49_chain20.json        # tag v1.0
get irsi_v11 cert_c49_chain20_s17.json $IRSI f685bfab116e36568aacc9558eb582376764b88d cert_c49_chain20_s17.json  # tag v1.1
get irsi_v12 cert_c49_chain20_v12.json $IRSI e478c81d8b57eba956fdaf621643eb6f0781f820 cert_c49_chain20_v12.json  # tag v1.2
get irsi_degkm cert_degkm_control.json $IRSI 4d23c9ad391c9a087cad03a3c3346b5f79badae0 cert_degkm_control.json  # DEGKM control
CLA=joshua17-create/sunflower-lower-bound
get cla_k5 k5_family_masks.txt $CLA 6e34d78ab482f4cfd5a5ece62a6d06141e3e8c8d k5_family_masks.txt
get cla_z11 z11_family_masks.txt $CLA 6e34d78ab482f4cfd5a5ece62a6d06141e3e8c8d z11_family_masks.txt
SYL=xSyllaa/sunflower-pendant-chain
get syl_F F.txt $SYL 914de1553ebdc6621425c4cb21d847dc1e0f5af7 data/F.txt
get syl_H H.txt $SYL 914de1553ebdc6621425c4cb21d847dc1e0f5af7 data/H.txt
get syl_B B.txt $SYL 914de1553ebdc6621425c4cb21d847dc1e0f5af7 data/B.txt
get syl_Bprime Bprime.txt $SYL 8e21d0287711b4f0bfc23e9443523e2aa472cca1 data/Bprime.txt
get kos two_pendant_extra_1.5678599.txt.xz DennisKoshta/sunflower3-49a-certificate \
    c4ee81737f9f6a4ae9793586a27004553a12d5c0 two_pendant_extra_1.5678599.txt.xz
# the 49a page as PR #322 writes it (the row's new sentence on the per-pendant count)
get pr322 49a.md maximusrufus/optimizationproblems 6bed46a4ba1e0d99f37f82fc1a110db22c0c2cd2 constants/49a.md
# xz decompression with the standard library (the .txt is what the page's hash 10ee18f9...4838 refers to)
python3 -I -c "import lzma, sys; open(sys.argv[2], 'wb').write(lzma.open(sys.argv[1]).read())" \
    kos/two_pendant_extra_1.5678599.txt.xz kos/two_pendant_extra_1.5678599.txt
shasum -a 256 -c "$here/SHA256SUMS"
