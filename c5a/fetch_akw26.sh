#!/bin/sh
# Download, into a new directory CERTS_ROOT/akw26, and check the SHA-256 of (data and text only; nothing here is run):
#   the twelve certificate files of John Akwei's Zenodo deposit 10.5281/zenodo.24138215 (v1.0, 11 Nov 2026),
#     sym_polya_m{128,256,512}_L8_{P,lam,W_end1}.txt and sym_polya_m1024_L6_{P,lam,W_end1}.txt: integers over 10^12
#     (64 MB; the m = 1024 boundary file alone is 41 MB);
#   the deposit's README_symmetric_certificates.md (the file format, the checks, the printed values) and its own
#     SHA256SUMS, read as text;
#   the V6 report, A_Data_Science_Analysis_of_Erdos_Problem_30_V6.pdf (johnakwei/Science @ f357556), read as text:
#     Appendix D, Theorem 5.4 and the text after it, Table 10.
# The deposit's programs (its two checkers, certifier and scripts) are not downloaded.
# Usage: ./fetch_akw26.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d akw26 ] && [ -n "$(ls -A akw26)" ]; then echo "refusing: $root/akw26 is not empty"; exit 1; fi
mkdir -p akw26
z="https://zenodo.org/api/records/24138215/files"
for f in README_symmetric_certificates.md SHA256SUMS \
         sym_polya_m1024_L6_P.txt sym_polya_m1024_L6_lam.txt sym_polya_m1024_L6_W_end1.txt \
         sym_polya_m128_L8_P.txt sym_polya_m128_L8_lam.txt sym_polya_m128_L8_W_end1.txt \
         sym_polya_m256_L8_P.txt sym_polya_m256_L8_lam.txt sym_polya_m256_L8_W_end1.txt \
         sym_polya_m512_L8_P.txt sym_polya_m512_L8_lam.txt sym_polya_m512_L8_W_end1.txt; do
  curl -s -f -L --retry 3 -o "akw26/$f" "$z/$f/content"
done
curl -s -f -L --retry 3 -o akw26/A_Data_Science_Analysis_of_Erdos_Problem_30_V6.pdf \
  "https://raw.githubusercontent.com/johnakwei/Science/f35755609cd36aad6bba7f315dca9e9765ba4441/A_Data_Science_Analysis_of_Erdos_Problem_30_V6.pdf"
shasum -a 256 -c "$here/SHA256SUMS_akw26"
