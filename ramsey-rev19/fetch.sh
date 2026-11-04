#!/bin/sh
# Download the witness files replayed here (data only, no code), each source into its own new directory, from the
# authors' own channels at pinned commits or records, and check their SHA-256 against SHA256SUMS.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
fresh() { if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi; mkdir -p "$1"; }
gh() {  # gh DIR FILE REPO COMMIT PATH
  curl -s -f -L -o "$1/$2" "https://raw.githubusercontent.com/$3/$4/$5"
}
# Nagda, Raghavan, Thakurta (arXiv:2603.09172): google-research/google-research, ramsey_number_bounds/improved_bounds,
# at 0d94f8c (22 Oct 2026, the folder's latest commit). Only the bounds that revision 18 does not list.
fresh nart
G=google-research/google-research; C=0d94f8c4d92678ce5f7a0eacae725f79174cfb03; P=ramsey_number_bounds/improved_bounds
gh nart R3_17_93.txt $G $C "$P/R%283%2C%2017%29%20%3E%3D%2093.txt"
gh nart R4_15_163.txt $G $C "$P/R%284%2C%2015%29%20%3E%3D%20163.txt"
gh nart R4_16_174.txt $G $C "$P/R%284%2C%2016%29%20%3E%3D%20174.txt"
gh nart R4_18_209.txt $G $C "$P/R%284%2C%2018%29%20%3E%3D%20209.txt"
gh nart R4_19_219.txt $G $C "$P/R%284%2C%2019%29%20%3E%3D%20219.txt"
gh nart R4_20_237.txt $G $C "$P/R%284%2C%2020%29%20%3E%3D%20237.txt"
gh nart R4_21_244.txt $G $C "$P/R%284%2C21%29%3E%3D244.txt"
gh nart R7_13_513.txt $G $C "$P/R%287%2C%2013%29%20%3E%3D%20513.txt"
gh nart R8_12_511.txt $G $C "$P/R%288%2C%2012%29%20%3E%3D%20511.txt"
gh nart R8_18_1052.txt $G $C "$P/R%288%2C%2018%29%20%3E%3D%201052.txt"
gh nart R9_11_617.txt $G $C "$P/R%289%2C%2011%29%20%3E%3D%20617.txt"
gh nart R10_12_856.txt $G $C "$P/R%2810%2C%2012%29%20%3E%3D%20856.txt"
# Yu (arXiv:2608.18169): the certificate in the paper's data record, Zenodo 21941603.
fresh yu
curl -s -f -L -o yu/CERTIFICATE_R4_20.txt "https://zenodo.org/api/records/21941603/files/CERTIFICATE_R4_20.txt/content"
# Coniglio, Furini, Ljubic, San Segundo, Thuerauf, Traversi (arXiv:2608.18769): fabiofurini/ramsey-number-lower-bounds
# at 0b53ca4 (2 Nov 2026), lower-bounds/certificates.
fresh coniglio
G=fabiofurini/ramsey-number-lower-bounds; C=0b53ca478bec9f9aa6c64e33c8b9b0088bafd9e8
for f in 24_graphSize_150 25_graphSize_159 26_graphSize_162 28_graphSize_180 29_graphSize_193 30_graphSize_198 \
         31_graphSize_209 32_graphSize_218 33_graphSize_229 34_graphSize_238 35_graphSize_252 36_graphSize_260 \
         37_graphSize_271 38_graphSize_280 39_graphSize_293 40_graphSize_304 41_graphSize_316 42_graphSize_326 \
         43_graphSize_339 44_graphSize_347 45_graphSize_363 46_graphSize_370 47_graphSize_384 48_graphSize_390 \
         49_graphSize_410; do
  gh coniglio "R_3_$f.clq" $G $C "lower-bounds/certificates/R_3_$f.clq"
done
# Cremer (arXiv:2610.12122): fritzcremer/ramsey-lower-bounds, release v1.1.0 (880b4fa), which the paper cites.
fresh cremer
G=fritzcremer/ramsey-lower-bounds; C=880b4fa5e683de842e514054da6edb8a5482ff4a
gh cremer n134-gpt-6.1-sol.txt $G $C R6_8/n134-gpt-6.1-sol.txt
gh cremer n344-gpt-6.1-sol-1.txt $G $C R8_10/n344-gpt-6.1-sol-1.txt
gh cremer n344-gpt-6.1-sol-2.txt $G $C R8_10/n344-gpt-6.1-sol-2.txt
# Tatarevic (arXiv:2608.06531): milostatarevic/ramsey-increment-theorem at 36c15ec, bounds/r1212_1640.graph.
fresh tatarevic
gh tatarevic r1212_1640.graph milostatarevic/ramsey-increment-theorem 36c15eca0c601be7b7f33e30ee6890f395d2e1ea \
   bounds/r1212_1640.graph
# Canli (Zenodo 23857057 and 23920916 point here): arhancanli/quantavane-ramsey-records at af0c27a.
fresh quantavane
G=arhancanli/quantavane-ramsey-records; C=af0c27aeb4eed353a8e3418174b7f4bdd3a4dbc1
gh quantavane R5_15_K278.txt $G $C colourings/R5_15_K278.txt
gh quantavane R8_12_K518.txt $G $C colourings/R8_12_K518.txt
shasum -a 256 -c "$here/SHA256SUMS"
