#!/bin/sh
# Download the data files replayed here (no code), each into its own new directory, and check their SHA-256:
# Griego's certificate (sebastian-griego/turan-c42-certificate, v1.0.0 = the commit the 42a page pins; against his
# .sha256 file), Pinho's 256-step profile (rpinho/c42-certificate at the commit the page pins; the value the page
# prints) and Koshta's tau = 41/200 certificate (DennisKoshta/c42-turan-power-sum-certificate at commit 8e1797e0; the
# value its note and SHA256SUMS state); and the three phase profiles of the Theorem A rows (thmA_check.py):
# Bakeas's (Bakeas/turan-c42-bound v1.1.0 = dc9c491c; the hash PR #302 and NOTES.md state), Pinho's
# (rpinho/c42-certificate at 009120a8, volterra-quadratic-knots/certificates; the hash its README and SHA256SUMS
# state) and Vatopoulos's (NikosVat/turan-c42 at 9668a60e, certificate/; the hash PR #327 states). Their checkers are
# neither downloaded nor run. A file already present is kept and checked again; a download goes only into a new or
# empty directory.
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
get() {  # directory, file, url
  [ -f "$1/$2" ] && return 0
  if [ -d "$1" ] && [ -n "$(ls -A "$1")" ]; then echo "refusing: $root/$1 is not empty"; exit 1; fi
  mkdir -p "$1"
  curl -s -f -L -o "$1/$2" "$3"
}
get griego turan42_certificate.json \
  "https://raw.githubusercontent.com/sebastian-griego/turan-c42-certificate/c8ddce14d9a5e898406d5dc6b8d08bb8a39507c7/certificate/turan42_certificate.json"
get pinho c42_256.json \
  "https://raw.githubusercontent.com/rpinho/c42-certificate/a8a9bc309bd1992413fa07595bd60c5e3fa1b65b/certificates/c42_256.json"
get koshta c42_tau41_200_M64_0.6873358063.json \
  "https://raw.githubusercontent.com/DennisKoshta/c42-turan-power-sum-certificate/8e1797e0d35454ad18cfebbbeb5858bf202b0748/certificates/c42_tau41_200_M64_0.6873358063.json"
get bakeas c42_profile_0.68637.json \
  "https://raw.githubusercontent.com/Bakeas/turan-c42-bound/dc9c491ccb0a9ee8de753739eb082fa84c25c091/c42_profile_0.68637.json"
get pinho_vqk C42_compact_best.json \
  "https://raw.githubusercontent.com/rpinho/c42-certificate/009120a811abd6840a9417e9f1c02cd490c5f608/volterra-quadratic-knots/certificates/C42_compact_best.json"
get vatopoulos c42_profile_0.686360751.json \
  "https://raw.githubusercontent.com/NikosVat/turan-c42/9668a60edf8877814c316a4bb1d901eff0c3988e/certificate/c42_profile_0.686360751.json"
shasum -a 256 -c "$here/SHA256SUMS"
