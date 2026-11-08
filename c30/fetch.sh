#!/bin/sh
# Download the LaTeX source of Norton, arXiv 2608.20292v3 (one paper.tex), and, for the candidate row of PR #317, its
# proof document and certificate JSON (rpinho/math-certificates @ 4f0bbf8, c30/round1: prose and data only, no code),
# each into its own new directory, and check their SHA-256 against SHA256SUMS. Usage: ./fetch.sh [CERTS_ROOT]
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d n2026 ] && [ -n "$(ls -A n2026)" ]; then echo "refusing: $root/n2026 is not empty"; exit 1; fi
mkdir -p n2026
curl -s -f -L -o n2026/eprint.tar.gz "https://arxiv.org/e-print/2608.20292v3"
(cd n2026 && tar xzf eprint.tar.gz && rm eprint.tar.gz)
if [ -d pin2026 ] && [ -n "$(ls -A pin2026)" ]; then echo "refusing: $root/pin2026 is not empty"; exit 1; fi
mkdir -p pin2026
P=https://raw.githubusercontent.com/rpinho/math-certificates/4f0bbf8d13183432dffd8addfd53811e830916c5/c30/round1
curl -s -f -L -o pin2026/FIRST_MOMENT_PROOF.md "$P/FIRST_MOMENT_PROOF.md"
curl -s -f -L -o pin2026/candidate_certificate.json "$P/artifacts/candidate_certificate.json"
shasum -a 256 -c "$here/SHA256SUMS"
