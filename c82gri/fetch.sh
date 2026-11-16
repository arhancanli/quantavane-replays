#!/bin/sh
# Download the C82 certificate of the 82a row [Gri26] (data only, no code) and check its SHA-256:
# maaxgrin/zhang-zagier-c82-bound @ cb88f92, runs/qstar_targeted_refined.jsonl (the refined records) and
# runs/qstar_cert_roots_262144_target253635.jsonl (the certified run, which names its candidate). The certificate is
# line 1 of the first file, the record whose label and bound_log equal the certified run's candidate_label and
# candidate_bound_log; its 'denominator' and 'numerator' fields are written to certificate.json (checked below).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -d gri26 ] && [ -n "$(ls -A gri26)" ]; then echo "refusing: $root/gri26 is not empty"; exit 1; fi
mkdir -p gri26
C=cb88f92434f649a1588b31ae91f67852830c6f77
for f in qstar_targeted_refined.jsonl qstar_cert_roots_262144_target253635.jsonl; do
  curl -s -f -L -o gri26/$f https://raw.githubusercontent.com/maaxgrin/zhang-zagier-c82-bound/$C/runs/$f
done
python3 -I -c "
import json
lines = open('gri26/qstar_targeted_refined.jsonl').read().splitlines()
d = json.loads(lines[0]); c = json.loads(open('gri26/qstar_cert_roots_262144_target253635.jsonl').read().splitlines()[0])
assert d['label'] == c['candidate_label'] and d['bound_log'] == c['candidate_bound_log'] and c['certified'] is True
json.dump({'denominator': d['denominator'], 'numerator': d['numerator']}, open('gri26/certificate.json', 'w'))
"
shasum -a 256 -c "$here/SHA256SUMS"
