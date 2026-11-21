#!/bin/sh
# Download the six colourings (data only) from the ancillary files of arXiv 2611.25143v1 into certs/ and check SHA-256.
set -e
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/certs" && cd "$here/certs"
for f in R453_37_points.json R463_70_points.json R553_98_points.json R564_200_points_rule.json R664_600_points_rule.json R665_83_points_short.json; do
  [ -s "$f" ] || { curl -s -f -L -o "$f" "https://arxiv.org/src/2611.25143v1/anc/$f"; sleep 3; }
done
shasum -a 256 -c "$here/SHA256SUMS"
