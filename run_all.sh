#!/bin/sh
# Replay every certificate fetched by fetch.sh; one output file per certificate in results/.
set -e
mkdir -p results
for job in "1c koshta_K24576_0.4103333.json" "1c julsgaard_c1c_K98304.json" "1c vidigoat_c1c_K98304_v2.json" "1c vidigoat_c1c_K98304_v3.json" \
           "1a vidigoat_c1a_n131072.json" "1a vidigoat_c1a_n262144.json" \
           "1d julsgaard_certificate.json" "1d kerem_c1d_n40960_1.4556234448.json" \
           "1d echul_c_f10240.json" "1d huiihao_construction.json"; do
  set -- $job
  python3 -I autocorr_check.py "$1" "certs/$2" --controls > "results/$2.txt"
done
python3 -I autocorr_check.py 1d certs/huiihao_construction.json --decimal > results/huiihao_construction.json.decimal.txt
tail -n +1 results/*.txt
