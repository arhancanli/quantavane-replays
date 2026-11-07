#!/bin/sh
# Download the certificates replayed here, at the commits their pull requests pin, and check their SHA-256.
set -e
mkdir -p certs && cd certs
get() { curl -s -f -L -o "$1" "$2"; }
get koshta_K24576_0.4103333.json https://raw.githubusercontent.com/DennisKoshta/c1c-autocorrelation-certificate/a4ee94ef29ff/K24576_0.4103333.json
get julsgaard_c1c_K98304.json https://raw.githubusercontent.com/Julsgaard/c1c-certificate/91e47b1aa03020fb4f9219ad16f2f8cfcec9cd1d/cert/certificate.json
get vidigoat_c1c_K98304_v2.json https://raw.githubusercontent.com/vidigoat/optimization-records/9f82508a30d74b6b3c59c6f495269bd3bf4f61e5/tao/certificates/c1c_v2/c1c_K98304_v2.json
get vidigoat_c1c_K98304_v3.json https://raw.githubusercontent.com/vidigoat/optimization-records/42d6a89fab4b1714af3a98a6dd43368615052b45/tao/certificates/c1c_v3/c1c_K98304_v3.json
get vidigoat_c1a_n131072.json https://raw.githubusercontent.com/vidigoat/optimization-records/9f82508a30d74b6b3c59c6f495269bd3bf4f61e5/tao/certificates/c1a/c1a_n131072.json
get vidigoat_c1a_n262144.json https://raw.githubusercontent.com/vidigoat/optimization-records/9f82508a30d74b6b3c59c6f495269bd3bf4f61e5/tao/certificates/c1a/c1a_n262144.json
get julsgaard_certificate.json https://raw.githubusercontent.com/Julsgaard/c1d-certificate/cc455ca0989f40b2b2a8a5248f507e9c78e95383/cert/certificate.json
get kerem_c1d_n40960_1.4556234448.json https://raw.githubusercontent.com/Keremozdemirra/tao-1c-1d-certificates/408f8d0b6e8d232e632354997f1cc7f223f93ff9/c1d_n40960_1.4556234448.json
get echul_c1d_records.tgz https://raw.githubusercontent.com/echul-seal/c1d-upper-bound/2be42dfdb6ecc60e2fe1345a38c4db21688f5726/c1d_records.tgz
get huiihao_construction.json https://raw.githubusercontent.com/huiihao/c1d-certificate/302817dc3227ce281e9cfb486f12ff0df7fe2876/construction.json
python3 -I -c "import tarfile; t = tarfile.open('echul_c1d_records.tgz'); m = t.getmember('ub/c_f10240.json'); assert m.isfile(); open('echul_c_f10240.json', 'wb').write(t.extractfile(m).read())"
# PR #326's two 1d rows: STEPS2026 (construction.json) and AER2026 (a Python file, parsed with ast as text and never run;
# its N = 600 list is copied verbatim into aer2026_n600.json)
get steps_construction.json https://raw.githubusercontent.com/steps-re/c1d-certificate/9a091494facc801eaf6b92b175b63c58b7c3b30b/construction.json
get aer_problem_4_signed_autoconvolution_bounds.py https://raw.githubusercontent.com/unclegareth-create/alphaevolve_repository_of_problems/4b3bbd2a44db92442e4ebb4e64e0bb90f229aad1/experiments/autocorrelation_problems/problem_4_signed_autoconvolution_bounds.py
python3 -I -c "import ast; src = open('aer_problem_4_signed_autoconvolution_bounds.py').read(); n = next(a for a in ast.parse(src).body if isinstance(a, ast.Assign) and getattr(a.targets[0], 'id', '') == 'STEP_FUNCTION_C4_PRIME_N600'); t = [ast.get_source_segment(src, e).strip() for e in n.value.elts]; assert len(t) == 600; open('aer2026_n600.json', 'w').write('{\"values\": [' + ', '.join(t) + ']}\\n')"
# the starred 1d row [CERT2026] (#224): the Extremal Values Registry copy the page cites (EV-000004, LF line endings)
get cert2026_construction.json https://raw.githubusercontent.com/extremal-values-registry/extremal-values-registry/50520b971e5d234bc9c4ccbdcd5a5107d5bb8791/certificates/third_autocorrelation_abs_inside/construction.json
shasum -a 256 -c ../SHA256SUMS
