#!/bin/sh
# Download the published step functions of eleven 1a and 1b rows, each from the source its row cites, into a new
# directory; copy their heights out as text (extract_values.py: ast, regex and json only; nothing downloaded is run);
# check SHA-256 of the files (and of the two values files cut from the live EinsteinArena responses, fetched 8 Nov 2026).
# Usage: ./fetch.sh [CERTS_ROOT]   (default: ./certs)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
mkdir -p "$root" && cd "$root"
if [ -n "$(ls -A . 2>/dev/null | grep -v ruvector.db)" ]; then echo "refusing: $root is not empty"; exit 1; fi
get() { curl -s -f -L -o "$1" "$2"; }
TOG=https://raw.githubusercontent.com/togethercomputer/EinsteinArena-new-SOTA/c388c6f7408c886311940896713339a1a70c2394
get tog_c1a_ours_2026.json $TOG/first-autocorrelation/solutions/ours_2026.json
get tog_c1b_together_ai_2026.py $TOG/erdos-minimum-overlap/solutions/together_ai_2026.py
get dm_mathematical_results.ipynb https://raw.githubusercontent.com/google-deepmind/alphaevolve_results/4226acbf237ff9ad10ba7673a2af127a2d8a5971/mathematical_results.ipynb
get aer_autocorrelation_problems.ipynb https://raw.githubusercontent.com/google-deepmind/alphaevolve_repository_of_problems/8f447457957deac61e28bf1676746f0753b3b2f8/experiments/autocorrelation_problems/autocorrelation_problems.ipynb
TTT=https://raw.githubusercontent.com/test-time-training/discover/6c40e82dab9d5de7416ac873ad5cd3106084aaed/results/mathematics
get ttt_ac1_sequence.json $TTT/ttt_ac1_sequence.json
get ttt_erdos_sequence.json $TTT/ttt_erdos_sequence.json
get atris_c1b_construction.json https://raw.githubusercontent.com/atrislabs/records/1fa857a021ddbb4064708daf336758b866f75cb4/1b-erdos-minimum-overlap/construction.json
get haugland_1609.08000v1.gz https://arxiv.org/e-print/1609.08000v1
gunzip -c haugland_1609.08000v1.gz > haugland_mop1.tex
get arena_problem1.json "https://einsteinarena.com/api/solutions/best?problem_id=1"
get arena_problem2.json "https://einsteinarena.com/api/solutions/best?problem_id=2"
# [PBV2026] (the 1a lower bound 1.292): the certificate JSON, read as data for its 200 multiplier coefficients only
get pbv2026_multiscale_arcsine_1292.json https://raw.githubusercontent.com/AndreiPiterbarg/sidon-autocorrelation/1ed5889100617ee086db1ce3456b75c1f8e5f6be/delsarte_dual/grid_bound_alt_kernel/certificates/multiscale_arcsine_1292.json
X="python3 -I $here/extract_values.py"
$X py tog_c1b_together_ai_2026.py h_values c1b_together_ai_2026.json
$X cell dm_mathematical_results.ipynb 53 step_function_heights_1 c1a_alphaevolve_2025.json
$X cell dm_mathematical_results.ipynb 78 best_sequence ae_reflect c1b_alphaevolve_2025.json
$X cell aer_autocorrelation_problems.ipynb 60 best_sequence c1a_alphaevolve_dec2025.json
$X tex51 haugland_mop1.tex c1b_haugland_2016.json
$X arena arena_problem2.json 2603 c1a_ea2026_2603.json
$X arena arena_problem1.json 2507 c1b_ea2026_2507.json
shasum -a 256 -c "$here/SHA256SUMS"
