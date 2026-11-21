#!/bin/sh
# 49a: Tilley's 54-point block (PR #379) with two more sets in each copy of B; program A (Python, on our 4 November
# checker ../../c49/sunflower_check.py) and program B (JavaScript), then red controls. Python 3.9+, node 18+; seconds.
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
mkdir -p results
status=0
CLAIM=1.5678672030558807
python3 -I tilley_a.py ../../c49/sunflower_check.py --extra-left extras.json --extra-right extras.json --claim $CLAIM > results/A.txt 2>&1 && echo "A PASS" || { echo "A FAIL"; status=1; }
TILLEY_EXTRAS="$(cat extras.json)" node tilley_b.mjs extras > results/B.txt 2>&1 && grep -q 'VARIATION extras (2 cells): valid true' results/B.txt && echo "B PASS" || { echo "B FAIL"; status=1; }
grep -q 'rho^(1/52) in \[1.567867203055880727' results/B.txt || { echo "B rate differs"; status=1; }
# red controls: (1) a third set (Pinho's 64768, a seven-point set inside pendant 1) added: not valid in A or B;
# (2) the claimed decimal raised by 1e-16: not a valid lower bound
echo '[268173364, 261940, 64768]' > results/extras_control.json
if python3 -I tilley_a.py ../../c49/sunflower_check.py --extra-left results/extras_control.json --extra-right results/extras_control.json > results/control_A_third_set.txt 2>&1; then echo "control A third set NOT refused"; status=1; else echo "control A third set refused"; fi
if TILLEY_EXTRAS="$(cat results/extras_control.json)" node tilley_b.mjs extras > results/control_B_third_set.txt 2>&1; then echo "control B third set NOT refused"; status=1; else echo "control B third set refused"; fi
if python3 -I tilley_a.py ../../c49/sunflower_check.py --extra-left extras.json --extra-right extras.json --claim 1.5678672030558808 > results/control_A_claim_up.txt 2>&1; then echo "control A claim+1e-16 NOT refused"; status=1; else echo "control A claim+1e-16 refused"; fi
# third program: the referee's (written separately from the note); its own red controls are inside c49_run.py
(cd referee && python3 -s -E c49_run.py > ../results/referee_run.txt 2>&1 && python3 -s -E c49_rate.py > ../results/referee_rate.txt 2>&1) \
  && grep -q 'X\*+X.-both: struct_ok=True sunflower=False sperner_violation=False' results/referee_run.txt && echo "referee program PASS" \
  || { echo "referee program FAIL"; status=1; }
[ $status = 0 ] && echo "ALL PASS" || echo FAIL
exit $status
