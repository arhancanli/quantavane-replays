#!/bin/sh
# 50a, the lower bound [PT21] 0.533 (Parekh and Thompson, arXiv 2105.05698, Theorem 2): the two product-state ratios of
# Lemma 18 (the Gharibian-Parekh ratio minimised over [-1, 1/3) >= 0.498766 and over [-5/9, 1/3) >= 0.557931) certified
# by two separately written programs (c50apt_a.py: hypergeometric series in exact rationals; c50apt_b.py: closed forms
# through arcsin in decimal intervals; both copied from ../c50a with a range option), then the closing linear program
# in exact rationals (c50apt_lp.py). Three red controls. Nothing is downloaded; no program of any source read or run.
# Usage: ./run_all.sh   Python 3.9+; about a minute.
here=$(cd "$(dirname "$0")" && pwd)
out="$here/results/c50apt.txt"
mkdir -p "$here/results"
{
  echo "== program A: min over [-1, 1/3) >= 0.498766 (B of Lemma 18) and over [-5/9, 1/3) >= 0.557931 (A)"
  python3 -I "$here/c50apt_a.py" --lo -1 --target 0.498766 | grep "cells\|RESULT\|VERDICT"
  python3 -I "$here/c50apt_a.py" --lo -5/9 --target 0.557931 | grep "cells\|RESULT\|VERDICT"
  echo "== program B: the same two minima"
  python3 -I "$here/c50apt_b.py" --lo -1 --target 0.498766 --cells 4000 | grep "cells\|RESULT\|VERDICT"
  python3 -I "$here/c50apt_b.py" --lo -0.5555555555555555555555555555555555556 --target 0.557931 --cells 2000 | grep "cells\|RESULT\|VERDICT"
  echo "== the linear program of the proof of Theorem 2"
  python3 -I "$here/c50apt_lp.py"
  echo "== red control 1: the full-range minimum claimed >= 0.498768 (it is 0.4987674...), A then B"
  python3 -I "$here/c50apt_a.py" --lo -1 --target 0.498768 | grep "does not\|RESULT\|VERDICT"
  python3 -I "$here/c50apt_b.py" --lo -1 --target 0.498768 --cells 4000 | grep "does not\|RESULT\|VERDICT"
  echo "== red control 2: the small-edge minimum claimed >= 0.557933 (it is 0.5579321...), A then B"
  python3 -I "$here/c50apt_a.py" --lo -5/9 --target 0.557933 | grep "does not\|RESULT\|VERDICT"
  python3 -I "$here/c50apt_b.py" --lo -0.5555555555555555555555555555555555556 --target 0.557933 --cells 2000 | grep "does not\|RESULT\|VERDICT"
  echo "== red control 3: the bound claimed as 0.5331"
  python3 -I "$here/c50apt_lp.py" 0.5331
} > "$out" 2>&1
pass=$(grep -c "VERDICT PASS" "$out"); ref=$(grep -c "VERDICT REFUSED" "$out"); bad=$(grep -c -i "traceback\|error" "$out")
if [ "$pass" = 5 ] && [ "$ref" = 5 ] && [ "$bad" = 0 ]; then
  echo "RESULT: ALL CHECKS PASS (both programs: min R >= 0.498766 on [-1, 1/3) and >= 0.557931 on [-5/9, 1/3); the linear program gives r* = 231411/434165 = 0.5330024... >= 0.533; red controls refused)" >> "$out"; tail -1 "$out"; exit 0
else
  echo "RESULT: FAILED (pass=$pass refused=$ref errors=$bad)" >> "$out"; tail -1 "$out"; exit 1
fi
