#!/usr/bin/env python3
"""50a [PT21] 0.533: the last step of Parekh and Thompson's Theorem 2 (arXiv 2105.05698, proof of Theorem 2, read
through the HTML rendering): with the product-state ratios A >= 0.557931 on small edges (the Gharibian-Parekh ratio
R(v) = (1 - F_3(v))/(1 - 3v) minimised over v in [-5/9, 1/3), Lemma 18) and B >= 0.498766 on large edges (the same ratio
over [-1, 1/3)), both certified by c50apt_a.py and c50apt_b.py, and the matching ratios 3/8 and 3/4 of Lemma 17 (read),
the algorithm keeps the better state, so its ratio is at least
    min over s in [0, 1] of max{ (3/8) s + (3/4)(1 - s), A s + B (1 - s) }.
The first line decreases and the second increases (A > B), so the minimum is at their crossing:
    s* = (3/4 - B)/(A - B + 3/8),  r* = 3/4 - (3/8) s*.
Computed in exact rationals and again in 50-digit decimals; r* is increasing in A and in B (A > 3/8), so the certified
lower bounds A >= 0.557931, B >= 0.498766 give r >= r*. Usage: c50apt_lp.py [CLAIM]   (default 0.533)
"""
import sys
from fractions import Fraction as Q
from decimal import Decimal as D, getcontext

getcontext().prec = 50
claim = sys.argv[1] if len(sys.argv) > 1 else "0.533"
A, B = Q("0.557931"), Q("0.498766")
assert A > B and A > Q(3, 8)
s = (Q(3, 4) - B) / (A - B + Q(3, 8)); r = Q(3, 4) - Q(3, 8) * s
assert 0 <= s <= 1 and r == A * s + B * (1 - s)
Ad, Bd = D("0.557931"), D("0.498766")
sd = (D("0.75") - Bd) / (Ad - Bd + D("0.375")); rd = D("0.75") - D("0.375") * sd
agree = abs(D(r.numerator) / D(r.denominator) - rd) < D(10) ** -40
ok = r >= Q(claim) and agree
print("s* = %s, r* = %s = %s (decimal route agrees: %s); the row's %s <= r*: %s"
      % (s, r, format(D(r.numerator) / D(r.denominator), ".15f"), agree, claim, r >= Q(claim)))
print("VERDICT " + ("PASS" if ok else "REFUSED"))
