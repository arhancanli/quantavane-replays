#!/usr/bin/env python3
"""For the note on the 1d page's [MV2010] row: Matolcsi and Vinuesa's signed example (arXiv 0907.1379v2, appendix,
n = 150, "||f*f||_inf = 1.45810...") gives, with c = v * v (the autoconvolution of the printed values, exact
rationals) and S = sum v, both
    max_m c_m * 2n / S^2      (max_t f*f / (int f)^2: the |max f*f| constant that the 1d page warns is different), and
    max_m |c_m| * 2n / S^2    (max_t |f*f| / (int f)^2: the 1d constant).
f*f is linear between the knots, so both maxima are attained at knots.  Exact rationals; the values are read as the
printed decimals.
usage: python3 -I signed_note.py FILE.json
"""
import json
import sys
from fractions import Fraction


def dec(fr, digits):
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    q = fr.numerator * 10 ** digits // fr.denominator
    s = str(q).rjust(digits + 1, "0")
    return sign + s[:-digits] + "." + s[-digits:]


d = json.load(open(sys.argv[1]), parse_float=Fraction)
v = [Fraction(x) for x in d["values"]]
n = len(v)
S = sum(v)
c = [sum(v[i] * v[m - i] for i in range(max(0, m - n + 1), min(m, n - 1) + 1)) for m in range(2 * n - 1)]
mx, mn = max(c), min(c)
k = Fraction(2 * n) / (S * S)
print(f"{sys.argv[1]}: n = {n}, {sum(1 for x in v if x < 0)} negative values, sum = {dec(S, 12)}")
print(f"  max_t f*f / (int f)^2   = {dec(mx * k, 15)}... (knot {c.index(mx)})")
print(f"  min_t f*f / (int f)^2   = {dec(mn * k, 15)}... (knot {c.index(mn)})")
print(f"  max_t |f*f| / (int f)^2 = {dec(max(mx, -mn) * k, 15)}...")
res = mx * k < Fraction(14582, 10000) and max(mx, -mn) * k > 2
print(f"RESULT note: the signed example gives {dec(mx * k, 7)}... for max f*f (the |max f*f| constant) but "
      f"{dec(max(mx, -mn) * k, 6)}... for max |f*f| (C_1d), so it is not a 1d certificate: {'CONFIRMED' if res else 'NOT CONFIRMED'}")
sys.exit(0 if res else 1)
