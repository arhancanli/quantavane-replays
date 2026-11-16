#!/usr/bin/env python3
"""13b, lower bound [Elek1994] 0.8257 ("use unit disk and equilateral triangle"), program A: the argument and its
numbers in exact rational interval arithmetic.

A convex set covering a congruent copy of every set of diameter 1 contains a disk D of diameter 1 (radius r = 1/2)
and an equilateral triangle T of side 1, so its area is at least m(v) = area conv(D u (T + v)) for some v (D is
rotation invariant, so the triangle's orientation may be fixed; centres: D at 0, T's centroid at v). We show
min_v m(v) >= pi/8 + sqrt(3)/4 = 0.8257117... >= 0.8257:
  - for a point p with |p| = d > r, conv(D u {p}) \\ D is a cap of area g(d) = r sqrt(d^2 - r^2) - r^2 arccos(r/d)
    (two tangent triangles minus the sector), lying in the angular sector of half-angle arccos(r/d) around p;
    g = 0 for d <= r; g'(d) = r sqrt(d^2 - r^2)/d >= 0 and g''(d) = r^3/(d^2 sqrt(d^2 - r^2)) > 0, so g is convex
    and non-decreasing on [0, oo), and L(v) = pi r^2 + sum_k g(|q_k + v|) (q_k the vertices of T) is convex in v;
  - L(Rv) = L(v) for the rotation R by 120 degrees, so L(v) >= L((v + Rv + R^2 v)/3) = L(0) for every v;
  - case (a), all |q_k + v| <= d1 = 7/10: two vertices p, p' with |p - p'| = 1 and |p|, |p'| <= d1 <= 1/sqrt 2 see
    each other from 0 at an angle gamma with cos gamma = (|p|^2 + |p'|^2 - 1)/(2|p||p'|) <= 1 - 1/(2 d1^2)
    <= 1/(2 d1^2) - 1 = cos(2 arccos(r/d1)) (the left side increases and the right side decreases in |p|, |p'|),
    so the caps' sectors are disjoint, the three caps are disjoint, and m(v) >= L(v) >= L(0);
  - case (b), some |q_k + v| > d1: m(v) >= area conv(D u {q_k + v}) >= pi r^2 + g(d1);
  - L(0) = pi/4 + 3 g(1/sqrt 3) = pi/4 + 3 (sqrt(3)/12 - pi/24) = pi/8 + sqrt(3)/4.
Checked here by exact rationals (pi by Machin's formula with alternating tails; square roots by integer square roots;
arccos(5/7) = arctan(sqrt(24)/5) by halving and the alternating series): the identity L(0) = pi/8 + sqrt(3)/4 in
intervals, d1^2 <= 1/2, pi/4 + g(d1) > pi/8 + sqrt(3)/4, and the printed 0.8257 <= pi/8 + sqrt(3)/4 (a truncation).
Red controls (--control): the bound claimed as 0.8258; d1 = 3/4 (above 1/sqrt 2, so case (a) is not proved); a
triangle of side 0.9 in place of 1 (L(0) changes; 0.8257 then fails). Standard library only.
"""
import sys
from fractions import Fraction as F
from math import isqrt

CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
K = 40   # digits for square roots


def atan_series(x, terms):
    """arctan of a rational 0 <= x <= 1/4, enclosed by consecutive partial sums"""
    s = F(0); p = x
    for k in range(terms):
        s += (-1) ** k * p / (2 * k + 1)
        p *= x * x
    nxt = p / (2 * terms + 1)
    return (s, s + nxt) if terms % 2 == 0 else (s - nxt, s)


def sqrt_iv(q):
    n_lo = q.numerator * 10 ** (2 * K) // q.denominator
    n_hi = -((-q.numerator * 10 ** (2 * K)) // q.denominator)
    lo = isqrt(n_lo); hi = isqrt(n_hi)
    if hi * hi < n_hi: hi += 1
    return F(lo, 10 ** K), F(hi, 10 ** K)


def atan_iv(lo, hi):
    """arctan on [lo, hi] (0 <= lo <= hi), monotone: halve via atan x = 2 atan(x/(1 + sqrt(1 + x^2)))"""
    out = []
    for x, side in ((lo, 0), (hi, 1)):
        mult = 1
        while x > F(1, 4):
            s = sqrt_iv(1 + x * x)
            x = x / (1 + s[0]) if side else x / (1 + s[1])   # larger x for the upper end, smaller for the lower
            mult *= 2
        a = atan_series(x, 40)
        out.append(a[side] * mult)
    return out[0], out[1]


a5 = atan_series(F(1, 5), 40); a239 = atan_series(F(1, 239), 20)
PI = (16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0])


def arccos_iv(xl, xh):
    """arccos on [xl, xh] in (0, 1): arccos x = arctan(sqrt(1 - x^2)/x), decreasing in x"""
    s_lo = sqrt_iv(1 - xh * xh)[0]; s_hi = sqrt_iv(1 - xl * xl)[1]
    return atan_iv(s_lo / xh, s_hi / xl)


def g_iv(d, r):
    """cap area g(d) for a rational d > r: r sqrt(d^2 - r^2) - r^2 arccos(r/d)"""
    s = sqrt_iv(d * d - r * r)
    ac = arccos_iv(r / d, r / d)
    return r * s[0] - r * r * ac[1], r * s[1] - r * r * ac[0]


def dec(q, k=12):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def main():
    assert F(314159265358979323846, 10 ** 20) < PI[0] and PI[1] < F(314159265358979323847, 10 ** 20)
    r = F(1, 2)
    side = F(9, 10) if CONTROL == "side" else F(1)
    d1 = F(3, 4) if CONTROL == "d1" else F(7, 10)
    printed = F("0.8258") if CONTROL == "claim" else F("0.8257")
    ok = True
    # L(0) = pi r^2 + 3 g(circumradius), circumradius side/sqrt 3: g at an irrational d, so enclose d^2 = side^2/3 exactly:
    # g(d) = r sqrt(d^2 - r^2) - r^2 arccos(r/d), with r/d = r sqrt 3/side enclosed
    d2 = side * side / 3
    sq = sqrt_iv(d2 - r * r)
    rd2 = r * r / d2                          # (r/d)^2
    rdl, rdh = sqrt_iv(rd2)
    ac = arccos_iv(rdl, rdh)
    g0 = (r * sq[0] - r * r * ac[1], r * sq[1] - r * r * ac[0])
    L0 = (PI[0] * r * r + 3 * g0[0], PI[1] * r * r + 3 * g0[1])
    s3 = sqrt_iv(F(3))
    closed = (PI[0] / 8 + s3[0] / 4, PI[1] / 8 + s3[1] / 4)
    same = side == 1 and abs(L0[0] - closed[0]) < F(1, 10 ** 30) and abs(L0[1] - closed[1]) < F(1, 10 ** 30)
    print("triangle side %s: L(0) = pi/4 + 3 g(side/sqrt 3) in [%s, %s]" % (side, dec(L0[0], 30), dec(L0[1], 30)))
    print("pi/8 + sqrt(3)/4 in [%s, %s]; the same number: %s" % (dec(closed[0], 30), dec(closed[1], 30), same))
    # the triangle's vertices must lie outside D and its edge midpoints inside (the cap picture at v = 0)
    pic = d2 > r * r and side * side / 12 < r * r
    print("at v = 0: vertices outside D (d^2 = %s > 1/4): %s; edge midpoints inside D (side^2/12 < 1/4): %s" % (d2, d2 > r * r, side * side / 12 < r * r))
    case_a = d1 * d1 <= F(1, 2) and d1 > r
    print("case (a): d1 = %s, d1^2 <= 1/2 (so the caps' sectors are disjoint): %s" % (d1, case_a))
    gb = g_iv(d1, r)
    case_b = PI[0] / 4 + gb[0] > L0[1]
    print("case (b): pi/4 + g(d1) >= %s > L(0): %s" % (dec(PI[0] / 4 + gb[0], 12), case_b))
    bound = L0[0]
    holds = printed <= bound
    tr = holds and bound < printed + F(1, 10 ** 4)
    print("so m(v) >= %s... for every v; printed lower bound %s: %s (truncation: %s)"
          % (dec(bound, 12), printed, "holds" if holds else "FAILS", tr))
    ok = same and pic and case_a and case_b and holds and tr
    print("key %s" % dec(L0[0], 25))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
