#!/usr/bin/env python3
"""1c, upper bound [BS2020]: C_1c <= 1/(2(1 + theta_0)) = 0.4107675..., theta_0 = -min_x sin(x)/x = 0.2172336...
Program A: the critical point, in exact rational interval arithmetic.

The row (constants/1c.md) gives the bound in closed form through theta_0 (the exact form is [MR2020, (1.2)]; the bound
itself is a theorem of Barnard and Steinerberger, read, not checked). We compute theta_0 from its definition:
  - sin(x)/x is even and equals 1 at 0; it is >= 0 on (0, pi]; for x >= 2 pi, sin(x)/x >= -1/x >= -1/(2 pi) > -0.16;
  - on (pi, 2 pi), (sin x/x)' = -f(x)/x^2 with f(x) = sin x - x cos x; f' = x sin x < 0 there, f(pi) = pi > 0 and
    f(3 pi/2) = -1 < 0, and f < 0 on [3 pi/2, 2 pi) (sin < 0 <= cos there... cos >= 0 on [3 pi/2, 2 pi)), so the only
    critical point is the root x_0 of f in (pi, 3 pi/2), and the minimum over x > 0 is sin(x_0)/x_0 (checked < -0.16);
  - at x_0, sin x_0 = x_0 cos x_0, so theta_0 = -cos x_0 = 1/sqrt(1 + x_0^2) (cos x_0 < 0), decreasing in x_0.
x_0 is enclosed by bisection on the sign of f (evaluated in exact rational intervals: sin and cos by Taylor series with
the Lagrange remainder, pi by Machin's formula with the alternating tails), then theta_0 and 1/(2(1 + theta_0)) are
enclosed with exact rational square roots (integer square roots of scaled numerators). The page prints theta_0 =
0.2172336... (a truncation) and the bound as 0.4107675...; the bound's value is 0.41076748818..., so that figure is a
correct rounding upward (valid for an upper bound), not its leading digits. Red controls (--control): theta_0 claimed
as 0.2172337...; the bound claimed as 0.4107674 (below the value); theta_0 taken as -max sin(x)/x. Standard library only.
"""
import sys
from fractions import Fraction as F
from math import isqrt

CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None


def atan_inv(n, terms):
    """arctan(1/n) enclosed: alternating series, partial sums bracket the value"""
    s = F(0); p = F(1, n); k = 0
    lo = hi = None
    for k in range(terms):
        s += (-1) ** k * p / (2 * k + 1)
        p /= n * n
    nxt = p / (2 * terms + 1)
    return (s, s + nxt) if terms % 2 == 0 else (s - nxt, s)


a5 = atan_inv(5, 40); a239 = atan_inv(239, 20)
PI = (16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0])


def sincos(x, terms=60):
    """sin and cos of a rational x with |x| < 6: Taylor to order 2*terms, remainder |x|^(n+1)/(n+1)!"""
    s = F(0); c = F(0); t = F(1)
    for n in range(2 * terms + 2):
        if n > 0:
            t = t * x / n
        if n % 4 == 0: c += t
        elif n % 4 == 1: s += t
        elif n % 4 == 2: c -= t
        else: s -= t
    r = abs(t * x / (2 * terms + 2))
    return (s - r, s + r), (c - r, c + r)


def f_interval(x):
    (s0, s1), (c0, c1) = sincos(x)
    lo = s0 - max(x * c0, x * c1)
    hi = s1 - min(x * c0, x * c1)
    return lo, hi


def sqrt_lo(q, k=60):
    n = q.numerator * 10 ** (2 * k) // q.denominator
    return F(isqrt(n), 10 ** k)


def sqrt_hi(q, k=60):
    n = -((-q.numerator * 10 ** (2 * k)) // q.denominator)
    r = isqrt(n)
    if r * r < n:
        r += 1
    return F(r, 10 ** k)


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def trunc_or_round(q_lo, q_hi, printed):
    d = F(printed); u = F(1, 10 ** (len(printed.split(".")[1])))
    tr = d <= q_lo and q_hi < d + u
    rd = abs(q_lo - d) <= u / 2 and abs(q_hi - d) <= u / 2
    return tr, rd


def main():
    ok = True
    pi_lo, pi_hi = PI
    assert F(314159265358979323846, 10 ** 20) < pi_lo and pi_hi < F(314159265358979323847, 10 ** 20)
    lo, hi = pi_hi, 3 * pi_hi / 2         # f(lo) > 0 > f(hi) is checked, not assumed
    flo, fhi = f_interval(lo), f_interval(hi)
    assert flo[0] > 0 and fhi[1] < 0, "f does not change sign on the starting bracket"
    for _ in range(120):
        m = (lo + hi) / 2
        m = F(round(m * 2 ** 160), 2 ** 160)
        fm = f_interval(m)
        if fm[0] > 0:
            lo = m
        elif fm[1] < 0:
            hi = m
        else:
            break
    print("x_0 (root of sin x = x cos x in (pi, 3 pi/2)) in [%s, %s], width %.1e" % (dec(lo, 30), dec(hi, 30), float(hi - lo)))
    th_lo = 1 / sqrt_hi(1 + hi * hi)
    th_hi = 1 / sqrt_lo(1 + lo * lo)
    if CONTROL == "max":
        th_lo, th_hi = F(-1), F(-1)        # theta_0 = -max sin(x)/x = -1: the bound becomes 1/(2(1 - 1)) undefined
    # the minimum value is below -1/(2 pi), so it is the global minimum over x > 0
    glob = -th_lo < -1 / (2 * pi_hi)
    print("min sin(x)/x = -theta_0 <= %s < -1/(2 pi): %s (so it is the global minimum)" % (dec(-th_lo, 12), glob))
    ok &= glob
    print("theta_0 in [%s, %s]" % (dec(th_lo, 40), dec(th_hi, 40)))
    if th_lo <= -1:
        print("bound undefined"); print("VERDICT REFUSED"); return
    b_lo = 1 / (2 * (1 + th_hi)); b_hi = 1 / (2 * (1 + th_lo))
    print("1/(2(1 + theta_0)) in [%s, %s]" % (dec(b_lo, 40), dec(b_hi, 40)))
    pt = "0.2172337" if CONTROL == "theta" else "0.2172336"
    pb = "0.4107674" if CONTROL == "bound" else "0.4107675"
    t1, r1 = trunc_or_round(th_lo, th_hi, pt)
    t2, r2 = trunc_or_round(b_lo, b_hi, pb)
    above = F(pb) >= b_hi
    print("page theta_0 = %s...: truncation %s; page bound %s...: truncation %s, correct rounding %s, at least the value"
          " (so valid as an upper bound) %s" % (pt, t1, pb, t2, r2, above))
    ok &= t1 and r2 and above
    print("key %s %s" % (dec(th_lo, 30), dec(b_lo, 30)))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
