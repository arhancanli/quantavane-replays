#!/usr/bin/env python3
"""7a, upper bound [ZZ2020] mu(pi) <= 7.103205334137..., program B (written separately from program A): the growth
constants through the saddle points, not through the cubic in N.

Zeilberger and Zudilin (arXiv 1912.06345, Proposition 2) give N_j = g(y_j), g(y) = y (y^2 + 6y + 25)^2/(y - 25)^3,
at the saddle points y_j. We take y_j as the stationary points of g: (log g)' = 1/y + 2(2y + 6)/(y^2 + 6y + 25)
- 3/(y - 25) = 0, which (multiplied out by hand and checked below at random rationals) is 2y^3 - 125y^2 - 500y - 625 = 0.
Its real root y_3 (66.33...) by bisection, the complex pair by deflation and the quadratic formula; then
N_3 = g(y_3) and |N_1| = |g(y_1)| in complex decimal arithmetic. Checks: the N_j satisfy the paper's cubic
108 N^3 - 2359989 N^2 + 138304 N - 2048 = 0 (to 10^-30 relative), and y_j match the paper's -1.91975076 -+ 1.01250889 i
and 66.33950152. Then mu = 1 + (log N_3 + K)/(-(log|N_1| + K)), K = -(5/2) log 2 + 4 - pi/(2 sqrt 3) + log(3 sqrt 3/4).
Arithmetic: decimal at 60 digits (ln, sqrt correctly rounded), pi by the Chudnovsky-free Gauss-Legendre iteration,
results moved outward by 10^-45. Red controls (--control): the claim 7.103205334136; the exponent 3 in (y - 25)^3
replaced by 2 in g (the stationary equation changes accordingly: y^2 terms re-derived, so N_3 changes); K with 2 log 2.
Standard library only.
"""
import sys
from decimal import Decimal as D, getcontext
from fractions import Fraction as F
import random

getcontext().prec = 60
PAD = D(10) ** -45
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
EXP = 2 if CONTROL == "g" else 3


def gauss_legendre_pi():
    a = D(1); b = 1 / D(2).sqrt(); t = D(1) / 4; p = D(1)
    for _ in range(8):
        an = (a + b) / 2; b = (a * b).sqrt(); t -= p * (a - an) ** 2; a = an; p *= 2
    return (a + b) ** 2 / (4 * t)


PI = gauss_legendre_pi()


def stationary_poly():
    """coefficients of y (y^2+6y+25)(y-25) [1/y + 2(2y+6)/(y^2+6y+25) - EXP/(y-25)], expanded (highest first)"""
    # (y^2+6y+25)(y-25) + 2(2y+6) y (y-25) - EXP y (y^2+6y+25)
    def mul(p, q):
        r = [0] * (len(p) + len(q) - 1)
        for i, x in enumerate(p):
            for j, y in enumerate(q):
                r[i + j] += x * y
        return r
    def add(*ps):
        n = max(len(p) for p in ps)
        out = [0] * n
        for p in ps:
            for i, x in enumerate(reversed(p)):
                out[n - 1 - i] += x
        return out
    q = [1, 6, 25]; l = [1, -25]; y = [1, 0]
    return add(mul(q, l), mul(mul([4, 12], y), l), [-EXP * c for c in mul(y, q)])


def main():
    P = stationary_poly()
    while P and P[0] == 0: P = P[1:]
    print("stationary polynomial of g (highest degree first): %s" % P)
    # check the expansion against the definition at random rationals
    rnd = random.Random(7)
    for _ in range(5):
        x = F(rnd.randint(-900, 900), rnd.randint(1, 97))
        if x == 0 or x == 25 or x * x + 6 * x + 25 == 0: continue
        lhs = x * (x * x + 6 * x + 25) * (x - 25) * (1 / x + 2 * (2 * x + 6) / (x * x + 6 * x + 25) - EXP / (x - 25))
        rhs = sum(c * x ** (len(P) - 1 - i) for i, c in enumerate(P))
        assert lhs == rhs, "expansion check"
    ev = lambda x: sum(D(c) * x ** (len(P) - 1 - i) for i, c in enumerate(P))
    lo, hi = D(30), D(200)
    assert ev(lo) * ev(hi) < 0
    for _ in range(200):
        m = (lo + hi) / 2
        if (ev(m) < 0) == (ev(lo) < 0): lo = m
        else: hi = m
    y3 = (lo + hi) / 2
    # deflate: P / (y - y3) = a y^2 + b y + c
    a = D(P[0]); b = D(P[1]) + a * y3; c = D(P[2]) + b * y3
    disc = b * b - 4 * a * c
    re = -b / (2 * a); im = (-disc).sqrt() / (2 * a) if disc < 0 else D(0)
    print("y_3 = %s; y_1,2 = %s -+ %s i (paper: 66.33950152, -1.91975076 -+ 1.01250889 i)"
          % (format(y3, ".20f"), format(re, ".20f"), format(im, ".20f")))

    def cmul(u, v): return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])
    def cdiv(u, v):
        n = v[0] * v[0] + v[1] * v[1]
        return ((u[0] * v[0] + u[1] * v[1]) / n, (u[1] * v[0] - u[0] * v[1]) / n)
    def g(y):
        q = (y[0] * y[0] - y[1] * y[1] + 6 * y[0] + 25, 2 * y[0] * y[1] + 6 * y[1])
        num = cmul(y, cmul(q, q))
        l = (y[0] - 25, y[1]); den = (D(1), D(0))
        for _ in range(EXP): den = cmul(den, l)
        return cdiv(num, den)
    N3 = g((y3, D(0)))[0]
    N1 = g((re, -im))
    absN1 = (N1[0] ** 2 + N1[1] ** 2).sqrt()
    cub = lambda n: ((108 * n[0] ** 3 - 3 * 108 * n[0] * n[1] ** 2) - 2359989 * (n[0] ** 2 - n[1] ** 2) + 138304 * n[0] - 2048,
                     (3 * 108 * n[0] ** 2 * n[1] - 108 * n[1] ** 3) - 2359989 * 2 * n[0] * n[1] + 138304 * n[1])
    r3 = cub((N3, D(0))); r1 = cub(N1)
    on_cubic = abs(r3[0]) < D(10) ** -30 * N3 ** 3 * 108 and abs(r1[0]) + abs(r1[1]) < D(10) ** -30
    print("N_3 = g(y_3) = %s; N_1 = g(y_1) = %s + %s i, |N_1| = %s" % (format(N3, ".20f"), format(N1[0], ".20f"),
          format(N1[1], ".20f"), format(absN1, ".20f")))
    print("N_1 and N_3 satisfy the paper's cubic (18): %s" % on_cubic)
    k2 = 2 if CONTROL == "log2" else D(5) / 2
    s3 = D(3).sqrt()
    K = -k2 * D(2).ln() + 4 - PI / (2 * s3) + (3 * s3 / 4).ln()
    sig = N3.ln() + K
    tau = -(absN1.ln() + K)
    mu = 1 + sig / tau
    mlo, mhi = mu - PAD, mu + PAD
    print("lim sup log|I'_n|/n = %s; lim log b'_n/n = %s" % (format(-tau, ".20f"), format(sig, ".20f")))
    print("1 + sigma/tau in [%s, %s]" % (format(mlo, ".30f"), format(mhi, ".30f")))
    printed = D("7.103205334136") if CONTROL == "claim" else D("7.103205334137")
    holds = mlo >= printed and mhi <= printed + D("1e-12")
    print("row %s: the leading digits of the value: %s" % (printed, holds))
    ok = on_cubic and holds
    print("key %s" % format(mlo, ".40f")[:27])
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
