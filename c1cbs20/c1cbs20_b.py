#!/usr/bin/env python3
"""1c, upper bound [BS2020]: C_1c <= 1/(2(1 + theta_0)), theta_0 = -min_x sin(x)/x.  Program B (written separately
from program A): the minimum by branch and bound, with no use of the critical-point equation.

g(x) = sin(x)/x is even. On (0, pi/2], g >= 2/pi > 0 (Jordan's inequality); for x >= 4 pi, g >= -1/(4 pi) > -0.08.
On [pi/2, 4 pi] we bound g on cells [m - h, m + h] by the second-order Taylor bound
    g(m + u) >= g(m) - |g'(m)| h - M2 h^2 / 2,   M2 = 1/a + 2/a^2 + 2/a^3 >= sup |g''| on the cell (a = m - h),
from g'' = -sin x/x - 2 cos x/x^2 + 2 sin x/x^3. Cells whose bound is above the best value found so far are dropped,
the others are bisected, until every cell's lower bound is within 10^-32 of the best value g(x^) (an upper bound for
the minimum). So min g is enclosed in [L, g(x^)] with g(x^) - L <= 10^-32, theta_0 = -min g, and the bound
1/(2(1 + theta_0)) follows. Arithmetic: Python decimal at 80 digits; every sin, cos and quotient is computed by our own
series (Taylor after reduction by 2 pi, pi from Machin's formula) and moved outward by 10^-70, far above the
accumulated rounding (under 10^-75). Red controls (--control): the page's theta_0 claimed as 0.2172337...; the bound
claimed as 0.4107674; the search restricted to [pi/2, pi] (where g >= 0, so the minimum is missed). Standard library only.
"""
import sys
from decimal import Decimal as D, getcontext

getcontext().prec = 80
PAD = D(10) ** -70
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None


def atan_inv(n):
    x = D(1) / n; s = D(0); p = x; k = 0
    while p > D(10) ** -85:
        s += (p / (2 * k + 1)) if k % 2 == 0 else -(p / (2 * k + 1))
        p = p * x * x; k += 1
    return s


PI = 16 * atan_inv(5) - 4 * atan_inv(239)


def sin_cos(x):
    k = (x / (2 * PI)).to_integral_value()
    y = x - 2 * PI * k
    s = D(0); c = D(0); t = D(1); n = 0
    while True:
        if n > 0:
            t = t * y / n
        r = n % 4
        if r == 0: c += t
        elif r == 1: s += t
        elif r == 2: c -= t
        else: s -= t
        n += 1
        if n > 10 and abs(t) < D(10) ** -85:
            break
    return s, c


def g_and_d(m):
    s, c = sin_cos(m)
    g = s / m
    d1 = c / m - s / (m * m)
    return g, d1


def main():
    assert D("3.14159265358979323846264338327950288") < PI < D("3.14159265358979323846264338327950289")
    lo0, hi0 = PI / 2, 4 * PI
    if CONTROL == "range":
        hi0 = PI
    # start: 3000 cells; best value from the cell centres
    cells = []
    n0 = 3000
    w = (hi0 - lo0) / n0
    best = D(10)
    for i in range(n0):
        a = lo0 + w * i
        cells.append((a, a + w))
    tol = D(10) ** -32
    xbest = None
    rounds = 0
    while True:
        rounds += 1
        nxt = []
        bounds = []
        for a, b in cells:
            m = (a + b) / 2; h = (b - a) / 2
            g, d1 = g_and_d(m)
            if g + PAD < best:
                best = g + PAD; xbest = m
            M2 = 1 / a + 2 / (a * a) + 2 / (a * a * a)
            lb = g - abs(d1) * h - M2 * h * h / 2 - PAD
            bounds.append((lb, a, b))
        for lb, a, b in bounds:
            if lb > best:
                continue                      # this cell cannot hold the minimum
            if best - lb <= tol:
                nxt.append((a, b, lb)); continue
            m = (a + b) / 2
            nxt.append((a, m, None)); nxt.append((m, b, None))
        if all(x[2] is not None for x in nxt):
            L = min(x[2] for x in nxt)
            break
        cells = [(a, b) for a, b, _ in nxt]
        assert rounds < 400 and len(cells) < 200000
    U = best
    print("min of sin(x)/x on [pi/2, %s]: in [%s, %s] (x near %s; %d final cells, %d rounds)"
          % ("pi" if CONTROL == "range" else "4 pi", format(L, ".36f"), format(U, ".36f"), format(xbest, ".20f"), len(nxt), rounds))
    # outside: g >= 2/pi > 0 on (0, pi/2]; g >= -1/(4 pi) for x >= 4 pi
    outside = -1 / (4 * PI)
    glob = U < outside
    print("below -1/(4 pi) = %s: %s, so this is the minimum over all x" % (str(outside)[:12], glob))
    th_lo, th_hi = -U, -L
    b_lo = 1 / (2 * (1 + th_hi)) - PAD; b_hi = 1 / (2 * (1 + th_lo)) + PAD
    print("theta_0 in [%s, %s]" % (format(th_lo, ".34f"), format(th_hi, ".34f")))
    print("1/(2(1 + theta_0)) in [%s, %s]" % (format(b_lo, ".34f"), format(b_hi, ".34f")))
    pt = D("0.2172337") if CONTROL == "theta" else D("0.2172336")
    pb = D("0.4107674") if CONTROL == "bound" else D("0.4107675")
    u = D("0.0000001")
    t1 = pt <= th_lo and th_hi < pt + u
    r2 = abs(b_lo - pb) <= u / 2 and abs(b_hi - pb) <= u / 2
    above = pb >= b_hi
    print("page theta_0 = %s...: truncation %s; page bound %s: correct rounding %s, at least the value %s"
          % (pt, t1, pb, r2, above))
    ok = glob and t1 and r2 and above
    print("key %s %s" % (format(th_lo, ".34f")[:32], format(b_lo, ".34f")[:32]))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
