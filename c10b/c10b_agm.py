#!/usr/bin/env python3
"""Program 2 for 10b (the complex Grothendieck constant): Haagerup's upper bound [H1987] by a second route that shares
no code with c10b_series.py: complete elliptic integrals by the arithmetic-geometric mean, in interval arithmetic on
Python's decimal module (every operation rounded outward: floor for lower ends, ceiling for upper ends).

With m = k^2, a_0 = 1, b_0 = sqrt(1 - m), c_0^2 = m, a_{n+1} = (a_n + b_n)/2, b_{n+1} = sqrt(a_n b_n),
c_{n+1} = (a_n - b_n)/2 and M = AGM(1, sqrt(1 - m)):
    K(m) = pi / (2M),   E(m) = K(m) (1 - sum_{n >= 0} 2^(n-1) c_n^2)        (Gauss, Legendre)
and int_0^{pi/2} cos^2 t (1 - m sin^2 t)^(-1/2) dt = (E - (1 - m) K) / m. Haagerup's equation
k int_0^{pi/2} cos^2 t (1 - k^2 sin^2 t)^(-1/2) dt = (pi/8)(k + 1) becomes, after pi cancels,
    F(k) = 4 (m - sum_n 2^(n-1) c_n^2) / (M k) - (k + 1) = 0,
and the bound is 8 / (pi (K_0 + 1)). After n >= 1, b_n <= M <= a_n; once c_{N+1} < 10^-30 the rest of the sum is at
most 2^(N+1) c_{N+1}^2 (each c_{n+1} = c_n^2 / (4 a_{n+1}) with a_{n+1} > 1/2). F is increasing; K_0 is bracketed by
bisection on F's interval. pi by Salamin-Brent: pi = 4 M^2 / (1 - sum_{n >= 1} 2^(n+1) c_n^2) for b_0 = 1/sqrt(2).

Also (information): Haagerup's conjectured value (int cos^2 t (1 + sin^2 t)^(-1/2) dt)^(-1) = (2K(-1) - E(-1))^(-1)
by the same formulas at m = -1 (AGM(1, sqrt 2)), and Pisier's e^(1 - gamma), gamma by Euler-Maclaurin at n = 10^5
(decimal's ln and exp are correctly rounded, so +-1 unit in the last place encloses them).

  c10b_agm.py [--claim-h DEC] [--claim-p DEC]
Standard library only; run with python3 -I.
"""
import sys
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING

PREC = 70
DN = Context(prec=PREC, rounding=ROUND_FLOOR)
UP = Context(prec=PREC, rounding=ROUND_CEILING)


class I:
    """A closed interval [lo, hi] of decimals."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = D(lo)
        self.hi = D(lo) if hi is None else D(hi)
        assert self.lo <= self.hi

    def __add__(s, o):
        o = o if isinstance(o, I) else I(o)
        return I(DN.add(s.lo, o.lo), UP.add(s.hi, o.hi))

    def __sub__(s, o):
        o = o if isinstance(o, I) else I(o)
        return I(DN.subtract(s.lo, o.hi), UP.subtract(s.hi, o.lo))

    def __mul__(s, o):
        o = o if isinstance(o, I) else I(o)
        ps = [(x, y) for x in (s.lo, s.hi) for y in (o.lo, o.hi)]
        return I(min(DN.multiply(x, y) for x, y in ps), max(UP.multiply(x, y) for x, y in ps))

    def __truediv__(s, o):
        o = o if isinstance(o, I) else I(o)
        assert o.lo > 0 or o.hi < 0, "division by an interval containing 0"
        ps = [(x, y) for x in (s.lo, s.hi) for y in (o.lo, o.hi)]
        return I(min(DN.divide(x, y) for x, y in ps), max(UP.divide(x, y) for x, y in ps))

    def sqrt(s):
        assert s.lo >= 0
        return I(DN.sqrt(s.lo), UP.sqrt(s.hi))

    def sq(s):
        if s.lo >= 0:
            return I(DN.multiply(s.lo, s.lo), UP.multiply(s.hi, s.hi))
        if s.hi <= 0:
            return I(DN.multiply(s.hi, s.hi), UP.multiply(s.lo, s.lo))
        m = max(s.lo.copy_negate(), s.hi)
        return I(0, UP.multiply(m, m))

    def hull(s, o):
        return I(min(s.lo, o.lo), max(s.hi, o.hi))


TINY = D("1e-30")


def agm_sums(a, b, c0sq):
    """Returns (M, S) with M = AGM(a, b) and S = sum_{n >= 0} 2^(n-1) c_n^2 (c_0^2 given), as intervals."""
    S = c0sq * I("0.5")
    n = 0
    while True:
        c = (a - b) * I("0.5")
        a, b = (a + b) * I("0.5"), (a * b).sqrt()
        n += 1
        S = S + c.sq() * I(D(2) ** (n - 1))
        cnext = (a - b) * I("0.5")
        if max(cnext.lo.copy_abs(), cnext.hi.copy_abs()) < TINY:
            cmax = max(cnext.lo.copy_abs(), cnext.hi.copy_abs())
            tail = I(0, UP.multiply(UP.multiply(D(2) ** (n + 1), cmax), cmax))
            S = S + tail
            lo, hi = min(a.lo, b.lo), max(a.hi, b.hi)     # the AGM lies between b_n and a_n
            return I(lo, hi), S


def pi_interval():
    a, b = I(1), I(1) / I(2).sqrt()
    S = I(0)
    n = 0
    while True:
        c = (a - b) * I("0.5")
        a, b = (a + b) * I("0.5"), (a * b).sqrt()
        n += 1
        S = S + c.sq() * I(D(2) ** (n + 1))
        cnext = (a - b) * I("0.5")
        if max(cnext.lo.copy_abs(), cnext.hi.copy_abs()) < TINY:
            cmax = max(cnext.lo.copy_abs(), cnext.hi.copy_abs())
            S = S + I(0, UP.multiply(UP.multiply(D(2) ** (n + 3), cmax), cmax))
            M = I(min(a.lo, b.lo), max(a.hi, b.hi))
            return I(4) * M.sq() / (I(1) - S)


def F(k):
    m = k.sq()
    M, S = agm_sums(I(1), (I(1) - m).sqrt(), m)
    return I(4) * (m - S) / (M * k) - (k + I(1))


def gamma_interval(n=100000):
    lo = hi = D(0)
    for i in range(1, n + 1):
        lo = DN.add(lo, DN.divide(1, i))
        hi = UP.add(hi, UP.divide(1, i))
    ln = Context(prec=PREC).ln(D(n))
    ulp = D(10) ** (ln.adjusted() - PREC + 1)
    lnn = I(DN.subtract(ln, ulp), UP.add(ln, ulp))
    nn = D(n)
    corr = (I(1) / I(2 * nn) * I(-1) + I(1) / I(12 * nn * nn)
            - I(1) / I(120 * nn ** 4) + I(1) / I(252 * nn ** 6))
    g = I(lo, hi) - lnn + corr
    return I(DN.subtract(g.lo, UP.divide(1, 240 * nn ** 8)), g.hi)


def exp_interval(x):
    lo = Context(prec=PREC).exp(x.lo)
    hi = Context(prec=PREC).exp(x.hi)
    return I(DN.subtract(lo, D(10) ** (lo.adjusted() - PREC + 1)), UP.add(hi, D(10) ** (hi.adjusted() - PREC + 1)))


def fmt(x, d=30):
    q = D(10) ** -d
    big_dn = Context(prec=200, rounding=ROUND_FLOOR)
    big_up = Context(prec=200, rounding=ROUND_CEILING)
    return "[%s, %s]" % (big_dn.quantize(x.lo, q), big_up.quantize(x.hi, q))


def main():
    claim_h = claim_p = None
    a = sys.argv[1:]
    while a:
        if a[0] == "--claim-h":
            claim_h = D(a[1]); a = a[2:]
        elif a[0] == "--claim-p":
            claim_p = D(a[1]); a = a[2:]
        else:
            raise SystemExit("unknown option " + a[0])
    # anchors of the elliptic routine: K(0) = E(0) = pi/2; at m = 1/2, K(1/2)^2 ... checked through Legendre's relation
    PI = pi_interval()
    print("pi in %s (Salamin-Brent)" % fmt(PI, 40))
    M, S = agm_sums(I(1), I("0.5").sqrt(), I("0.5"))
    Kh = PI / (I(2) * M)
    Eh = Kh * (I(1) - S)
    leg = I(2) * Kh * Eh - Kh.sq()              # Legendre at m = 1/2: 2 K E - K^2 = pi/2
    print("Legendre's relation at m = 1/2: 2KE - K^2 - pi/2 in %s" % fmt(leg - PI * I("0.5"), 45))
    lo, hi = D("0.5"), D("0.99")
    assert F(I(lo)).hi < 0 < F(I(hi)).lo
    while hi - lo > D("1e-36"):
        mid = Context(prec=200, rounding=ROUND_FLOOR).quantize(DN.divide(DN.add(lo, hi), 2), D("1e-40"))
        v = F(I(mid))
        if v.lo > 0:
            hi = mid
        elif v.hi < 0:
            lo = mid
        else:
            break
    assert F(I(lo)).hi < 0 < F(I(hi)).lo
    K0 = I(lo, hi)
    print("K_0 in %s" % fmt(K0))
    H = I(8) / (PI * (K0 + I(1)))
    print("Haagerup's bound 8/(pi(K_0+1)) in %s" % fmt(H))
    ok = True
    if claim_h is not None:
        holds = H.hi <= claim_h
        print("printed upper bound %s [H1987]: %s (printed minus value >= %.3e)"
              % (claim_h, "holds" if holds else "FAILS", float(claim_h - H.hi)))
        ok &= holds
    # information: Haagerup's conjecture, m = -1: int cos^2 t (1 + sin^2 t)^(-1/2) dt = 2K(-1) - E(-1)
    M1, S1 = agm_sums(I(1), I(2).sqrt(), I(-1))
    K1 = PI / (I(2) * M1)
    E1 = K1 * (I(1) - S1)
    conj = I(1) / (I(2) * K1 - E1)
    print("information: Haagerup's conjectured value (int cos^2 t (1 + sin^2 t)^(-1/2) dt)^(-1) in %s" % fmt(conj))
    g = gamma_interval()
    print("information: Euler's gamma in %s" % fmt(g))
    e = exp_interval(I(1) - g)
    print("information: e^(1 - gamma) [P1978] in %s" % fmt(e))
    if claim_p is not None:
        holds = e.hi <= claim_p
        print("printed upper bound %s [P1978]: %s" % (claim_p, "holds" if holds else "FAILS"))
        ok &= holds
    print("VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
