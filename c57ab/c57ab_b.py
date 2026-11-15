#!/usr/bin/env python3
"""57a and 57b, program B (written separately from program A): the same two closed forms (see program A).

Route B: for 0 < s <= 1 and N = 120,
    Gamma(s) = N^s e^(-N) sum_{k>=0} N^k / (s (s+1) ... (s+k))  +  E,   0 <= E <= N^(s-1) e^(-N)
(the lower incomplete gamma function as a positive series, and the upper one bounded by t^(s-1) <= N^(s-1) on
[N, oo)); the series is cut once the ratio N/(s+k+1) is below 1/2, the tail bounded by the last term. Decimal
intervals at 110 digits; exp and ln of the decimal module, correctly rounded and widened by one unit in the last
place. Self-tests: Gamma(1/2)^2 = pi, and the reflection Gamma(1/12) Gamma(11/12) = pi / sin(pi/12) with
sin(pi/12) = (sqrt 6 - sqrt 2)/4. Standard library only.
Usage: c57ab_b.py [--print-ag D] [--print-rad D] [--variant ag-no-factor|rad-swapped]
"""
import sys
from fractions import Fraction as Q
from decimal import Decimal as D, localcontext, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

DPREC = 110
def _ctx(rnd, prec=DPREC):
    return localcontext(Context(prec=prec, rounding=rnd))

class DI:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        self.lo = lo; self.hi = lo if hi is None else hi
    @staticmethod
    def q(v):
        v = Q(v)
        with _ctx(ROUND_FLOOR): lo = D(v.numerator) / D(v.denominator)
        with _ctx(ROUND_CEILING): hi = D(v.numerator) / D(v.denominator)
        return DI(lo, hi)
    def __add__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        with _ctx(ROUND_FLOOR): lo = self.lo + o.lo
        with _ctx(ROUND_CEILING): hi = self.hi + o.hi
        return DI(lo, hi)
    __radd__ = __add__
    def __neg__(self):
        return DI(self.hi.copy_negate(), self.lo.copy_negate())
    def __sub__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        return self + (-o)
    def __rsub__(self, o):
        return DI.q(o) - self
    def __mul__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        a, b, c, d = self.lo, self.hi, o.lo, o.hi
        with _ctx(ROUND_FLOOR): lo = min(a * c, a * d, b * c, b * d)
        with _ctx(ROUND_CEILING): hi = max(a * c, a * d, b * c, b * d)
        return DI(lo, hi)
    __rmul__ = __mul__
    def recip(self):
        if not (self.lo > 0 or self.hi < 0):
            raise ArithmeticError("division by an interval containing 0")
        with _ctx(ROUND_FLOOR): lo = D(1) / self.hi
        with _ctx(ROUND_CEILING): hi = D(1) / self.lo
        return DI(lo, hi)
    def __truediv__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        return self * o.recip()
    def mag(self):
        return max(self.lo.copy_abs(), self.hi.copy_abs())
    def hull(self, o):
        return DI(min(self.lo, o.lo), max(self.hi, o.hi))
    def mid(self):
        return (self.lo + self.hi) / 2

def _cr(f, x):
    """a correctly rounded (half-even) decimal function, widened by one ulp each way"""
    with _ctx(ROUND_HALF_EVEN) as c:
        v = f(x, c)
        return v.next_minus(c), v.next_plus(c)

def exp_di(x):
    lo, _ = _cr(lambda v, c: v.exp(c), x.lo)
    _, hi = _cr(lambda v, c: v.exp(c), x.hi)
    return DI(lo, hi)
def ln_di(x):
    assert x.lo > 0
    lo, _ = _cr(lambda v, c: v.ln(c), x.lo)
    _, hi = _cr(lambda v, c: v.ln(c), x.hi)
    return DI(lo, hi)
def sqrt_di(x):
    assert x.lo >= 0
    lo, _ = _cr(lambda v, c: v.sqrt(c), x.lo)
    _, hi = _cr(lambda v, c: v.sqrt(c), x.hi)
    return DI(max(lo, D(0)), hi)

def machin_pi_exact(digits):
    """a rational within 10^-digits of pi (Machin, alternating series cut where the term is below 10^-(digits+5))"""
    def at(m):
        s = Q(0); k = 0
        while True:
            t = Q(1, (2 * k + 1) * m ** (2 * k + 1))
            if t < Q(1, 10 ** (digits + 5)):
                return s
            s += t if k % 2 == 0 else -t; k += 1
    return 16 * at(5) - 4 * at(239)
_pi = machin_pi_exact(DPREC + 5)
PI = DI.q(_pi - Q(1, 10 ** (DPREC + 3))).hull(DI.q(_pi + Q(1, 10 ** (DPREC + 3))))

NG = 120

def gamma(s):
    """Gamma(s) for a rational 0 < s <= 1"""
    s = Q(s)
    S = DI.q(s)
    term = S.recip()                     # k = 0: 1/s
    tot = term
    k = 0
    while True:
        k += 1
        term = term * NG / (S + k)
        tot = tot + term
        if k > 2 * NG and term.hi < D(10) ** -100:
            tot = tot + DI(D(0), term.hi)          # remaining ratio < 1/2: tail <= the last term
            break
    lnN = ln_di(DI.q(NG))
    pref = exp_di(S * lnN - NG)          # N^s e^(-N)
    err = exp_di((S - 1) * lnN - NG)     # N^(s-1) e^(-N)
    return pref * tot + DI(D(0), err.hi)

def rounds_to(v, printed):
    d = len(printed.split(".")[1])
    p = D(printed); h = D(5) / D(10) ** (d + 1)
    with localcontext(Context(prec=80)):
        return p - h <= v.lo and v.hi <= p + h

def main():
    a = sys.argv[1:]
    opt = dict(zip(a[0::2], a[1::2]))
    variant = opt.get("--variant")
    pag = opt.get("--print-ag", "0.4719"); prad = opt.get("--print-rad", "0.5433")
    if variant:
        print("CONTROL: variant", variant)
    print("57a/57b program B (incomplete gamma series)")
    g12 = gamma(Q(1, 2))
    G = {q: gamma(q) for q in (Q(1, 3), Q(1, 4), Q(11, 12), Q(5, 6), Q(1, 6), Q(1, 12))}
    s2, s6 = sqrt_di(DI.q(2)), sqrt_di(DI.q(6))
    refl = G[Q(1, 12)] * G[Q(11, 12)] - PI / ((s6 - s2) / 4)
    sq = g12 * g12 - PI
    t1 = sq.lo.copy_abs() < D(10) ** -45 and sq.hi.copy_abs() < D(10) ** -45
    t2 = refl.lo.copy_abs() < D(10) ** -45 and refl.hi.copy_abs() < D(10) ** -45
    print("self-test Gamma(1/2)^2 = pi: %s; Gamma(1/12) Gamma(11/12) = pi/sin(pi/12): %s" % (t1, t2))
    if not (t1 and t2):
        print("VERDICT FAILED (self-test)"); return
    for q in sorted(G):
        print("Gamma(%s) = %s..." % (q, str(G[q].lo)[:32]))
    s3 = sqrt_di(DI.q(3))
    fac = sqrt_di(s3 + 1).recip() if variant != "ag-no-factor" else DI.q(1)
    ag = G[Q(1, 3)] * G[Q(11, 12)] / G[Q(1, 4)] * fac
    rad = G[Q(1, 3)] * G[Q(5, 6)] / G[Q(1, 6)] if variant != "rad-swapped" else G[Q(1, 3)] * G[Q(1, 6)] / G[Q(5, 6)]
    r1 = rounds_to(ag, pag); r2 = rounds_to(rad, prad)
    print("row 57a [AG1937] = %s...: ~ %s: %s" % (str(ag.lo)[:32], pag, r1))
    print("row 57b [Rad1943] = %s...: ~ %s: %s" % (str(rad.lo)[:32], prad, r2))
    print("VERDICT PASS" if (r1 and r2) else "VERDICT REFUSED: a printed value does not match")

if __name__ == "__main__":
    main()
