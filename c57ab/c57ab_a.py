#!/usr/bin/env python3
"""57a and 57b, program A: the closed forms of [AG1937] (Bloch) and [Rad1943] (Landau) against the page's decimals.

  57a [AG1937]: B <= Gamma(1/3) Gamma(11/12) / (sqrt(1 + sqrt 3) Gamma(1/4)) ~ 0.4719 (Ahlfors-Grunsky).
  57b [Rad1943]: L <= Gamma(1/3) Gamma(5/6) / Gamma(1/6) ~ 0.5433 (Rademacher).
Closed forms named in the rows (the extremal-function analysis behind them is not checked): each value is enclosed
and compared with the page's decimal (rounding to the printed digits).

Route A: Gamma(x) = exp(log Gamma(x + 40)) / (x (x+1) ... (x+39)), log Gamma by Stirling's series with 20 Bernoulli
terms (for real z > 0 the remainder is bounded by the first omitted term); integer fixed point at 2^-200 with this
program's own exp, log, sqrt and pi. Self-tests: Gamma(1/2)^2 = pi, Gamma(1) = 1, and the duplication-reflection
identity Gamma(1/6) = 2^(-1/3) (3/pi)^(1/2) Gamma(1/3)^2. Standard library only.
Usage: c57ab_a.py [--print-ag D] [--print-rad D] [--variant ag-no-factor|rad-swapped]
"""
import sys
from fractions import Fraction as Q
from math import isqrt, comb

F = 200
ONE = 1 << F

def fl(a):
    return a >> F                         # floor of a / 2^F
def ce(a):
    return -((-a) >> F)

class FI:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        self.lo = lo; self.hi = lo if hi is None else hi
    @staticmethod
    def q(v):
        v = Q(v)
        return FI((v.numerator << F) // v.denominator, -((-(v.numerator << F)) // v.denominator))
    def __add__(self, o):
        if not isinstance(o, FI): o = FI.q(o)
        return FI(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__
    def __neg__(self):
        return FI(-self.hi, -self.lo)
    def __sub__(self, o):
        if not isinstance(o, FI): o = FI.q(o)
        return FI(self.lo - o.hi, self.hi - o.lo)
    def __rsub__(self, o):
        return FI.q(o) - self
    def __mul__(self, o):
        if isinstance(o, int):
            return FI(self.lo * o, self.hi * o) if o >= 0 else FI(self.hi * o, self.lo * o)
        if not isinstance(o, FI): o = FI.q(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return FI(fl(min(p)), ce(max(p)))
    __rmul__ = __mul__
    def recip(self):
        if self.lo > 0:
            return FI((ONE << F) // self.hi, -((-(ONE << F)) // self.lo))
        if self.hi < 0:
            return -((-self).recip())
        raise ArithmeticError("1/x with 0 in x")
    def __truediv__(self, o):
        if isinstance(o, int) and o > 0:
            return FI(self.lo // o, -((-self.hi) // o))
        if not isinstance(o, FI): o = FI.q(o)
        return self * o.recip()
    def hull(self, o):
        return FI(min(self.lo, o.lo), max(self.hi, o.hi))
    def mag(self):
        return max(abs(self.lo), abs(self.hi))
    def f(self):
        return self.lo / ONE

ZERO = FI(0)

# ---------------- elementary functions (this program's own) ----------------
def _exp_point(a):
    """exp(a / 2^F) as an interval"""
    s = max(0, abs(a).bit_length() - F + 12)
    y = FI(a >> s, -((-a) >> s))           # |y| < 2^-11
    tot = FI(ONE); term = FI(ONE)
    for k in range(1, 26):
        term = term * y / k
        tot = tot + term
    r = 2 * (y.mag() ** 26 // (ONE ** 25)) // 1 + 1   # 2|y|^26/26! < 2|y|^26 (scaled, rounded up)
    tot = tot + FI(-r, r)
    for _ in range(s):
        tot = tot * tot
    return tot
def exp_fi(x):
    return FI(_exp_point(x.lo).lo, _exp_point(x.hi).hi)

def _atanh_series(u):
    """2 atanh(u) for an interval 0 <= u <= 1/3"""
    tot = ZERO; u2 = u * u; p = u; k = 0
    while True:
        tot = tot + p / (2 * k + 1)
        p = p * u2; k += 1
        if p.hi * 9 < 1 << 10:             # remaining tail <= p/(1 - u^2) < 2 p
            tot = tot + FI(0, 2 * p.hi + 2)
            return tot * 2
LOG2 = _atanh_series(FI.q(Q(1, 3)))
def _log_point(a):
    assert a > 0
    e = a.bit_length() - 1 - F            # a/2^F = 2^e y, 1 <= y < 2
    yl = a << (-e) if e < 0 else a >> e
    yh = yl if e <= 0 else -((-a) >> e)
    y = FI(yl, yh)
    u = (y - ONE * 0 - FI(ONE)) / (y + FI(ONE))
    u = FI(max(u.lo, 0), u.hi)
    return LOG2 * e + _atanh_series(u)
def log_fi(x):
    return FI(_log_point(x.lo).lo, _log_point(x.hi).hi)

def sqrt_fi(x):
    lo = isqrt(max(x.lo, 0) << F)
    hi = isqrt(x.hi << F)
    if hi * hi < x.hi << F: hi += 1
    return FI(lo, hi)

def _arctan_inv(m):
    tot = ZERO; k = 0; mm = m * m; den = m
    while True:
        t = FI.q(Q(1, (2 * k + 1) * den))
        if t.hi < 4:
            return tot + FI(-4, 4)
        tot = tot + (t if k % 2 == 0 else -t)
        den *= mm; k += 1
PI = _arctan_inv(5) * 16 - _arctan_inv(239) * 4
SQRT2 = sqrt_fi(FI(2 * ONE))
INV_SQRT_2PI = sqrt_fi(PI * 2).recip()
LOG_SQRT_2PI = log_fi(PI * 2) / 2


def bernoulli(n):
    B = [Q(0)] * (n + 1); B[0] = Q(1)
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, j) * B[j] for j in range(m)) / (m + 1)
    return B
KST = 20
BER = bernoulli(2 * KST + 2)
LOG_2PI = log_fi(PI * 2)

def gamma(x):
    """Gamma(x) for a rational 0 < x <= 1"""
    x = Q(x); z = x + 40
    Z = FI.q(z)
    lg = (Z - FI.q(Q(1, 2))) * log_fi(Z) - Z + LOG_2PI / 2
    for k in range(1, KST + 1):
        lg = lg + FI.q(BER[2 * k] / ((2 * k) * (2 * k - 1) * z ** (2 * k - 1)))
    r = FI.q(abs(BER[2 * KST + 2]) / ((2 * KST + 2) * (2 * KST + 1) * z ** (2 * KST + 1)))
    lg = lg + FI(-r.hi, r.hi)
    prod = Q(1)
    for j in range(40):
        prod *= x + j
    e = exp_fi(lg)                       # divide by the exact rational prod = n/d: (e * d) / n with outward rounding
    n, d = prod.numerator, prod.denominator
    return FI((e.lo * d) // n, -((-(e.hi * d)) // n))

def rounds_to(v, printed):
    d = len(printed.split(".")[1])
    p = FI.q(Q(printed)); h = FI.q(Q(5, 10 ** (d + 1)))
    return (p - h).hi <= v.lo and v.hi <= (p + h).lo

def show(x, n=30):
    s = str((x.lo * 10 ** n) >> F)
    return (s[:-n] or "0") + "." + s[-n:]

def main():
    a = sys.argv[1:]
    opt = dict(zip(a[0::2], a[1::2]))
    variant = opt.get("--variant")
    pag = opt.get("--print-ag", "0.4719"); prad = opt.get("--print-rad", "0.5433")
    if variant:
        print("CONTROL: variant", variant)
    print("57a/57b program A (Stirling)")
    g12 = gamma(Q(1, 2)); g1 = gamma(1)
    G = {q: gamma(q) for q in (Q(1, 3), Q(1, 4), Q(11, 12), Q(5, 6), Q(1, 6))}
    t1 = abs((g12 * g12 - PI).f()) < 1e-45 and abs((g1 - FI(ONE)).f()) < 1e-45
    dup = exp_fi(-(LOG2 / 3)) * sqrt_fi(PI.recip() * 3) * G[Q(1, 3)] * G[Q(1, 3)]
    t2 = abs((dup - G[Q(1, 6)]).f()) < 1e-45
    print("self-test Gamma(1/2)^2 = pi, Gamma(1) = 1: %s; Gamma(1/6) = 2^(-1/3) (3/pi)^(1/2) Gamma(1/3)^2: %s" % (t1, t2))
    if not (t1 and t2):
        print("VERDICT FAILED (self-test)"); return
    for q in sorted(G):
        print("Gamma(%s) = %s..." % (q, show(G[q])))
    s3 = sqrt_fi(FI(3 * ONE))
    fac = sqrt_fi(s3 + FI(ONE)).recip() if variant != "ag-no-factor" else FI(ONE)
    ag = G[Q(1, 3)] * G[Q(11, 12)] / G[Q(1, 4)] * fac
    rad = G[Q(1, 3)] * G[Q(5, 6)] / G[Q(1, 6)] if variant != "rad-swapped" else G[Q(1, 3)] * G[Q(1, 6)] / G[Q(5, 6)]
    r1 = rounds_to(ag, pag); r2 = rounds_to(rad, prad)
    print("row 57a [AG1937] Gamma(1/3) Gamma(11/12)/(sqrt(1 + sqrt 3) Gamma(1/4)) = %s...: ~ %s: %s" % (show(ag), pag, r1))
    print("row 57b [Rad1943] Gamma(1/3) Gamma(5/6)/Gamma(1/6) = %s...: ~ %s: %s" % (show(rad), prad, r2))
    print("VERDICT PASS" if (r1 and r2) else "VERDICT REFUSED: a printed value does not match")

if __name__ == "__main__":
    main()
