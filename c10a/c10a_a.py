#!/usr/bin/env python3
"""10a, program A: values behind five rows of the real Grothendieck constant page, from their definitions.

  * [Dav1984], [Ree1991] K_DR ~ 1.67696...: the Davie-Reeds bound as defined in Jones and Malavolta (arXiv 2603.30039,
    Lemma unique-c and the Davie-Reeds lemma): C* is the unique root in (0, 1) of H(C) = 4 phi(C)^2 - 4 Phi(-C) + 1
    (H' = 4 phi (1 - 2 C phi) > 0 there, since 2 C phi(C) < 2 phi(0) < 1), lambda* = 2 C* phi(C*), and
    K_DR = 1/R(C*), R(C) = (4 phi^2 - 2 C phi (4 Phi(-C) - 1))/(1 - 2 C phi). At the root 4 Phi(-C*) - 1 = 4 phi(C*)^2,
    so R(C*) = 4 phi(C*)^2 and K_DR = (pi/2) exp(C*^2); this program uses that identity (program B does not).
  * [G1953] sinh(pi/2) ~ 2.30130 (upper) and pi/2 ~ 1.57080 (lower); [K1979] pi/(2 log(1 + sqrt 2)) ~ 1.782214;
    [SLXCKKM26] 6 pi/11 ~ 1.71360: closed forms, checked against the page's decimals (rounding to the printed digits).
Arithmetic: integer fixed point at 2^-200, interval ends rounded outward; exp, log, sqrt, pi and the erf series are
this program's own. A value from its definition for K_DR; the other four are closed forms. Standard library only.
Usage: c10a_a.py [--print-kdr D] [--print-sinh D] [--variant phi-plus]
"""
import sys
from fractions import Fraction as Q
from math import isqrt

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

def phi_fi(x):
    x2 = x * x
    if x.lo < 0 < x.hi:
        x2 = FI(0, x2.hi)
    return exp_fi(-x2 / 2) * INV_SQRT_2PI

def _Psi_series(y):
    """Psi(y) = (1 - erf(y/sqrt 2))/2 for a point 0 <= y <= 4: the alternating erf series"""
    x = FI(y) / SQRT2
    x2 = x * x
    term = x; tot = ZERO; n = 0
    while True:
        tot = tot + term / (2 * n + 1)
        nxt = -(term * x2) / (n + 1)
        n += 1
        if n > x2.hi // ONE + 2 and nxt.mag() < 1 << 8:
            tot = tot + FI(-nxt.mag(), nxt.mag())          # alternating, decreasing: |tail| <= first omitted term
            break
        term = nxt
    erf = tot * 2 * sqrt_fi(PI).recip()
    return (FI(ONE) - erf) / 2


def phi_fi(x):
    x2 = x * x
    if x.lo < 0 < x.hi:
        x2 = FI(0, x2.hi)
    return exp_fi(-x2 / 2) * INV_SQRT_2PI

def Psi(c):
    """Psi(c) = P(Z > c) for an exact point 0 <= c <= 4 (scaled integer)"""
    return _Psi_series(c)

def H_of(c, variant):
    """H(C) = 4 phi(C)^2 - 4 Phi(-C) + 1 at an exact point (Phi(-C) = Psi(C)); control: Phi(C) in its place"""
    ph = phi_fi(FI(c))
    Pm = Psi(c) if variant != "phi-plus" else FI(ONE) - Psi(c)
    return ph * ph * 4 - Pm * 4 + FI(ONE)

def rounds_to(val, printed):
    """the interval val lies within half a unit of the printed decimal's last digit"""
    d = len(printed.split(".")[1]) if "." in printed else 0
    p = FI.q(Q(printed)); h = FI.q(Q(5, 10 ** (d + 1)))
    return (p - h).hi <= val.lo and val.hi <= (p + h).lo

def show(x, n=30):
    """x as a decimal string (the lower end, n digits)"""
    s = str((x.lo * 10 ** n) >> F)
    return (s[:-n] or "0") + "." + s[-n:]

def main():
    a = sys.argv[1:]
    opt = dict(zip(a[0::2], a[1::2]))
    variant = opt.get("--variant")
    pk = opt.get("--print-kdr", "1.67696"); ps = opt.get("--print-sinh", "2.30130")
    if variant:
        print("CONTROL: H with Phi(C) in place of Phi(-C)")
    print("10a program A")
    # self-tests
    t1 = abs(PI.f() - 3.141592653589793) < 1e-15
    p1 = Psi(ONE); t2 = abs(p1.lo / ONE - 0.15865525393145705) < 1e-15
    t3 = abs(_exp_point(ONE).lo / ONE - 2.718281828459045) < 1e-15 and abs(LOG2.f() - 0.6931471805599453) < 1e-15
    print("self-test pi, Psi(1), e, log 2: %s" % (t1 and t2 and t3))
    if not (t1 and t2 and t3):
        print("VERDICT FAILED (self-test)"); return
    ok = True
    # C*: bisection on [0, 1] with exact dyadic points
    lo, hi = 0, ONE
    hl, hh = H_of(lo, variant), H_of(hi, variant)
    print("H(0) in [%.12f, %.12f]; H(1) in [%.12f, %.12f]" % (hl.f(), hl.hi / ONE, hh.f(), hh.hi / ONE))
    if not (hl.hi < 0 < hh.lo):
        print("no sign change of H on [0, 1]")
        print("VERDICT REFUSED: C* not found"); return
    for _ in range(150):
        mid = (lo + hi) // 2
        hm = H_of(mid, variant)
        if hm.hi < 0: lo = mid
        elif hm.lo > 0: hi = mid
        else: break
    C = FI(lo, hi)
    print("C* in [%s, %s]" % (show(FI(lo)), show(FI(hi))))
    lam = C * phi_fi(C) * 2
    kdr = PI / 2 * exp_fi(C * C)
    print("information: lambda* = 2 C* phi(C*) = %s" % show(lam, 12))
    print("K_DR = (pi/2) exp(C*^2) in [%s, %s]" % (show(kdr), show(FI(kdr.hi))))
    r = rounds_to(kdr, pk); ok &= r
    print("row [Dav1984], [Ree1991]: K_DR ~ %s: %s" % (pk, r))
    s2 = sqrt_fi(FI(2 * ONE))
    vals = [
        ("[G1953] upper, sinh(pi/2)", (exp_fi(PI / 2) - exp_fi(-(PI / 2))) / 2, ps),
        ("[G1953] lower, pi/2", PI / 2, "1.57080"),
        ("[K1979] pi/(2 log(1 + sqrt 2))", PI / (log_fi(s2 + FI(ONE)) * 2), "1.782214"),
        ("[SLXCKKM26] 6 pi/11", PI * 6 / 11, "1.71360"),
    ]
    for name, v, printed in vals:
        r = rounds_to(v, printed); ok &= r
        print("row %s = %s...: ~ %s: %s" % (name, show(v, 25), printed, r))
    print("VERDICT PASS" if ok else "VERDICT REFUSED: a printed value does not match")

if __name__ == "__main__":
    main()
