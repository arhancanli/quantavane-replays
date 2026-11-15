#!/usr/bin/env python3
"""10a, program B (written separately from program A): the same five rows (see program A for the statements).

Route B: decimal intervals at 60 digits (exp, ln, sqrt of the decimal module, correctly rounded, widened by one unit
in the last place; Phi by its everywhere-convergent series at 110 digits; pi by Machin in rationals). C* by bisection
on H(C) = 4 phi(C)^2 - 4 Phi(-C) + 1 in decimal points, and K_DR = 1/R(C) evaluated over the final C-interval with
the FULL formula R(C) = (4 phi^2 - 2 C phi (4 Phi(-C) - 1))/(1 - 2 C phi) (no use of the identity at the root).
sinh through exp, log(1 + sqrt 2) through ln. Standard library only.
Usage: c10a_b.py [--print-kdr D] [--print-sinh D] [--variant phi-plus]
"""
import sys
from fractions import Fraction as Q
from decimal import Decimal as D, localcontext, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

DPREC = 60
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
_pi_hp = machin_pi_exact(130)
INV_SQRT_2PI = sqrt_di(PI * 2).recip()

def Phi_point(x):
    """Phi(x) for a decimal x: 1/2 + phi(x) sum_n x^(2n+1)/(2n+1)!!, at 110 digits; tail <= 2 |last term|"""
    P = 110
    with _ctx(ROUND_HALF_EVEN, P) as c:
        x2 = x * x
        term = x; s = x; n = 0
        while True:
            n += 1
            term = term * x2 / (2 * n + 1)
            s += term
            if 2 * n + 3 > 2 * x2 and abs(term) < D(10) ** -(P - 5):
                break
        pi = D(_pi_hp.numerator) / D(_pi_hp.denominator)       # within 10^-125 of pi
        ph = (-x2 / 2).exp() / (2 * pi).sqrt()
        v = D(1) / 2 + ph * s
        err = ph * abs(term) * 2 + abs(v) * D(10) ** -(P - 10) + D(10) ** -(P - 10)
    with _ctx(ROUND_FLOOR): lo = v - err
    with _ctx(ROUND_CEILING): hi = v + err
    return DI(lo, hi)


def phi_at(x):
    """phi at a decimal interval x with 0 <= x"""
    x2 = x * x
    return exp_di(DI(x2.hi.copy_negate(), x2.lo.copy_negate()) / 2) * INV_SQRT_2PI

def Phi_neg(x):
    """Phi(-x) for a decimal interval x (decreasing in x)"""
    return DI(Phi_point(x.hi.copy_negate()).lo, Phi_point(x.lo.copy_negate()).hi)

def H(x, variant):
    ph = phi_at(x)
    P = Phi_neg(x) if variant != "phi-plus" else DI.q(1) - Phi_neg(x)
    return ph * ph * 4 - P * 4 + 1

def R(x):
    ph = phi_at(x)
    lam = x * ph * 2
    num = ph * ph * 4 - lam * (Phi_neg(x) * 4 - 1)
    return num / (DI.q(1) - lam)

def rounds_to(v, printed):
    d = len(printed.split(".")[1])
    p = D(printed); h = D(5) / D(10) ** (d + 1)
    with localcontext(Context(prec=80)):
        return p - h <= v.lo and v.hi <= p + h

def main():
    a = sys.argv[1:]
    opt = dict(zip(a[0::2], a[1::2]))
    variant = opt.get("--variant")
    pk = opt.get("--print-kdr", "1.67696"); ps = opt.get("--print-sinh", "2.30130")
    if variant:
        print("CONTROL: H with Phi(C) in place of Phi(-C)")
    print("10a program B")
    with _ctx(ROUND_HALF_EVEN, 80):
        p1 = Phi_point(D(1)); t1 = abs(p1.lo - D("0.841344746068542948585232545632")) < D(10) ** -29
        t2 = abs(PI.lo - D("3.14159265358979323846264338327950288")) < D(10) ** -34
    print("self-test Phi(1), pi: %s" % (t1 and t2))
    if not (t1 and t2):
        print("VERDICT FAILED (self-test)"); return
    lo, hi = D(0), D(1)
    hl, hh = H(DI(lo), variant), H(DI(hi), variant)
    print("H(0) in [%.12f, %.12f]; H(1) in [%.12f, %.12f]" % (hl.lo, hl.hi, hh.lo, hh.hi))
    if not (hl.hi < 0 < hh.lo):
        print("no sign change of H on [0, 1]")
        print("VERDICT REFUSED: C* not found"); return
    for _ in range(170):
        with _ctx(ROUND_HALF_EVEN, 80):
            mid = (lo + hi) / 2
        hm = H(DI(mid), variant)
        if hm.hi < 0: lo = mid
        elif hm.lo > 0: hi = mid
        else: break
    C = DI(lo, hi)
    print("C* in [%s, %s]" % (str(lo)[:34], str(hi)[:34]))
    kdr = R(C).recip()
    print("K_DR = 1/R(C*) in [%s, %s]" % (str(kdr.lo)[:34], str(kdr.hi)[:34]))
    ok = rounds_to(kdr, pk)
    print("row [Dav1984], [Ree1991]: K_DR ~ %s: %s" % (pk, ok))
    e = exp_di(PI / 2); s2 = sqrt_di(DI.q(2))
    vals = [
        ("[G1953] upper, sinh(pi/2)", (e - e.recip()) / 2, ps),
        ("[G1953] lower, pi/2", PI / 2, "1.57080"),
        ("[K1979] pi/(2 log(1 + sqrt 2))", PI / (ln_di(s2 + 1) * 2), "1.782214"),
        ("[SLXCKKM26] 6 pi/11", PI * 6 / 11, "1.71360"),
    ]
    for name, v, printed in vals:
        r = rounds_to(v, printed); ok &= r
        print("row %s = %s...: ~ %s: %s" % (name, str(v.lo)[:27], printed, r))
    print("VERDICT PASS" if ok else "VERDICT REFUSED: a printed value does not match")

if __name__ == "__main__":
    main()
