#!/usr/bin/env python3
"""Program B for 48a: C_48 = c_0^2 from the system of Davis and Power (arXiv 2604.03170v1, (3)-(6)), by a second
route that shares no code with c48a_a.py: interval arithmetic on Python's decimal module at 80 digits, every operation
rounded outward (exp, ln and sqrt are correctly rounded in decimal, so one unit in the last place encloses them).

    t_0 = sqrt(2 ln 2);  A = t_0 + 2 Q(t_0);  B = A/2;  H(a) = 2a e^(-a^2/2) + 2 Q(a) = B (a > t_0);
    p_0 = 2 e^(-a^2/2);  Q(z)/sqrt(2 pi) = p_0;  c_0 = B sqrt(2 pi) e^(z^2/2);  C_48 = c_0^2,
with Q(x) = int_x^inf e^(-t^2/2) dt = sqrt(pi/2) - e^(-x^2/2) sum_{n >= 0} x^(2n+1) / (1 * 3 * ... * (2n+1)) (positive
terms; once x^2/(2n+3) <= 1/2 the rest is at most the last term), and pi by the Bailey-Borwein-Plouffe series (positive
terms after grouping; the tail after k terms is below 16^-k). a and z by interval bisection.

  c48a_b.py [--claim-round DEC] [--variant t0=sqrt(ln2) | p0-without-2 | H-with-Q-once]
Standard library only; run with python3 -I.
"""
import sys
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING

P = 80
DN = Context(prec=P, rounding=ROUND_FLOOR)
UP = Context(prec=P, rounding=ROUND_CEILING)
NEAR = Context(prec=P)


class J:
    """A closed interval of decimals."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = D(lo)
        self.hi = D(lo) if hi is None else D(hi)
        assert self.lo <= self.hi

    def __add__(s, o):
        o = o if isinstance(o, J) else J(o)
        return J(DN.add(s.lo, o.lo), UP.add(s.hi, o.hi))

    def __sub__(s, o):
        o = o if isinstance(o, J) else J(o)
        return J(DN.subtract(s.lo, o.hi), UP.subtract(s.hi, o.lo))

    def __mul__(s, o):
        o = o if isinstance(o, J) else J(o)
        c = [(x, y) for x in (s.lo, s.hi) for y in (o.lo, o.hi)]
        return J(min(DN.multiply(x, y) for x, y in c), max(UP.multiply(x, y) for x, y in c))

    def __truediv__(s, o):
        o = o if isinstance(o, J) else J(o)
        assert o.lo > 0 or o.hi < 0
        c = [(x, y) for x in (s.lo, s.hi) for y in (o.lo, o.hi)]
        return J(min(DN.divide(x, y) for x, y in c), max(UP.divide(x, y) for x, y in c))


def ulp(x):
    return D(10) ** (x.adjusted() - P + 1)


def j_exp(x):
    """e^x, increasing: correctly rounded values at both ends, widened by one ulp."""
    lo, hi = NEAR.exp(x.lo), NEAR.exp(x.hi)
    return J(DN.subtract(lo, ulp(lo)), UP.add(hi, ulp(hi)))


def j_sqrt(x):
    return J(DN.sqrt(x.lo), UP.sqrt(x.hi))


def j_sq(x):
    assert x.lo >= 0
    return J(DN.multiply(x.lo, x.lo), UP.multiply(x.hi, x.hi))


def pi_bbp():
    s = J(0)
    k = 0
    while True:
        term = (J(4) / J(8 * k + 1) - J(2) / J(8 * k + 4) - J(1) / J(8 * k + 5) - J(1) / J(8 * k + 6)) / J(16 ** k)          # 16^k as an exact integer
        s = s + term
        k += 1
        if 16 ** k > 10 ** (P - 2):
            # tail after k terms: each later term is positive and below 16^-j, so the tail is below 16^-k * 16/15
            return s + J(0, UP.divide(16, UP.multiply(15, 16 ** k)))


PI = pi_bbp()
SQ_PI_2 = j_sqrt(PI / J(2))
SQ_2PI = j_sqrt(PI * J(2))


def Q(x):
    """int_x^inf e^(-t^2/2) dt for an interval x >= 0 (decreasing in x)."""
    def at(v):
        vv = J(v)
        x2 = j_sq(vv)
        term = vv
        s = J(0)
        n = 0
        while True:
            s = s + term
            n += 1
            ratio_hi = UP.divide(x2.hi, 2 * n + 1)
            term = term * x2 / J(2 * n + 1)
            if ratio_hi <= D("0.5") and term.hi < D(10) ** (-(P - 5)):
                s = s + J(0, UP.multiply(term.hi, 2))      # the rest is at most twice the current term
                break
        return SQ_PI_2 - j_exp(J(0) - x2 * J("0.5")) * s
    a, b = at(x.hi), at(x.lo)
    return J(a.lo, b.hi)


def H(x):
    return x * J(2) * j_exp(J(0) - j_sq(x) * J("0.5")) + Q(x) * J(2)


def bisect(f, lo, hi, target):
    """f decreasing; returns [lo, hi] with f(lo) > target > f(hi), certified."""
    lo, hi = D(lo), D(hi)
    while UP.subtract(hi, lo) > D("1e-70"):
        mid = DN.divide(DN.add(lo, hi), 2)
        v = f(J(mid))
        if v.lo > target.hi:
            lo = mid
        elif v.hi < target.lo:
            hi = mid
        else:
            break
    assert f(J(lo)).lo > target.hi and f(J(hi)).hi < target.lo, "bracket not certified"
    return J(lo, hi)


def show(x, d=30):
    q = D(10) ** -d
    return "[%s, %s]" % (Context(prec=200, rounding=ROUND_FLOOR).quantize(x.lo, q),
                         Context(prec=200, rounding=ROUND_CEILING).quantize(x.hi, q))


def main():
    claim = variant = None
    a_ = sys.argv[1:]
    while a_:
        if a_[0] == "--claim-round":
            claim = a_[1]; a_ = a_[2:]
        elif a_[0] == "--variant":
            variant = a_[1]; a_ = a_[2:]; print("CONTROL variant " + variant)
        else:
            raise SystemExit("unknown option " + a_[0])
    l2 = NEAR.ln(D(2))
    LN2 = J(DN.subtract(l2, ulp(l2)), UP.add(l2, ulp(l2)))
    t0 = j_sqrt(LN2 * J(2)) if variant != "t0=sqrt(ln2)" else j_sqrt(LN2)
    A = t0 + Q(t0) * J(2)
    B = A * J("0.5")
    print("pi in %s" % show(PI))
    print("t_0 in %s" % show(t0))
    print("A in %s" % show(A))
    Hf = (lambda v: H(v) - Q(v)) if variant == "H-with-Q-once" else H
    a = bisect(Hf, t0.hi, D(10), B)
    print("a in %s" % show(a))
    p0 = j_exp(J(0) - j_sq(a) * J("0.5")) * J(1 if variant == "p0-without-2" else 2)
    print("p_0 in %s" % show(p0))
    z = bisect(lambda v: Q(v) / SQ_2PI, D(0), D(3), p0)
    print("z in %s" % show(z))
    c0 = B * SQ_2PI * j_exp(j_sq(z) * J("0.5"))
    C48 = j_sq(c0)
    print("c_0 in %s" % show(c0))
    print("C_48 = c_0^2 in %s" % show(C48))
    ok = True
    if claim is not None:
        d = len(claim.split(".")[1])
        half = D(5) * D(10) ** (-(d + 1))
        cl = D(claim)
        holds = DN.subtract(cl, half) <= C48.lo and C48.hi < UP.add(cl, half)
        print("c_0^2 rounds to %s at %d decimals: %s" % (claim, d, "yes" if holds else "NO"))
        ok = ok and holds
    for name, val, printed in (("a", a, "1.80334"), ("p_0", p0, "0.39342"), ("z", z, "0.27041"), ("c_0", c0, "2.30952")):
        d = len(printed.split(".")[1]); half = D(5) * D(10) ** (-(d + 1)); cl = D(printed)
        r = DN.subtract(cl, half) <= val.lo and val.hi < UP.add(cl, half)
        print("information: %s rounds to the paper's %s: %s" % (name, printed, "yes" if r else "no"))
    print("VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
