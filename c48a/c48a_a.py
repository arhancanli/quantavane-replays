#!/usr/bin/env python3
"""Program A for 48a (the one-dimensional convex sub-Gaussian comparison constant): C_48 = c_0^2 from the system of
equations of Davis and Power (arXiv 2604.03170v1, (3)-(6) and Theorem 1), in ball arithmetic on exact rationals
(centre and radius, rounded to 2^-PREC with the rounding added to the radius).

The system (with Q(x) = int_x^inf e^(-t^2/2) dt, phi the standard normal density, Phibar(z) = Q(z)/sqrt(2 pi)):
    t_0 = sqrt(2 ln 2);  A = t_0 + 2 Q(t_0);  B = A/2;
    H(x) = 2x e^(-x^2/2) + 2 Q(x) for x >= t_0, decreasing; a > t_0 the root of H(a) = B;
    p_0 = 2 e^(-a^2/2);  z the root of Phibar(z) = p_0;  c_0 = B / phi(z);  C_48 = c_0^2.
Route (ours): Q(x) = sqrt(pi/2) - sum_n (-1)^n x^(2n+1) / (2^n n! (2n+1)) (alternating; error below the first omitted
term once the terms decrease); pi by Machin's formula; ln 2 = 2 atanh(1/3); exp by its Taylor series with a geometric
tail; square roots by integer square roots. a and z by bisection on rationals, keeping a bracket whose two ends are
certified on opposite sides.

  c48a_a.py [--claim-round DEC] [--variant t0=sqrt(ln2) | p0-without-2 | H-with-Q-once]
--claim-round D checks that c_0^2 rounds to D at D's number of decimals (the page prints "c_0^2 ~ 5.33386");
--variant changes the system (red controls). Standard library only; run with python3 -I.
"""
import sys
from fractions import Fraction as Fr
from math import isqrt

PREC = 220
ONE = 1 << PREC
EPS = Fr(1, 1 << (PREC - 8))


class Ball:
    __slots__ = ("c", "r")

    def __init__(self, c, r=0):
        c, r = Fr(c), Fr(r)
        if c.denominator > ONE:
            cc = Fr((c.numerator * ONE) // c.denominator, ONE)
            r += c - cc
            c = cc
        if r.denominator > ONE:
            r = Fr(-((-r.numerator * ONE) // r.denominator), ONE)
        self.c, self.r = c, r

    def lo(self):
        return self.c - self.r

    def hi(self):
        return self.c + self.r

    def __add__(s, o):
        o = o if isinstance(o, Ball) else Ball(o)
        return Ball(s.c + o.c, s.r + o.r)
    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, Ball) else Ball(o)
        return Ball(s.c - o.c, s.r + o.r)

    def __rsub__(s, o):
        return Ball(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, Ball) else Ball(o)
        return Ball(s.c * o.c, abs(s.c) * o.r + abs(o.c) * s.r + s.r * o.r)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, Ball) else Ball(o)
        m = abs(o.c) - o.r
        assert m > 0
        inv = Ball(1 / o.c, o.r / (abs(o.c) * m))
        return s * inv


def sqrt_ball(x):
    lo, hi = x.lo(), x.hi()
    assert lo >= 0
    a = Fr(isqrt((lo.numerator * ONE * ONE) // lo.denominator), ONE)
    b = Fr(isqrt(-((-hi.numerator * ONE * ONE) // hi.denominator)) + 1, ONE)
    return Ball((a + b) / 2, (b - a) / 2)


def atan_inv(n):
    s, k = Fr(0), 0
    while True:
        t = Fr(1, (2 * k + 1) * n ** (2 * k + 1))
        if t < EPS:
            return Ball(s, t)          # alternating, decreasing: the error is below the first omitted term
        s += t if k % 2 == 0 else -t
        k += 1


PI = atan_inv(5) * 16 - atan_inv(239) * 4


def ln2():
    s, k = Fr(0), 0
    while True:
        t = Fr(1, (2 * k + 1) * 3 ** (2 * k + 1))
        tail = t * Fr(9, 8)
        if tail < EPS:
            return Ball(s, tail) * 2
        s += t
        k += 1


def exp_point(x):
    """e^x for a rational x with |x| < 20."""
    s, t, k = Fr(0), Fr(1), 0
    while not (k > 2 * abs(x) + 1 and abs(t) < EPS):
        s += t
        k += 1
        t = t * x / k
    return Ball(s, 2 * abs(t))


def exp_ball(x):
    v = exp_point(x.c)
    return v + Ball(0, x.r * 3 ** (int(abs(x.c) + x.r) + 1))


def P0(x):
    """int_0^x e^(-t^2/2) dt for a rational x >= 0 (alternating series)."""
    s, n, p, fact2 = Fr(0), 0, x, Fr(1)        # p = x^(2n+1), fact2 = 2^n n!
    while True:
        t = p / (fact2 * (2 * n + 1))
        if n > x * x and t < EPS:
            return Ball(s, t)
        s += t if n % 2 == 0 else -t
        n += 1
        p *= x * x
        fact2 *= 2 * n


SQ_PI_2 = sqrt_ball(PI / 2)


def Q_ball(x):
    """Q over a ball x >= 0: Q(centre) plus radius * max|Q'| (Q' = -e^(-x^2/2), |Q'| <= 1)."""
    return SQ_PI_2 - P0(x.c) + Ball(0, x.r)


def H_ball(x):
    e = exp_ball(Ball(-(x.c * x.c), 2 * abs(x.c) * x.r + x.r * x.r) * Fr(1, 2))
    return x * e * 2 + Q_ball(x) * 2


def bisect(f_ball, lo, hi, target, decreasing=True):
    """Root of f = target (f monotone) bracketed on rationals; target is a Ball."""
    for _ in range(PREC - 30):
        mid = Fr(((lo + hi) / 2 * ONE).__floor__(), ONE)
        v = f_ball(Ball(mid))
        above = v.lo() > target.hi()
        below = v.hi() < target.lo()
        if not (above or below):
            break
        if above == decreasing:
            lo = mid
        else:
            hi = mid
    va, vb = f_ball(Ball(lo)), f_ball(Ball(hi))
    if decreasing:
        assert va.lo() > target.hi() and vb.hi() < target.lo(), "bracket not certified"
    return Ball((lo + hi) / 2, (hi - lo) / 2)


def show(x, d=30):
    lo, hi = x.lo(), x.hi()
    f = lambda v: (v.numerator * 10 ** d) // v.denominator
    g = lambda v: -((-v.numerator * 10 ** d) // v.denominator)
    def s(q):
        t = str(abs(q)).rjust(d + 1, "0")
        return ("-" if q < 0 else "") + t[:-d] + "." + t[-d:]
    return "[%s, %s]" % (s(f(lo)), s(g(hi)))


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
    L2 = ln2()
    t0 = sqrt_ball(L2 * 2) if variant != "t0=sqrt(ln2)" else sqrt_ball(L2)
    A = t0 + Q_ball(t0) * 2
    B = A * Fr(1, 2)
    print("pi in %s" % show(PI))
    print("t_0 in %s" % show(t0))
    print("A in %s" % show(A))
    if variant == "H-with-Q-once":
        Hb = lambda x: H_ball(x) - Q_ball(x)
    else:
        Hb = lambda x: H_ball(x)
    a = bisect(Hb, t0.hi(), Fr(10), B, decreasing=True)       # H(t_0) = A > B and H(10) < B
    print("a in %s" % show(a))
    p0 = exp_ball(Ball(-(a.c * a.c), 2 * abs(a.c) * a.r + a.r * a.r) * Fr(1, 2)) * (1 if variant == "p0-without-2" else 2)
    print("p_0 in %s" % show(p0))
    sq2pi = sqrt_ball(PI * 2)
    Phibar = lambda z: Q_ball(z) / sq2pi
    z = bisect(Phibar, Fr(0), Fr(3), p0, decreasing=True)
    print("z in %s" % show(z))
    # c_0 = B / phi(z) = B sqrt(2 pi) e^(z^2/2), increasing in z > 0: evaluate at both ends of z's bracket
    def c0_at(zz):
        return B * sq2pi * exp_ball(Ball(zz * zz / 2))
    lo_v, hi_v = c0_at(z.lo()), c0_at(z.hi())
    c0 = Ball((lo_v.lo() + hi_v.hi()) / 2, (hi_v.hi() - lo_v.lo()) / 2)
    C48 = Ball((c0.lo() ** 2 + c0.hi() ** 2) / 2, (c0.hi() ** 2 - c0.lo() ** 2) / 2)
    print("c_0 in %s" % show(c0))
    print("C_48 = c_0^2 in %s" % show(C48))
    ok = True
    if claim is not None:
        d = len(claim.split(".")[1])
        half = Fr(1, 2 * 10 ** d)
        cl = Fr(claim)
        holds = cl - half <= C48.lo() and C48.hi() < cl + half
        print("c_0^2 rounds to %s at %d decimals: %s" % (claim, d, "yes" if holds else "NO"))
        ok &= holds
    for name, val, printed in (("a", a, "1.80334"), ("p_0", p0, "0.39342"), ("z", z, "0.27041"), ("c_0", c0, "2.30952")):
        d = len(printed.split(".")[1]); half = Fr(1, 2 * 10 ** d); cl = Fr(printed)
        r = cl - half <= val.lo() and val.hi() < cl + half
        print("information: %s rounds to the paper's %s: %s" % (name, printed, "yes" if r else "no"))
    print("VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
