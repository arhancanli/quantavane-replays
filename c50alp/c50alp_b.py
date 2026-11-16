#!/usr/bin/env python3
"""50a [LP24] 0.595, program B (copied from c50apt_b.py, itself from c50a/c50a_b.py: closed forms of F_3 through arcsin
in decimal intervals), with the ratio replaced by Lee and Parekh's
    rho_p(x) = [p (1 - F_3(x)) + (1 - p)(1 + (6/5) max(0, -(1 + 3x)))] / (1 - 3x)
(arXiv 2401.03616, Theorem 2) and p = 0.674 (--p P). Cells as in program A. Original header follows.

50a, program B (written separately from program A): [GP19] 0.498 for Quantum Max Cut product states
(arXiv 1909.08846v1, Theorem 4), recomputed from its definition (a value, not a certificate).

The ratio is min over t in [-1, 1/r) of (1 - F_r(t))/(1 - r t), F_r(t) = (2/r)(Gamma((r+1)/2)/Gamma(r/2))^2 t
2F1(1/2, 1/2; r/2 + 1; t^2) (Briet, de Oliveira Filho and Vallentin). Route B uses closed forms instead of the series:
Gauss's contiguous relation c(c-1)(z-1) F(c-1) + c(c-1-(2c-a-b-1)z) F(c) + (c-a)(c-b) z F(c+1) = 0 with
a = b = 1/2, c = 3/2, 2F1(1/2,1/2;1/2;z) = (1-z)^(-1/2) and 2F1(1/2,1/2;3/2;t^2) = arcsin(t)/t gives
    F_3(t) = (2 / (pi t)) [ sqrt(1 - t^2) - (1 - 2 t^2) arcsin(t) / t ],   F_3(-1) = -1,   F_3(0) = 0,
and F_1(t) = (2/pi) arcsin t (the Goemans-Williamson anchor). Arithmetic: decimal intervals at 50 digits with
directed rounding; sqrt of the decimal module widened by one unit in the last place; arctan by halving the argument
(atan x = 2 atan(x / (1 + sqrt(1 + x^2)))) and its alternating series; pi = 4 atan 1. A uniform grid of cells on
[-1, 0]; F_r increases, so on [a, b] the ratio is at least (1 - F_r(b)) / (1 - r a); on [0, 1/r) it is at least
1 - F_r(1/r).
Usage: c50a_b.py [--r 1|3] [--target T] [--cells M] [--g-unsquared]
"""
import sys
from decimal import Decimal as D, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

getcontext().prec = 50

def lo(f):
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        return f()
def hi(f):
    with localcontext() as c:
        c.rounding = ROUND_CEILING
        return f()

class J:
    """interval [a, b] of decimals"""
    def __init__(self, a, b=None):
        self.a, self.b = a, (a if b is None else b)
        assert self.a <= self.b
    def __add__(s, o):
        return J(lo(lambda: s.a + o.a), hi(lambda: s.b + o.b))
    def __sub__(s, o):
        return J(lo(lambda: s.a - o.b), hi(lambda: s.b - o.a))
    def __mul__(s, o):
        ps = [(s.a, o.a), (s.a, o.b), (s.b, o.a), (s.b, o.b)]
        return J(min(lo(lambda: x * y) for x, y in ps), max(hi(lambda: x * y) for x, y in ps))
    def __truediv__(s, o):
        assert o.a > 0 or o.b < 0
        qs = [(s.a, o.a), (s.a, o.b), (s.b, o.a), (s.b, o.b)]
        return J(min(lo(lambda: x / y) for x, y in qs), max(hi(lambda: x / y) for x, y in qs))
    def neg(s):
        return J(-s.b, -s.a)

def jsqrt(x):
    assert x.a >= 0
    with localcontext() as c:
        c.rounding = ROUND_HALF_EVEN
        ra, rb = x.a.sqrt(), x.b.sqrt()
    return J(max(D(0), ra.next_minus()), rb.next_plus())

ONE = J(D(1))

def jatan_small(x):
    """atan for |x| <= 0.1: alternating series sum (-1)^k x^(2k+1)/(2k+1) with the first omitted term as error"""
    s = J(D(0))
    p = x
    x2 = x * x
    k = 0
    while True:
        term = p / J(D(2 * k + 1))
        s = s + term if k % 2 == 0 else s - term
        k += 1
        p = p * x2
        nxt = p / J(D(2 * k + 1))
        bound = max(abs(nxt.a), abs(nxt.b))
        if bound < D("1e-48"):
            return J(lo(lambda: s.a - bound), hi(lambda: s.b + bound))

def jatan(x):
    """atan of an interval x, by halving until |x| <= 0.1 (atan is increasing, so the ends are treated apart)"""
    def at_point(v):
        y = J(v)
        doublings = 0
        while max(abs(y.a), abs(y.b)) > D("0.1"):
            y = y / (ONE + jsqrt(ONE + y * y))
            doublings += 1
        r = jatan_small(y)
        for _ in range(doublings):
            r = r + r
        return r
    ra, rb = at_point(x.a), at_point(x.b)
    return J(ra.a, rb.b)

PI = J(D(4)) * jatan(ONE)

def jarcsin(t):
    """arcsin for a point t in (-1, 1): atan(t / sqrt(1 - t^2))"""
    T = J(t)
    return jatan(T / jsqrt(ONE - T * T))

args = sys.argv[1:]
r, target, M, gunsq, LOW, PM = 3, D("0.595"), 4000, False, D(-1), D("0.674")
while args:
    k = args.pop(0)
    if k == "--r":
        r = int(args.pop(0))
    elif k == "--target":
        target = D(args.pop(0))
    elif k == "--cells":
        M = int(args.pop(0))
    elif k == "--g-unsquared":
        gunsq = True
    elif k == "--lo":
        LOW = D(args.pop(0))
    elif k == "--p":
        PM = D(args.pop(0))
    else:
        raise SystemExit("bad argument")

def F(t):
    """an interval containing F_r(t) for a decimal point t in [-1, 1]"""
    if r == 1:
        if t == -1:
            return J(D(-1))
        return J(D(2)) * jarcsin(t) / PI
    if r == 3:
        if t == 0:
            return J(D(0))
        if t == -1:
            val = J(D(-1))
        else:
            T = J(t)
            bracket = jsqrt(ONE - T * T) - (ONE - J(D(2)) * T * T) * jarcsin(t) / T
            val = J(D(2)) * bracket / (PI * T)
        if gunsq:   # the factor (Gamma(2)/Gamma(3/2))^2 = 4/pi replaced by 2/sqrt(pi): multiply by sqrt(pi)/2
            val = val * jsqrt(PI) / J(D(2))
        return val
    raise SystemExit("r = 1 or 3")

print("program B (c50a_b.py): closed forms through arcsin, decimal intervals")
if gunsq:
    print("CONTROL variant: (Gamma(2)/Gamma(3/2)) not squared")
print("pi in [%s, %s]" % (format(PI.a, ".30f"), format(PI.b, ".30f")))
if r == 3:
    print("F_3(-1/2) in [%s, %s]" % (format(F(D("-0.5")).a, ".20f"), format(F(D("-0.5")).b, ".20f")))
if r == 1:
    # arcsin is convex on [0, 1] with arcsin 0 = 0, arcsin 1 = pi/2, so F_1(t) <= t and the ratio is >= 1 there
    pos_low = D(1)
else:
    F1r = F(D(1) / D(r))
    pos_low = lo(lambda: PM * (1 - F1r.b) + (1 - PM))
print("p = %s; t in [0, 1/%d): ratio >= p (1 - F(1/%d)) + (1 - p) >= %s" % (PM, r, r, format(pos_low, ".6f")))
ok = pos_low >= target
least, where = None, None
step = (0 - LOW) / D(M)
queue = [(LOW + step * i, LOW + step * (i + 1)) for i in range(M)]
FC = {}
def Fcached(t):
    if t not in FC:
        FC[t] = F(t)
    return FC[t]
nc = 0; failed = None
while queue:
    a_, b_ = queue.pop()
    Fb = Fcached(b_)
    den = hi(lambda: 1 - r * a_)
    num = lo(lambda: PM * (1 - Fb.b) + (1 - PM) * (1 + D(6) / D(5) * max(D(0), -(1 + 3 * b_))))
    lb = lo(lambda: num / den)
    if lb >= target:
        nc += 1
        if least is None or lb < least:
            least, where = lb, (a_, b_)
        continue
    if b_ - a_ < D(2) ** -30:
        failed = (a_, b_, lb); break
    mid = (a_ + b_) / 2
    queue += [(a_, mid), (mid, b_)]
if failed:
    ok = False
    print("t in [%s, 0]: a cell of width < 2^-30 at %s does not clear the target (lower bound %s)" % (LOW, format(failed[0], ".9f"), format(failed[2], ".9f")))
else:
    ok = ok and least >= target
    print("t in [%s, 0], %d equal cells refined adaptively to %d: least cell lower bound %s on [%s, %s]"
          % (LOW, M, nc, format(least, ".9f"), format(where[0], ".9f"), format(where[1], ".9f")))
print("RESULT p = %s: min over [%s, 1/3) of rho_p >= %s: %s" % (PM, LOW, target, ok))
print("VERDICT", "PASS" if ok else "REFUSED")
