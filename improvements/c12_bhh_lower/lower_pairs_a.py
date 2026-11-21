#!/usr/bin/env python3
"""MST constant in the plane: certified lower bound with the isolated-pair term (program A).

Standard library only. Every real number is an interval [lo, hi] of integers at scale 2^-P, with every operation
rounded outward. Elementary functions are Taylor series with explicit remainder bounds.

Bound (see lower_pairs_README.md):
  beta >= H(a) + T1(a) + T2(a),
  H(a)  = int_0^a chi(r) dr = erf(a sqrt(pi)/2)/2 + (a/2) exp(-pi a^2/4),        chi = (1 - pi r^2/4) e^{-pi r^2/4}
  T1(a) = int_a^inf e^{-pi r^2} dr = erfc(a sqrt(pi))/2                           (isolated points)
  T2(a) = int_a^inf n2(r) dr = pi int_0^1 u G(B(u)) du                            (isolated pairs)
  B(u)  = pi + 2 asin(u/2) + (u/2) sqrt(4 - u^2)   (area of the union of two unit discs at distance u)
  G(B)  = a e^{-B a^2}/(2B) + sqrt(pi) erfc(a sqrt(B)) / (4 B^{3/2})  (= int_a^inf r^2 e^{-B r^2} dr)
B is increasing in u and G is decreasing in B, so on a cell [u0, u1]: int u G(B(u)) du >= G(B(u1)) (u1^2 - u0^2)/2.

Usage: python3 -I lower_pairs_a.py [--a 2679/2500] [--cells 20000] [--claim 0.63299] [--controls]
Exit status 0 iff the certified lower bound is >= the claim (and, with --controls, every red control is refused).
"""
import sys
from fractions import Fraction

P = 320
ONE = 1 << P


def fdiv(a, b):  # floor(a/b) for integers, b > 0
    return a // b


def cdiv(a, b):
    return -((-a) // b)


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi
        assert self.lo <= self.hi

    @staticmethod
    def frac(q):
        q = Fraction(q)
        return I(fdiv(q.numerator * ONE, q.denominator), cdiv(q.numerator * ONE, q.denominator))

    def __add__(s, o):
        o = o if isinstance(o, I) else I.frac(o)
        return I(s.lo + o.lo, s.hi + o.hi)

    __radd__ = __add__

    def __neg__(s):
        return I(-s.hi, -s.lo)

    def __sub__(s, o):
        o = o if isinstance(o, I) else I.frac(o)
        return I(s.lo - o.hi, s.hi - o.lo)

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        o = o if isinstance(o, I) else I.frac(o)
        p = [s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi]
        return I(fdiv(min(p), ONE), cdiv(max(p), ONE))

    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, I) else I.frac(o)
        assert o.lo > 0, "division by an interval that is not positive"
        p = [(s.lo, o.lo), (s.lo, o.hi), (s.hi, o.lo), (s.hi, o.hi)]
        return I(min(fdiv(x * ONE, y) for x, y in p), max(cdiv(x * ONE, y) for x, y in p))

    def __rtruediv__(s, o):
        return I.frac(o) / s

    def width(s):
        return Fraction(s.hi - s.lo, ONE)

    def lof(s):
        return Fraction(s.lo, ONE)

    def hif(s):
        return Fraction(s.hi, ONE)


def isqrt_floor(n):
    import math
    return math.isqrt(n)


def sqrt(x):
    assert x.lo >= 0
    lo = isqrt_floor(x.lo * ONE)
    hi = isqrt_floor(x.hi * ONE)
    if hi * hi < x.hi * ONE:
        hi += 1
    return I(lo, hi)


def exp_pos(x):
    """exp on an interval inside [0, 16]: Taylor series at each endpoint, remainder 2 x^N / N! once N > 2x."""
    assert x.lo >= 0 and x.hi <= 16 * ONE

    def at(v, up):
        # v: Fraction-free integer scaled value; sum_{n<N} v^n/n! with outward rounding, plus tail bound when up
        s = ONE
        term = ONE
        n = 0
        while True:
            n += 1
            term = (cdiv(term * v, ONE * n) if up else fdiv(term * v, ONE * n))
            s += term
            if n > 2 * (v // ONE + 1) and term < 2:  # term below 2^-P and series past its peak
                break
        if up:
            s += 2 * term + 2  # tail <= 2 * (last term) (ratio <= 1/2 from here) plus rounding slack
        return s

    return I(at(x.lo, False), at(x.hi, True))


def exp_neg(x):
    """exp(-x) for an interval x inside [0, 16]."""
    return I(ONE) / exp_pos(x)


def _alt_series(first_lo, first_hi, ratio_num, ratio_den, divisor, tiny):
    """Bracket of sum_n (-1)^n t_n / divisor(n) where t_0 = first, t_n = t_{n-1} * ratio_num(n) / ratio_den(n), all
    scaled integers. Terms are carried as lower and upper bounds; the limit lies between the last two partial sums
    once the terms decrease (checked by the caller's stopping rule)."""
    tl, tu = first_lo, first_hi
    slo = shi = 0
    hist = []
    n = 0
    prev = None
    while True:
        d = divisor(n)
        ql, qu = fdiv(tl, d), cdiv(tu, d)
        if n % 2 == 0:
            slo += ql; shi += qu
        else:
            slo -= qu; shi -= ql
        hist.append((slo, shi))
        if prev is not None and qu < prev and qu < tiny:
            break
        prev = qu
        n += 1
        num, den = ratio_num(n), ratio_den(n)
        tl, tu = fdiv(tl * num, den), cdiv(tu * num, den)
    (a_lo, a_hi), (b_lo, b_hi) = hist[-1], hist[-2]
    return min(a_lo, b_lo), max(a_hi, b_hi)


def atan_small(q):
    """atan(q), rational 0 < q <= 1/5, alternating series."""
    q = Fraction(q)
    tl, tu = fdiv(q.numerator * ONE, q.denominator), cdiv(q.numerator * ONE, q.denominator)
    q2n, q2d = q.numerator ** 2, q.denominator ** 2
    lo, hi = _alt_series(tl, tu, lambda n: q2n, lambda n: q2d, lambda n: 2 * n + 1, 4)
    return I(lo, hi)


def pi_interval():
    return 16 * atan_small(Fraction(1, 5)) - 4 * atan_small(Fraction(1, 239))


PI = pi_interval()


def _erf_series_point(V):
    """sum_n (-1)^n V^(2n+1)/(n!(2n+1)) for a scaled integer V >= 0 (V/ONE <= 4); returns a bracket."""
    V2 = V * V  # scale ONE^2
    return _alt_series(V, V, lambda n: V2, lambda n: n * ONE * ONE, lambda n: 2 * n + 1, 4)


def erf_pos(z):
    """erf on an interval z inside [0, 4], increasing: 2/sqrt(pi) times the alternating series at each end."""
    assert z.lo >= 0 and z.hi <= 4 * ONE
    c = I(2 * ONE) / sqrt(PI)
    l0, _ = _erf_series_point(z.lo)
    _, h1 = _erf_series_point(z.hi)
    return I((c * I(l0, l0)).lo, (c * I(h1, h1)).hi)


def asin_half(w):
    """asin on an interval w inside [0, 1/2]: sum C(2n,n)/(4^n (2n+1)) w^(2n+1), positive terms; the coefficient
    ratio is below 1, so the tail after the last term t is at most t w^2/(1 - w^2) <= t/3. asin is increasing."""
    assert w.lo >= 0 and 2 * w.hi <= ONE

    def series(V, up):
        p = V  # coef * w^(2n+1), scaled
        s = 0
        n = 0
        while True:
            term = cdiv(p, 2 * n + 1) if up else fdiv(p, 2 * n + 1)
            s += term
            if n > 2 and term < 4:
                break
            n += 1
            num = (2 * n - 1) * (2 * n) * V * V
            den = 4 * n * n * ONE * ONE
            p = cdiv(p * num, den) if up else fdiv(p * num, den)
        return s + (term + 4 if up else 0)

    return I(series(w.lo, False), series(w.hi, True))


def B_of(u):
    """Area of the union of two unit discs at distance u in [0, 1] (an interval)."""
    return PI + 2 * asin_half(u * Fraction(1, 2)) + (u * Fraction(1, 2)) * sqrt(I.frac(4) - u * u)


def G_of(a, Bv):
    """int_a^inf r^2 e^{-B r^2} dr for a point a (Fraction) and an interval B (positive)."""
    aI = I.frac(a)
    t1 = aI * exp_neg(Bv * (aI * aI)) / (2 * Bv)
    sB = sqrt(Bv)
    erfc = I(ONE) - erf_pos(aI * sB)
    t2 = sqrt(PI) * erfc / (4 * Bv * sB)
    return t1 + t2


def bound(a, cells):
    a = Fraction(a)
    aI = I.frac(a)
    sp = sqrt(PI)
    H = erf_pos(aI * sp * Fraction(1, 2)) * Fraction(1, 2) + aI * Fraction(1, 2) * exp_neg(PI * aI * aI * Fraction(1, 4))
    T1 = (I(ONE) - erf_pos(aI * sp)) * Fraction(1, 2)
    T2lo = Fraction(0)
    for i in range(cells):
        u0 = Fraction(i, cells)
        u1 = Fraction(i + 1, cells)
        Bv = B_of(I.frac(u1))
        Bup = I(Bv.hi)  # B(u) <= Bup on the cell, G decreasing in B
        g = G_of(a, Bup)
        w = (u1 * u1 - u0 * u0) / 2
        T2lo += Fraction((PI * g * w).lo, ONE)
    lo = H.lof() + T1.lof() + T2lo
    return lo, H, T1, T2lo


def main(argv):
    a = Fraction(2679, 2500)
    cells = 20000
    claim = Fraction("0.63299")
    controls = False
    i = 0
    while i < len(argv):
        if argv[i] == "--a":
            a = Fraction(argv[i + 1]); i += 2
        elif argv[i] == "--cells":
            cells = int(argv[i + 1]); i += 2
        elif argv[i] == "--claim":
            claim = Fraction(argv[i + 1]); i += 2
        elif argv[i] == "--controls":
            controls = True; i += 1
        else:
            raise SystemExit("unknown argument " + argv[i])
    # self-tests against known values
    def inside(x, lo, hi):
        return Fraction(lo) < x.lof() and x.hif() < Fraction(hi)
    assert PI.width() < Fraction(1, 10 ** 80) and inside(PI, "3.14159265358979323846264338327", "3.14159265358979323846264338328")
    assert inside(erf_pos(I(ONE)), "0.842700792949714869341220635082", "0.842700792949714869341220635083")
    s1 = asin_half(I.frac(Fraction(1, 2)))  # asin(1/2) = pi/6
    assert s1.lof() <= PI.hif() / 6 and PI.lof() / 6 <= s1.hif() and s1.width() < Fraction(1, 10 ** 60)
    assert inside(B_of(I.frac(1)), "5.0548156085708", "5.0548156085709")  # 4 pi/3 + sqrt(3)/2
    assert inside(exp_neg(I(ONE)), "0.367879441171442321595523770161", "0.367879441171442321595523770162")
    assert inside(erf_pos(I(3 * ONE)), "0.999977909503001414558627223870", "0.999977909503001414558627223871")
    print("self-tests: pi, erf(1), asin(1/2) = pi/6, B(1) = 4pi/3 + sqrt(3)/2, exp(-1): ok")
    lo, H, T1, T2 = bound(a, cells)
    print(f"a = {a}, cells = {cells}")
    print(f"H  = int_0^a chi        in [{float(H.lof()):.15f}, {float(H.hif()):.15f}]")
    print(f"T1 = isolated points    in [{float(T1.lof()):.15f}, {float(T1.hif()):.15f}]")
    print(f"T2 = isolated pairs     >= {float(T2):.15f}")
    print(f"certified: beta_MST >= {float(lo):.12f}  (exact lower end {lo.numerator.bit_length()}-bit rational)")
    ok = lo >= claim
    print(f"claim {float(claim)}: {'HOLDS' if ok else 'NOT IMPLIED'}")
    if controls:
        refused = 0
        # 1: claim one unit above in the 6th decimal of the certified value
        c1 = Fraction(int(lo * 10 ** 6) + 1, 10 ** 6)
        r1 = not (lo >= c1)
        print(f"red control 1 (claim {float(c1)}): {'refused' if r1 else 'ACCEPTED (bad)'}"); refused += r1
        # 2: a coarse grid (20 cells) must lose the claim
        lo2, *_ = bound(a, 20)
        r2 = not (lo2 >= claim)
        print(f"red control 2 (20 cells gives {float(lo2):.9f}): {'refused' if r2 else 'ACCEPTED (bad)'}"); refused += r2
        # 3: without the pair term the bound is the 14 Nov one (0.6315...), below the claim
        r3 = not (H.lof() + T1.lof() >= claim)
        print(f"red control 3 (no pair term, {float(H.lof() + T1.lof()):.9f}): {'refused' if r3 else 'ACCEPTED (bad)'}"); refused += r3
        # 4: cutoff far from the crossing (a = 1/2) loses the claim
        lo4, *_ = bound(Fraction(1, 2), 400)
        r4 = not (lo4 >= claim)
        print(f"red control 4 (a = 1/2 gives {float(lo4):.9f}): {'refused' if r4 else 'ACCEPTED (bad)'}"); refused += r4
        if refused != 4:
            print("RED CONTROLS FAILED"); return 2
        print("red controls: 4 of 4 refused")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
