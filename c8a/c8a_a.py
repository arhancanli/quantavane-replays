#!/usr/bin/env python3
"""8a: the lower bound C_8 >= 2/log(gamma_1) ~ 0.755106, program A.

The page: a zero of zeta at 1/2 + i gamma with gamma >= 2 forces 1/2 <= 1 - 1/(R log gamma), i.e. R >= 2/log gamma; the
row uses gamma_1 = 14.1347... . This program proves that zeta has a zero 1/2 + i gamma with gamma in a short interval
(t_lo, t_hi), by a sign change of Hardy's real function Z(t) = e^{i theta(t)} zeta(1/2 + it), and then encloses
2/log gamma. (A zero at any ordinate gamma >= 2 gives the bound R >= 2/log gamma; that it is the lowest zero is not
needed for the bound and is not checked here.)

- zeta(s) by Euler-Maclaurin with N = 40 and 30 Bernoulli terms; remainder bounded by
  |B_2p|/(2p)! |s(s+1)...(s+2p-1)| N^(1-sigma-2p)/(sigma+2p-1) (the periodic Bernoulli function is at most |B_2p|).
- e^{i theta(t)} = pi^(-it/2) Gamma(1/4 + it/2)/|Gamma(1/4 + it/2)|: only Im log Gamma is needed, by Stirling's series at
  w = z + 20 (20 terms; remainder |B_42|/(42*41 |w|^41) sec^42(arg(w)/2)) and the shift log Gamma(z) = log Gamma(z+20)
  - sum log(z+j). The branch does not matter (only e^{i theta} enters).
- Arithmetic: closed intervals of integers at scale 2^-320 with outward rounding; our own sqrt (isqrt), ln (atanh
  series), atan (Taylor, |x| < 1/2), cos and sin (reduction mod 2 pi, Taylor with remainder), pi by Machin.
Self-test: Im(e^{i theta} zeta(1/2 + it)) must contain 0 at every point (Z is real). Standard library only.
"""
import sys
from fractions import Fraction as Fr
from math import isqrt

P = 320
ONE = 1 << P
EPS_STOP = Fr(1, 1 << (P - 8))     # series stop once a term is below 2^-312 (outward rounding keeps 1 unit)


class R:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo; self.hi = lo if hi is None else hi
        assert self.lo <= self.hi

    @staticmethod
    def q(x):
        x = Fr(x); n = x.numerator * ONE
        return R(n // x.denominator, -((-n) // x.denominator))

    def __add__(a, b): b = b if isinstance(b, R) else R.q(b); return R(a.lo + b.lo, a.hi + b.hi)
    def __radd__(a, b): return a + b
    def __sub__(a, b): b = b if isinstance(b, R) else R.q(b); return R(a.lo - b.hi, a.hi - b.lo)
    def __neg__(a): return R(-a.hi, -a.lo)

    def __mul__(a, b):
        if not isinstance(b, R):
            b = Fr(b)
            if b.denominator == 1:
                k = b.numerator
                return R(min(a.lo * k, a.hi * k), max(a.lo * k, a.hi * k))
            b = R.q(b)
        c = (a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi)
        return R(min(c) >> P, -((-max(c)) >> P))
    __rmul__ = __mul__

    def inv(a):
        assert a.lo > 0 or a.hi < 0, "division by an interval containing 0"
        lo = ONE * ONE // a.hi if True else 0
        # 1/[lo, hi] = [1/hi, 1/lo]
        return R((ONE * ONE) // a.hi, -((-ONE * ONE) // a.lo))

    def __truediv__(a, b):
        if not isinstance(b, R):
            b = Fr(b)
            if b.denominator == 1 and b.numerator > 0:
                k = b.numerator
                return R(a.lo // k, -((-a.hi) // k))
            b = R.q(b)
        return a * b.inv()

    def width(a): return Fr(a.hi - a.lo, ONE)
    def mid_fr(a): return Fr(a.lo + a.hi, 2 * ONE)
    def abs_max(a): return Fr(max(abs(a.lo), abs(a.hi)), ONE)
    def contains0(a): return a.lo <= 0 <= a.hi
    def widen(a, e):   # e: Fraction >= 0
        k = -((-e.numerator * ONE) // e.denominator)
        return R(a.lo - k, a.hi + k)


def r_sqrt(x):          # x: positive Fraction
    x = Fr(x)
    lo = isqrt(x.numerator * ONE * ONE // x.denominator)
    hi = isqrt(-((-x.numerator * ONE * ONE) // x.denominator)) + 1
    return R(lo, hi)


def atan_series(x):     # x: R with |x| <= 1/2
    s = R(0); p = x; x2 = x * x; k = 0
    while True:
        term = p / (2 * k + 1)
        s = s + term if k % 2 == 0 else s - term
        k += 1
        p = p * x2
        if p.abs_max() < EPS_STOP:
            break
    return s.widen(p.abs_max())          # alternating tail below the next |term|


def r_pi():
    a = atan_series(R.q(Fr(1, 5))); b = atan_series(R.q(Fr(1, 239)))
    return 16 * a - 4 * b


PI = r_pi()


def ln_q(x):            # ln of a positive Fraction, as an R
    x = Fr(x); k = 0
    while x >= 2: x /= 2; k += 1
    while x < 1: x *= 2; k -= 1
    def atanh2(u):      # 2 atanh(u), u in [0, 1/3]
        s = R(0); p = R.q(u); u2 = R.q(u * u); j = 0
        while True:
            s = s + p / (2 * j + 1)
            j += 1
            p = p * u2
            if p.abs_max() < EPS_STOP:
                break
        tail = p.abs_max() / (1 - u * u)
        return (s * 2).widen(2 * tail)
    v = atanh2((x - 1) / (x + 1))
    if k:
        v = v + atanh2(Fr(1, 3)) * k
    return v


def ln_r(x):            # ln of a positive interval: monotone
    a = ln_q(Fr(x.lo, ONE)); b = ln_q(Fr(x.hi, ONE))
    return R(a.lo, b.hi)


def cos_sin(x):         # x: R, any size
    twopi = PI * 2
    k = round(x.mid_fr() / Fr(twopi.lo, ONE))
    y = x - twopi * k
    c = R.q(1); s = R(0); term = R.q(1); n = 0
    while True:
        n += 1
        term = term * y / n
        if n % 4 == 1: s = s + term
        elif n % 4 == 2: c = c - term
        elif n % 4 == 3: s = s - term
        else: c = c + term
        if term.abs_max() < EPS_STOP and n > 4:
            break
    e = term.abs_max() * y.abs_max()     # next term bound (|y| <= 4 here, terms decreasing)
    return c.widen(e), s.widen(e)


class C:
    __slots__ = ("re", "im")

    def __init__(self, re, im=None): self.re = re if isinstance(re, R) else R.q(re); self.im = R(0) if im is None else (im if isinstance(im, R) else R.q(im))
    def __add__(a, b): b = b if isinstance(b, C) else C(b); return C(a.re + b.re, a.im + b.im)
    def __sub__(a, b): b = b if isinstance(b, C) else C(b); return C(a.re - b.re, a.im - b.im)
    def __mul__(a, b):
        if isinstance(b, C):
            return C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
        return C(a.re * b, a.im * b)
    def conj(a): return C(a.re, -a.im)
    def norm2(a): return a.re * a.re + a.im * a.im
    def inv(a):
        n = a.norm2()
        return C(a.re / n, -a.im / n)
    def __truediv__(a, b):
        if isinstance(b, C): return a * b.inv()
        return C(a.re / b, a.im / b)
    def abs_max(a): return a.re.abs_max() + a.im.abs_max()


def bernoulli(m):
    """B_0..B_m by sum_{j<=k} C(k+1, j) B_j = 0"""
    B = [Fr(1)]
    from math import comb
    for k in range(1, m + 1):
        B.append(-sum(comb(k + 1, j) * B[j] for j in range(k)) / (k + 1))
    return B


BN = bernoulli(64)


def zeta_em(t, N=40, p=30):
    """zeta(1/2 + it) by Euler-Maclaurin; t a Fraction"""
    sigma = Fr(1, 2)
    s = C(R.q(sigma), R.q(t))
    acc = C(0)
    for n in range(1, N):
        c, si = cos_sin(ln_q(n) * t)
        r = r_sqrt(Fr(1, n))
        acc = acc + C(c * r, -(si * r))          # n^-s = n^-1/2 (cos(t ln n) - i sin(t ln n))
    cN, sN = cos_sin(ln_q(N) * t)
    Nms = C(cN * r_sqrt(Fr(1, N)), -(sN * r_sqrt(Fr(1, N))))   # N^-s
    acc = acc + Nms * N / (s - 1)                # N^(1-s)/(s-1)
    acc = acc + Nms / 2
    # s(s+1)...(s+2k-2) is kept EXACT (t is rational), multiplied by the exact coefficient, and only then rounded:
    # rounding the tiny coefficients first would multiply one rounding unit by |s(s+1)...|, which reaches 10^87
    from math import factorial
    pr, pi_ = sigma, Fr(t)                       # the exact complex number s
    for k in range(1, p + 1):
        if k > 1:
            for add in (2 * k - 3, 2 * k - 2):
                ar, ai = sigma + add, Fr(t)
                pr, pi_ = pr * ar - pi_ * ai, pr * ai + pi_ * ar
        coef = BN[2 * k] / factorial(2 * k) / Fr(N) ** (2 * k - 1)
        acc = acc + C(R.q(pr * coef), R.q(pi_ * coef)) * Nms
    # remainder: |B_2p|/(2p)! |s...(s+2p-1)| N^(1-sigma-2p)/(sigma+2p-1)
    mod = Fr(1)
    for j in range(2 * p):
        mod *= Fr(isqrt(int((Fr(1, 2) + j) ** 2 * 10 ** 40 + t * t * 10 ** 40)) + 1, 10 ** 20)
    rem = abs(BN[2 * p]) / factorial(2 * p) * mod / (Fr(N) ** (2 * p - 1) * Fr(isqrt(N * 10 ** 40), 10 ** 20)) / (sigma + 2 * p - 1)
    return C(acc.re.widen(rem), acc.im.widen(rem)), rem


def im_lngamma(z_re, z_im, K=20, m=20):
    """Im log Gamma(z) (a branch), z = z_re + i z_im with z_re > 0 rational"""
    # shift: log Gamma(z) = log Gamma(z + K) - sum_{j<K} log(z + j); Im log(x + iy) = atan(y/x) for x > 0
    def arg(x, y):    # x > 0, y > 0 rational: atan(y/x) = 2 atan(y/(x + sqrt(x^2+y^2))), halved until < 1/4
        r = r_sqrt(x * x + y * y)
        q = R.q(y) / (R.q(x) + r)
        mult = 2
        while q.abs_max() >= Fr(1, 4):
            # atan(q) = 2 atan(q / (1 + sqrt(1 + q^2))), q as an interval: monotone in q, use the end points
            def half(v):
                return R.q(v) / (R.q(1) + r_sqrt(1 + v * v))
            a = half(Fr(q.lo, ONE)); b = half(Fr(q.hi, ONE))
            q = R(a.lo, b.hi)
            mult *= 2
        return atan_series(q) * mult
    shift = R(0)
    for j in range(K):
        shift = shift + arg(z_re + j, z_im)
    wx, wy = z_re + K, z_im
    lnabs = ln_q(wx * wx + wy * wy) / 2
    argw = arg(wx, wy)
    # Im[(w - 1/2) log w - w] = (wx - 1/2) arg w + wy ln|w| - wy
    v = argw * (wx - Fr(1, 2)) + lnabs * wy - R.q(wy)
    w = C(R.q(wx), R.q(wy))
    winv = w.inv()
    wpow = winv                                 # w^-(2k-1)
    w2inv = winv * winv
    for k in range(1, m + 1):
        v = v + (wpow * R.q(BN[2 * k] / (2 * k * (2 * k - 1)))).im
        wpow = wpow * w2inv
    absw = Fr(isqrt(int((wx * wx + wy * wy) * 10 ** 40)), 10 ** 20)
    # sec(arg w / 2)^(2m+2) <= (1/cos(1/2 * 0.4))^(2m+2) since arg w < 0.4 here; cos(0.2) > 0.98
    assert wy / wx < Fr(2, 5)
    rem = abs(BN[2 * m + 2]) / ((2 * m + 2) * (2 * m + 1)) / absw ** (2 * m + 1) / Fr(98, 100) ** (2 * m + 2)
    return (v - shift).widen(rem)


def Z(t, ctrl=None):
    """Hardy's Z(t) = Re(e^{i theta} zeta(1/2+it)); returns (Z interval, Im part interval, zeta remainder)"""
    th = im_lngamma(Fr(1, 4), t / 2) - ln_q(Fr(PI.lo, ONE)).widen(PI.width() * 4) * (t / 2)
    if ctrl == "theta-sign":
        th = -th
    c, s = cos_sin(th)
    zt, rem = zeta_em(t)
    e = C(c, s)
    v = e * zt
    return v.re, v.im, rem


def main():
    a = sys.argv[1:]
    ctrl = a[a.index("--control") + 1] if "--control" in a else None
    printed = "0.755106"
    t_lo = Fr("14.13472514173469378")
    t_hi = Fr("14.13472514173469380")
    if ctrl == "bracket":
        t_lo, t_hi = Fr("14.13472514173469381"), Fr("14.13472514173469383")
    print("program A (Euler-Maclaurin N = 40, p = 30; Stirling shift 20, 20 terms; fixed point 2^-%d)%s" %
          (P, (" CONTROL " + ctrl) if ctrl else ""))
    # self-tests: pi, ln 2, cos^2 + sin^2, and zeta(1/2 + it) Im(e^{i theta} zeta) ~ 0 at t = 10
    assert Fr(314159265358979323846, 10 ** 20) < Fr(PI.lo, ONE) and Fr(PI.hi, ONE) < Fr(314159265358979323847, 10 ** 20)
    c, s = cos_sin(R.q(Fr(37, 3)))
    one = c * c + s * s
    assert one.lo <= ONE <= one.hi
    zr, zi, rem = Z(Fr(10), ctrl)
    print("self-test at t = 10: Z = %.15f, Im part contains 0: %s, EM remainder %.1e" %
          (float(zr.mid_fr()), zi.contains0(), float(rem)))
    ok = zi.contains0()
    out = []
    for t in (t_lo, t_hi):
        zr, zi, rem = Z(t, ctrl)
        ok = ok and zi.contains0()
        out.append(zr)
        print("Z(%s) in [%.6e, %.6e] (width %.1e), Im part contains 0: %s" %
              (t, float(Fr(zr.lo, ONE)), float(Fr(zr.hi, ONE)), float(zr.width()), zi.contains0()))
    change = (out[0].hi < 0 < out[1].lo) or (out[0].lo > 0 > out[1].hi)
    print("sign change of Z on [%s, %s]: %s" % (t_lo, t_hi, change))
    if change:
        L_hi = ln_q(t_hi); L_lo = ln_q(t_lo)
        b_lo = R.q(2) / L_hi                     # 2/log(t_hi): a lower bound for 2/log gamma
        b_hi = R.q(2) / L_lo
        lo, hi = Fr(b_lo.lo, ONE), Fr(b_hi.hi, ONE)
        print("so zeta has a zero 1/2 + i gamma with gamma in (%s, %s), and R >= 2/log gamma in [%s, %s]" %
              (t_lo, t_hi, dec(lo, 25), dec(hi, 25)))
        d = Fr(printed); u = Fr(1, 10 ** 6)
        tr = d <= lo and hi < d + u
        rd = abs(lo - d) <= u / 2 and abs(hi - d) <= u / 2
        rel = ", ".join(x for x, y in (("truncation", tr), ("correct rounding", rd)) if y) or "NONE"
        print("row 8a (lower bound, no reference): 2/log gamma_1 = %s...; page %s: %s; printed <= 2/log gamma: %s" %
              (dec(lo, 20), printed, rel, d <= lo))
        good = ok and rd and d <= lo
    else:
        good = False
    print("VERDICT " + (("REFUSED" if not good else "NOT REFUSED") if ctrl else ("PASS" if good else "FAILED")))


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


if __name__ == "__main__":
    main()
