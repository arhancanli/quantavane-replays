#!/usr/bin/env python3
"""8a, pull request #373: the finite inputs of the deduction C_8 <= 2/log(gamma_1) ~ 0.755107*, program A.

The row (OpenAI's 7/8 quasi-Riemann hypothesis [OAI2026] with Platt and Trudgian's verification of RH to 3e12) needs:
  (1) a zero 1/2 + i gamma_1 of zeta with gamma_1 known closely: Hardy's Z(t) changes sign on
      (14.13472514173469378, 14.13472514173469380) (the method of our folder c8a);
  (2) that this is the lowest zero: there is no zero beta + i gamma with 0 < gamma < gamma_1 (so among zeros with
      |gamma| >= 2 the bound 2/log|gamma| is largest at gamma_1). We show N(T0) = 1 at T0 = 89/5 = 17.8 through the
      Riemann-von Mangoldt identity N(T) = theta(T)/pi + 1 + S(T) (read, not checked): if Re zeta(sigma + i T0) > 0 for
      1/2 <= sigma <= 2 (and Re zeta(2 + it) > 0 always, since |zeta(2 + it) - 1| <= zeta(2) - 1 < 1), then
      |S(T0)| < 1/2, and N(T0) is the only integer within 1/2 of theta(T0)/pi + 1. With N(T0) = 1 and the zero of (1)
      below T0, gamma_1 is the lowest ordinate and its zero is on the critical line;
  (3) the arithmetic: 2/log gamma_1 <= 0.755107 (the printed bound), gamma_1^4 < 4 * 10^4 (the description's remark),
      and 8/log(3 * 10^12) < 0.279 < 2/log gamma_1 (the part above Platt-Trudgian's height).
Not checked: [OAI2026]'s theorem (zeta(s) != 0 for Re s > 7/8; unrefereed, 199 pages) and [PT2021]'s verification
of RH to 3e12 (refereed); both enter the deduction as stated in the pull request.

Re zeta(sigma + i T0) on [1/2, 2] is enclosed on subintervals of sigma (sigma as an interval: n^-sigma =
exp(-sigma ln n) is monotone in sigma) by Euler-Maclaurin with N = 40 and 30 Bernoulli terms; the remainder bound
|B_2p|/(2p)! |s(s+1)...(s+2p-1)| N^(1-sigma-2p)/(sigma+2p-1) is taken at the worst end of each subinterval. Intervals
are bisected until the lower end of Re zeta is positive. theta(T0) = Im log Gamma(1/4 + i T0/2) - (T0/2) log pi by
Stirling's series (shift 30, 20 terms; shift 20 for Z in (1)). Arithmetic: closed intervals of integers at scale 2^-320 with outward rounding,
our own ln, exp, sqrt, cos, sin, atan; pi by Machin. Standard library only.
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




def exp_q(x):
    """exp of a Fraction x with |x| <= 64, as an R: exp(x/2^m)^(2^m) with |x/2^m| <= 1/16, Taylor with remainder"""
    x = Fr(x); m = 0
    while abs(x) > Fr(1, 16):
        x /= 2; m += 1
    s = R.q(1); term = R.q(1); k = 0; xr = R.q(x)
    while True:
        k += 1
        term = term * xr / k
        s = s + term
        if term.abs_max() < EPS_STOP:
            break
    s = s.widen(2 * term.abs_max())          # tail of a series with |x| <= 1/16 is below twice the last term
    for _ in range(m):
        s = s * s
    return s


def npow_neg(n, slo, shi, lncache={}):
    """n^-sigma for sigma in [slo, shi] (rationals), n >= 1: decreasing in sigma"""
    if n == 1:
        return R.q(1)
    if n not in lncache:
        lncache[n] = ln_q(n)
    L = lncache[n]
    a = exp_q(-shi * Fr(L.hi, ONE))           # smallest: largest exponent magnitude
    b = exp_q(-slo * Fr(L.lo, ONE))
    return R(a.lo, b.hi)


def zeta_strip(slo, shi, t, N=40, p=30, cs_cache={}):
    """zeta(sigma + it) for sigma in [slo, shi] (1/2 <= slo <= shi), t a Fraction: Euler-Maclaurin, a complex interval"""
    from math import factorial
    sig = R(R.q(slo).lo, R.q(shi).hi)
    s = C(sig, R.q(t))
    acc = C(0)
    for n in range(1, N + 1):
        key = (n, t)
        if key not in cs_cache:
            cs_cache[key] = cos_sin(ln_q(n) * t) if n > 1 else (R.q(1), R(0))
        c, si = cs_cache[key]
        r = npow_neg(n, slo, shi)
        term = C(c * r, -(si * r))            # n^-s
        if n < N:
            acc = acc + term
        else:
            Nms = term
    acc = acc + Nms * N / (s - 1)              # N^(1-s)/(s-1)
    acc = acc + Nms / 2
    # Bernoulli terms: B_2k/(2k)! s(s+1)...(s+2k-2) N^(1-s-2k): the product as a complex interval, multiplied by the
    # integer numerator of B_2k/((2k)! N^(2k-1)) and then divided by its integer denominator (one rounding per step)
    prod = C(sig, R.q(t))
    for k in range(1, p + 1):
        if k > 1:
            for add in (2 * k - 3, 2 * k - 2):
                prod = prod * C(sig + add, R.q(t))
        coef = BN[2 * k] / factorial(2 * k) / Fr(N) ** (2 * k - 1)
        num, den = coef.numerator, coef.denominator
        tk = prod * Nms
        tk = C(tk.re * num, tk.im * num)
        tk = C(tk.re / den, tk.im / den) if den > 0 else tk
        acc = acc + tk
    # remainder at the worst end: |s+j| largest at shi, N^(1-sigma-2p)/(sigma+2p-1) largest at slo
    mod = Fr(1)
    for j in range(2 * p):
        mod *= Fr(isqrt(int((shi + j) ** 2 * 10 ** 40 + t * t * 10 ** 40)) + 1, 10 ** 20)
    npw = exp_q(-(slo + 2 * p - 1) * Fr(ln_q(N).lo, ONE))          # N^-(sigma + 2p - 1) <= exp(-(slo+2p-1) ln N)
    rem = abs(BN[2 * p]) / factorial(2 * p) * mod * Fr(npw.hi, ONE) / (slo + 2 * p - 1)
    return C(acc.re.widen(rem), acc.im.widen(rem)), rem


def theta(t):
    return im_lngamma(Fr(1, 4), t / 2, K=30) - ln_q(Fr(PI.lo, ONE)).widen(PI.width() * 4) * (t / 2)


def re_positive_on_segment(t, slo=Fr(1, 2), shi=Fr(2), pieces=24, maxdepth=14, flip=False):
    """bisect [slo, shi] until every piece has Re zeta > 0 (or report the first failing piece)"""
    stack = [(slo + (shi - slo) * k / pieces, slo + (shi - slo) * (k + 1) / pieces, 0) for k in range(pieces)]
    done = 0; worst = None
    while stack:
        a, b, d = stack.pop()
        z, rem = zeta_strip(a, b, t)
        re = -z.re if flip else z.re
        if re.lo > 0:
            done += 1
            lo = Fr(re.lo, ONE)
            worst = lo if worst is None or lo < worst else worst
            continue
        if re.hi < 0 or d >= maxdepth:
            return False, (a, b, float(Fr(re.lo, ONE)), float(Fr(re.hi, ONE))), done
        m = (a + b) / 2
        stack += [(a, m, d + 1), (m, b, d + 1)]
    return True, worst, done


def main():
    a = sys.argv[1:]
    ctrl = a[a.index("--control") + 1] if "--control" in a else None
    T0 = Fr(89, 5)
    claim_N = 1
    printed = "0.755107"
    t_lo, t_hi = Fr("14.13472514173469378"), Fr("14.13472514173469380")
    if ctrl == "height22":
        T0 = Fr(22)                              # above the second zero, 21.022...: N(22) = 2, the claim N = 1 must fail
    if ctrl == "printed":
        printed = "0.755106"                     # below 2/log gamma_1: refused as an upper bound
    if ctrl == "bracket":
        t_lo, t_hi = Fr("14.13472514173469381"), Fr("14.13472514173469383")   # just above gamma_1: no sign change
    print("program A (Euler-Maclaurin N = 40, p = 30 with sigma as an interval; Stirling shift 20; fixed point 2^-%d)%s" %
          (P, (" CONTROL " + ctrl) if ctrl else ""))
    assert Fr(314159265358979323846, 10 ** 20) < Fr(PI.lo, ONE) and Fr(PI.hi, ONE) < Fr(314159265358979323847, 10 ** 20)
    e1 = exp_q(Fr(1))
    assert Fr(27182818284590452353, 10 ** 19) < Fr(e1.lo, ONE) and Fr(e1.hi, ONE) < Fr(27182818284590452354, 10 ** 19)
    # self-test: the interval Euler-Maclaurin at a thin sigma = 2, t = 0 would be zeta(2); use t = 1/1000 instead
    z2, _ = zeta_strip(Fr(2), Fr(2), Fr(1, 10 ** 30))
    pi2_6 = PI * PI / 6
    assert z2.re.lo - 2 ** 200 <= pi2_6.hi and pi2_6.lo <= z2.re.hi + 2 ** 200, "zeta(2) self-test"
    print("self-tests: pi, e, zeta(2 + 10^-30 i) = pi^2/6 within 2^-120: ok")
    ok = True
    # (1) the zero on the line (as in folder c8a)
    vals = []
    for t in (t_lo, t_hi):
        zr, zi, rem = Z(t)
        ok = ok and zi.contains0()
        vals.append(zr)
    change = (vals[0].hi < 0 < vals[1].lo) or (vals[0].lo > 0 > vals[1].hi)
    print("(1) Z changes sign on (%s, %s): %s, so zeta(1/2 + i gamma_1) = 0 there" % (t_lo, t_hi, change))
    ok = ok and change
    # (2) Re zeta > 0 on the segment at height T0, then N(T0)
    good, info, pieces = re_positive_on_segment(T0)
    if good:
        print("(2) Re zeta(sigma + %s i) > 0 for 1/2 <= sigma <= 2: yes (%d subintervals, smallest lower end %.6f), so |S(T0)| < 1/2"
              % (T0, pieces, float(info)))
    else:
        print("(2) Re zeta(sigma + %s i) > 0 for 1/2 <= sigma <= 2: NOT SHOWN (piece [%s, %s]: Re in [%.3e, %.3e])"
              % (T0, info[0], info[1], info[2], info[3]))
    th = theta(T0)
    x = th / PI
    lo, hi = Fr(x.lo, ONE) + 1, Fr(x.hi, ONE) + 1
    cands = [n for n in range(-2, 10) if lo - Fr(1, 2) < n < hi + Fr(1, 2)]
    print("    theta(%s)/pi + 1 in [%.12f, %.12f]; integers within 1/2: %s" % (T0, float(lo), float(hi), cands))
    nT = cands[0] if good and len(cands) == 1 else None
    print("    N(%s) = %s; the claim N = %d: %s" % (T0, nT, claim_N, "holds" if nT == claim_N else "FAILS"))
    ok = ok and good and nT == claim_N
    # (3) arithmetic
    L_lo, L_hi = ln_q(t_lo), ln_q(t_hi)
    b_lo = Fr((R.q(2) / L_hi).lo, ONE); b_hi = Fr((R.q(2) / L_lo).hi, ONE)
    d = Fr(printed)
    print("(3) 2/log gamma_1 in [%s, %s]; printed upper bound %s: %s" % (dec(b_lo, 22), dec(b_hi, 22), printed,
          "holds" if b_hi <= d else "FAILS"))
    ok = ok and b_hi <= d
    g4 = t_hi ** 4
    print("    gamma_1^4 < %s^4 = %.4f < 40000: %s" % (t_hi, float(g4), g4 < 40000))
    ok = ok and g4 < 40000
    l3 = ln_q(3 * 10 ** 12)
    e8 = Fr((R.q(8) / l3).hi, ONE)
    print("    8/log(3*10^12) <= %s < 0.279 < 2/log gamma_1: %s" % (dec(e8, 12), e8 < Fr("0.279") < b_lo))
    ok = ok and e8 < Fr("0.279") < b_lo
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if ctrl else ("PASS" if ok else "FAILED")))


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


if __name__ == "__main__":
    main()
