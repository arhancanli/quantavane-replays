#!/usr/bin/env python3
"""Independent replay of the asterisked lower bound C_70 >= 1.20194336* ([Pin2026], R. Pinho; PR #259 of
teorth/optimizationproblems, constant 70a) from the certificate's data file.

Written from the 70a page, the paper (paper/c70.tex at 680e074) and the certificate's README. The claimant's
programs (verify_c70_arb.py, independent/verify_c70_rational.py, reproduce.sh) were neither read nor run.
Standard library only (Python 3.9+); run with python3 -I. No floating-point number decides anything: every
decision is a comparison of exact integers, exact rationals, or the ends of rigorous intervals.

  c70_check.py CERTS_ROOT [--controls]

CERTS_ROOT holds pin/certificates/c70_candidate.json (rpinho/c70-certificate @680e074), as fetch.sh downloads it.

What is computed
  Patankar's criterion (70a page; paper, Lemma 1): for the cube K = [-a/2, a/2]^n and the unit ball B,
    C70 >= (sum_k C(n,k) a^(n-k) kappa_k)^(1/n) / (a + kappa_n^(1/n)),   kappa_j = pi^(j/2) / Gamma(1 + j/2).
  Keeping only the term k, with a = s kappa_n^(1/n), gives C70 >= R_{n,k}(s),
    log R(s) = log C(n,k)/n + (1-t) log s - log(1+s) + (log kappa_k - t log kappa_n)/n,   t = k/n,
  and at s = (n-k)/k, for even n and k (paper, (2)), with L(m) = log m!,
    log R = [L(n) - L(k) - L(n-k) - L(k/2) + t L(n/2)]/n + (1-t) log(1-t) + t log t    (pi cancels).
  The certificate: n = 10^12, k = 367879441172, s = 632120558828/367879441172; the page's claim 1.20194336.
How it is decided
  Fixed-point intervals at scale 2^-400 (integers lo <= 2^400 x <= hi). log of an integer N: N = 2^e r with
  1 <= r < 2, log r = 2 atanh((N - 2^e)/(N + 2^e)), every series term rounded outward and the tail bounded;
  log 2 = 2 atanh(1/3). pi from Machin's formula as exact rationals (alternating series: the value lies between
  consecutive partial sums), checked against Euler's atan(1/2) + atan(1/3) and 50 known digits. L(m) by the
  Stirling series with Bernoulli numbers B_2..B_8 (recomputed exactly) and the remainder between 0 and the first
  omitted term (B_10 term), as holds for real m > 0 (DLMF 5.11(ii)); the implementation is tested against
  log m! computed from the exact integer m! for m = 1..60, 500, 3000. exp by outward-rounded Taylor sums.
  Cross-checks: the general-s form with explicit pi; the whole formula at n = 20000 against exact factorials;
  the algebra (criterion -> factorial form) at small n in 70-digit decimal; a 90-digit decimal evaluation with
  pi by the AGM; a double-precision lgamma value. After our run, our enclosures were compared with the
  claimant's published Arb balls (verification/c70_arb.log, an output file).
What rests on mathematics rather than computation (read, not computed): the criterion (Gram-Schmidt / Szasz
principal-minor inequality, AM-GM, intrinsic volumes of a parallelotope, Steiner's formula; Patankar's
reduction, re-proved in the paper), the positivity of the dropped Steiner terms, the Stirling remainder theorem,
and the separate limit C70 >= exp(1/(2e)) (Stirling asymptotics; only its decimal value is enclosed here).
"""
import hashlib
import json
import math
import os
import re
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction

F_CAND = "pin/certificates/c70_candidate.json"
H_CAND = "03ff0d6df3dfa6ea06f901a9319c968eb7c013f5d0cc3f1f00cfa58e9490c0cc"   # SHA256SUMS @680e074

PAGE_CLAIM = "1.20194336"
PAGE_DIGITS = "1.20194336844315056699"                 # 70a page: "enclosed near ..."
README_DIGITS = "1.2019433684431505669948027797355"     # README / paper / FINAL_LIVE_CHECK (Arb)
LIMIT_CLAIM = "1.20194336847"
LIMIT_DIGITS = "1.2019433684703144671942411393813"      # FINAL_LIVE_CHECK (Arb)
PI50 = "3.14159265358979323846264338327950288419716939937510"
# the claimant's Arb balls, verification/c70_arb.log @680e074 (an output file, read after our run): (midpoint, radius)
ARB_VALUE = ("1.20194336844315056699480277973552413784930981371792977194313655532718767431", "5.01e-75")
ARB_LOG = ("0.18393972056312117737731489748866704280577155234802617264248072878198077755", "1.65e-75")
ARB_LIMIT = ("1.201943368470314467194241139381297080440187153935169095630989013831578451122", "3.43e-76")

P = 400
ONE = 1 << P
G = 40                     # guard bits inside series


# ----------------------------------------------------------------------------------------------------------------
# Fixed-point intervals: (lo, hi) integers with lo <= 2^P x <= hi
# ----------------------------------------------------------------------------------------------------------------

def cdiv(a, b):
    return -((-a) // b)


def qiv(fr):
    fr = Fraction(fr)
    n, d = fr.numerator, fr.denominator
    return (n << P) // d, cdiv(n << P, d)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def sub(a, b):
    return a[0] - b[1], a[1] - b[0]


def qmul(iv, fr):
    """iv times the exact rational fr, rounded outward."""
    fr = Fraction(fr)
    p, q = fr.numerator, fr.denominator
    lo, hi = iv
    if p >= 0:
        return lo * p // q, cdiv(hi * p, q)
    return hi * p // q, cdiv(lo * p, q)


def contains(iv, fr):
    fr = Fraction(fr)
    return Fraction(iv[0], ONE) <= fr <= Fraction(iv[1], ONE)


def width(iv):
    return Fraction(iv[1] - iv[0], ONE)


def dec_floor(v, digits, scale=P):
    """floor(v / 2^scale) to `digits` decimals, as a string."""
    q = (v * 10 ** digits) >> scale
    s = str(abs(q)).rjust(digits + 1, "0")
    return ("-" if q < 0 else "") + s[:-digits] + "." + s[-digits:]


def dec_ceil(v, digits, scale=P):
    q = cdiv(v * 10 ** digits, 1 << scale)
    s = str(abs(q)).rjust(digits + 1, "0")
    return ("-" if q < 0 else "") + s[:-digits] + "." + s[-digits:]


def enc(iv, digits=40):
    return "[%s, %s]" % (dec_floor(iv[0], digits), dec_ceil(iv[1], digits))


# ----------------------------------------------------------------------------------------------------------------
# Logarithms, pi, Bernoulli numbers, Stirling, exp
# ----------------------------------------------------------------------------------------------------------------

def atanh_iv(a, b):
    """Interval at scale 2^P for atanh(a/b), integers 0 <= a, 3a <= b."""
    assert 0 <= a and 3 * a <= b
    if a == 0:
        return 0, 0
    Q = P + G
    zl = (a << Q) // b
    zh = cdiv(a << Q, b)
    z2l = (zl * zl) >> Q
    z2h = cdiv(zh * zh, 1 << Q)
    lo, t, k = 0, zl, 1
    while t:
        lo += t // k
        t = (t * z2l) >> Q
        k += 2
    hi, t, k = 0, zh, 1
    while t > 1:
        hi += cdiv(t, k)
        t = cdiv(t * z2h, 1 << Q)
        k += 2
    hi += 2 * t + 2                       # the omitted terms: at most t (1/k)(1/(1 - z^2)) <= 2t, plus rounding
    return lo >> G, cdiv(hi, 1 << G)


_L2 = atanh_iv(1, 3)
LN2 = (2 * _L2[0], 2 * _L2[1])


def ln_int(N):
    """log N for an integer N >= 1."""
    assert N >= 1
    e = N.bit_length() - 1
    z = atanh_iv(N - (1 << e), N + (1 << e))
    return e * LN2[0] + 2 * z[0], e * LN2[1] + 2 * z[1]


def ln_q(fr):
    fr = Fraction(fr)
    return sub(ln_int(fr.numerator), ln_int(fr.denominator))


def ln_iv(iv):
    """log of a positive quantity given as an interval at scale 2^P (log is increasing)."""
    lo = ln_int(iv[0])
    hi = ln_int(iv[1])
    return lo[0] - P * LN2[1], hi[1] - P * LN2[0]


def atan_inv(m, bits):
    """Exact rationals lo <= atan(1/m) <= hi, hi - lo < 2^-bits (alternating series, decreasing terms)."""
    S = Fraction(0)
    i = 0
    m2 = m * m
    pw = m
    while True:
        a = Fraction(1, (2 * i + 1) * pw)
        S2 = S + a if i % 2 == 0 else S - a
        if i >= 1 and a < Fraction(1, 2 ** bits):
            return min(S, S2), max(S, S2)
        S = S2
        i += 1
        pw *= m2


def pi_machin(bits):
    a = atan_inv(5, bits)
    b = atan_inv(239, bits)
    return 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]


def pi_euler(bits):
    a = atan_inv(2, bits)
    b = atan_inv(3, bits)
    return 4 * (a[0] + b[0]), 4 * (a[1] + b[1])


def frac_iv(lo, hi):
    return (lo.numerator << P) // lo.denominator, cdiv(hi.numerator << P, hi.denominator)


PI = frac_iv(*pi_machin(P + 20))
LNPI = ln_iv(PI)
LN2PI = add(LN2, LNPI)


def bernoulli(nmax):
    """B_0..B_nmax exactly, from sum_{k=0}^{m} C(m+1, k) B_k = 0 (B_1 = -1/2)."""
    B = [Fraction(1)]
    for m in range(1, nmax + 1):
        B.append(-sum(math.comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    return B


BERN = bernoulli(12)


def lnfact(m, N=4, flip=False):
    """log m! for an integer m >= 1: (m + 1/2) log m - m + log(2 pi)/2 + sum_{j=1}^{N} B_2j / (2j(2j-1) m^(2j-1))
    + R, with R between 0 and the first omitted term (flip=True reverses that side: a control)."""
    lm = ln_int(m)
    v = qmul(lm, Fraction(2 * m + 1, 2))
    v = (v[0] - m * ONE, v[1] - m * ONE)
    v = add(v, qmul(LN2PI, Fraction(1, 2)))
    for j in range(1, N + 1):
        v = add(v, qiv(BERN[2 * j] / (2 * j * (2 * j - 1) * Fraction(m) ** (2 * j - 1))))
    T = BERN[2 * N + 2] / ((2 * N + 2) * (2 * N + 1) * Fraction(m) ** (2 * N + 1))
    if flip:
        T = -T
    r = qiv(T)
    return v[0] + min(0, r[0]), v[1] + max(0, r[1])


def exp_iv(iv):
    """exp of an interval at scale 2^P (exp is increasing: each end separately, rounded outward)."""
    return _exp_end(iv[0], False), _exp_end(iv[1], True)


def _exp_end(x, up):
    """A lower (up=False) or upper (up=True) bound for 2^P exp(x / 2^P)."""
    if x < 0:
        y = _exp_end(-x, not up)             # exp(x) = 1/exp(-x): a lower end needs an upper exp(-x)
        return cdiv(ONE << P, y) if up else (ONE << P) // y
    r = 0
    while (x >> r) >= ONE // 2:              # reduce to [0, 1/2), then square r times
        r += 1
    xr = cdiv(x, 1 << r) if up else x >> r
    y = _exp_small(xr, up)
    for _ in range(r):
        y = cdiv(y * y, ONE) if up else (y * y) >> P
    return y


def _exp_small(x, up):
    assert 0 <= x < ONE
    Q = P + G
    xx = x << G
    s, t, i = 1 << Q, 1 << Q, 1
    while True:
        t = cdiv(t * xx, i << Q) if up else (t * xx) // (i << Q)
        if t <= 1:
            break
        s += t
        i += 1
    if up:
        s += 2 * t + 2                        # omitted terms at most t / (1 - x/(i+1)) <= 2t, plus rounding
        return cdiv(s, 1 << G)
    return s >> G


def e_bounds(N=80):
    """Exact rationals lo < e < hi: sum_{i<=N} 1/i! < e < that + 1/(N! N)."""
    s = Fraction(0)
    f = 1
    for i in range(N + 1):
        if i:
            f *= i
        s += Fraction(1, f)
    return s, s + Fraction(1, f * N)


# ----------------------------------------------------------------------------------------------------------------
# The retained Steiner term
# ----------------------------------------------------------------------------------------------------------------

def log_R(n, k, lf=lnfact):
    """log R_{n,k} at s = (n-k)/k by the factorial form (even n, k)."""
    if n % 2 or k % 2 or not (0 < k < n):
        raise ValueError("the factorial form needs even n and k with 0 < k < n (got n = %d, k = %d)" % (n, k))
    t = Fraction(k, n)
    num = sub(sub(sub(lf(n), lf(k)), lf(n - k)), lf(k // 2))
    num = add(num, qmul(lf(n // 2), t))
    first = qmul(num, Fraction(1, n))
    lt = sub(ln_int(k), ln_int(n))
    l1t = sub(ln_int(n - k), ln_int(n))
    return add(first, add(qmul(l1t, 1 - t), qmul(lt, t)))


def ln_kappa(j, lf=lnfact):
    """log kappa_j = (j/2) log pi - log((j/2)!) for even j >= 2."""
    if j % 2:
        raise ValueError("odd j needs Gamma at a half-integer")
    return sub(qmul(LNPI, Fraction(j, 2)), lf(j // 2))


def log_R_general(n, k, s, lf=lnfact):
    """log R_{n,k}(s) = log C(n,k)/n + (1-t) log s - log(1+s) + (log kappa_k - t log kappa_n)/n, explicit pi."""
    t = Fraction(k, n)
    s = Fraction(s)
    binom = sub(sub(lf(n), lf(k)), lf(n - k))
    kap = sub(ln_kappa(k, lf), qmul(ln_kappa(n, lf), t))
    v = qmul(add(binom, kap), Fraction(1, n))
    return add(v, sub(qmul(ln_q(s), 1 - t), ln_q(1 + s)))


def decimal_pi():
    """pi by the Gauss-Legendre AGM in the current decimal context (not rigorous; a cross-check)."""
    a, b, t, p = Decimal(1), Decimal(1) / Decimal(2).sqrt(), Decimal(1) / 4, Decimal(1)
    for _ in range(12):
        an = (a + b) / 2
        b = (a * b).sqrt()
        t -= p * (a - an) ** 2
        a = an
        p *= 2
    return (a + b) ** 2 / (4 * t)


def decimal_R(n, k):
    """R_{n,k} by the factorial form in decimal arithmetic, Stirling with 6 Bernoulli terms (a cross-check)."""
    pi = decimal_pi()
    half_l2pi = (2 * pi).ln() / 2

    def lf(m):
        M = Decimal(m)
        v = (M + Decimal(1) / 2) * M.ln() - M + half_l2pi
        for j in range(1, 7):
            b = BERN[2 * j]
            v += Decimal(b.numerator) / Decimal(b.denominator) / (2 * j * (2 * j - 1)) / M ** (2 * j - 1)
        return v
    t = Decimal(k) / Decimal(n)
    x = (lf(n) - lf(k) - lf(n - k) - lf(k // 2) + t * lf(n // 2)) / n + (1 - t) * (1 - t).ln() + t * t.ln()
    return x.exp()


def lnfact_exact(m):
    return ln_int(math.factorial(m)) if m > 1 else (0, 0)


# ----------------------------------------------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------------------------------------------

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Report:
    def __init__(self, title):
        self.fails = []
        self.t0 = time.time()
        print(title, flush=True)

    def info(self, key, val):
        print("  %s: %s" % (key, val), flush=True)

    def check(self, desc, ok, detail=""):
        print("  check %s: %s%s" % (desc, "PASS" if ok else "FAIL", (" (" + detail + ")") if detail else ""),
              flush=True)
        if not ok:
            self.fails.append(desc)
        return ok

    def done(self):
        print("  run time: %.1f s" % (time.time() - self.t0))
        print("  RESULT: %s" % ("ALL CHECKS PASS" if not self.fails else "FAILED: " + "; ".join(self.fails)))
        return not self.fails


def digits_agree(iv, s):
    """Both ends of iv truncate to the decimal string s (so every point of iv begins with s)."""
    d = len(s.split(".")[1])
    return dec_floor(iv[0], d) == s and dec_floor(iv[1], d) == s


def digits_round(iv, s):
    """Every point of iv rounds to the decimal string s at its number of decimals (|x - s| < half a unit)."""
    d = len(s.split(".")[1])
    half = Fraction(1, 2 * 10 ** d)
    return Fraction(s) - half <= Fraction(iv[0], ONE) and Fraction(iv[1], ONE) < Fraction(s) + half


def inside_ball(iv, ball):
    """Our enclosure lies inside the ball [mid +- rad]."""
    mid, rad = Fraction(ball[0]), Fraction(ball[1])
    return mid - rad <= Fraction(iv[0], ONE) and Fraction(iv[1], ONE) <= mid + rad


def selftests(rep):
    # logarithms against the correctly rounded decimal ln
    getcontext().prec = 140
    worst = 0
    for N in (2, 3, 7, 10, 255, 256, 257, 10 ** 12, 367879441172, 632120558828, 2 ** 89 - 1, 3 ** 70):
        iv = ln_int(N)
        ref = Decimal(N).ln()
        lo = Decimal(iv[0]) / Decimal(2) ** P
        hi = Decimal(iv[1]) / Decimal(2) ** P
        ok = lo - Decimal(10) ** -118 <= ref <= hi + Decimal(10) ** -118
        worst = max(worst, iv[1] - iv[0])
        if not ok:
            rep.check("log %d contains the decimal ln" % N, False)
            return
    rep.check("rigorous log of 12 integers (2 ... 3^70) contains the correctly rounded decimal ln (to 1e-118)", True,
              "widest enclosure %.1e" % float(Fraction(worst, ONE)))
    pm = pi_machin(P + 20)
    pe = pi_euler(P + 20)
    rep.check("pi: Machin and Euler enclosures overlap, each narrower than 2^-400, and both begin with the 50 known "
              "digits", max(pm[0], pe[0]) <= min(pm[1], pe[1]) and pm[1] - pm[0] < Fraction(1, 2 ** 400) and
              pe[1] - pe[0] < Fraction(1, 2 ** 400) and
              all(Fraction(PI50) <= x < Fraction(PI50) + Fraction(1, 10 ** 50) for x in (pm[0], pe[1])))
    rep.check("Bernoulli numbers B_2..B_12 recomputed: 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730",
              [BERN[i] for i in (2, 4, 6, 8, 10, 12)] == [Fraction(1, 6), Fraction(-1, 30), Fraction(1, 42),
                                                           Fraction(-1, 30), Fraction(5, 66), Fraction(-691, 2730)])
    bad = []
    for m in list(range(1, 61)) + [500, 3000]:
        ex = lnfact_exact(m)
        st = lnfact(m)
        if not (st[0] <= ex[0] and ex[1] <= st[1]):
            bad.append(m)
    rep.check("Stirling enclosure (4 Bernoulli terms, remainder between 0 and the B_10 term) contains log m! from "
              "the exact integer m! for m = 1..60, 500, 3000", not bad, "fails at m = %s" % bad if bad else "")
    flipped = [m for m in range(1, 61) if not (lnfact(m, flip=True)[0] <= lnfact_exact(m)[0] and
                                                lnfact_exact(m)[1] <= lnfact(m, flip=True)[1])]
    rep.check("control: the same test catches a remainder put on the wrong side (fails at %d of 60 m)" % len(flipped),
              len(flipped) > 0)
    x = qiv(Fraction(1, 7))
    ex = exp_iv(x)
    getcontext().prec = 140
    ref = (Decimal(1) / 7).exp()
    ok = Decimal(ex[0]) / Decimal(2) ** P <= ref <= Decimal(ex[1]) / Decimal(2) ** P
    for v in (Fraction(-1, 7), Fraction(-27, 10), Fraction(31, 10)):
        e2 = exp_iv(qiv(v))
        ref = (Decimal(v.numerator) / v.denominator).exp()
        ok = ok and Decimal(e2[0]) / Decimal(2) ** P <= ref <= Decimal(e2[1]) / Decimal(2) ** P and \
            e2[1] - e2[0] < 2 ** 30
    rep.check("rigorous exp at 1/7, -1/7, -2.7, 3.1 contains the correctly rounded decimal exp", ok)


def algebra_check(rep):
    """The criterion's single term evaluated directly (60-digit decimal, kappa from pi and factorials) against the
    factorial form, at small even n, k; and the full criterion sum >= the term."""
    getcontext().prec = 70
    pi = Decimal(PI[0]) / Decimal(2) ** P

    def kappas(nmax):       # kappa_j = 2 pi kappa_{j-2} / j, kappa_0 = 1, kappa_1 = 2
        ks = [Decimal(1), Decimal(2)]
        for i in range(2, nmax + 1):
            ks.append(2 * pi * ks[i - 2] / i)
        return ks
    worst = Decimal(0)
    sums_ok = True
    for n, k in ((4, 2), (6, 2), (6, 4), (10, 4), (20, 8), (100, 36), (400, 148)):
        kap = kappas(n)
        kn = kap[n]
        kk = kap[k]
        root = (kn.ln() / n).exp()
        s = Decimal(n - k) / Decimal(k)
        a = s * root
        term = Decimal(math.comb(n, k)) * a ** (n - k) * kk
        direct = (term.ln() / n).exp() / (a + root)
        t = Decimal(k) / Decimal(n)
        lf = [Decimal(math.factorial(m)).ln() for m in (n, k, n - k, k // 2, n // 2)]
        form = ((lf[0] - lf[1] - lf[2] - lf[3] + t * lf[4]) / n + (1 - t) * (1 - t).ln() + t * t.ln()).exp()
        worst = max(worst, abs(direct - form))
        total = sum(Decimal(math.comb(n, j)) * a ** (n - j) * kap[j] for j in range(n + 1))
        full = (total.ln() / n).exp() / (a + root)
        sums_ok = sums_ok and full > direct
    rep.check("algebra: the criterion's single term at s = (n-k)/k equals the factorial form (pi cancels) at 7 small "
              "(n, k), to 1e-55", worst < Decimal(10) ** -55, "largest difference %.1e" % worst)
    rep.check("the full criterion sum exceeds the single retained term at those (n, k) (the dropped terms are "
              "positive)", sums_ok)


def run(root, do_controls):
    rep = Report("ROW [Pin2026] (R. Pinho, PR #259): C70 >= 1.20194336*, rpinho/c70-certificate @680e074")
    path = os.path.join(root, F_CAND)
    h = sha256(path)
    rep.info("sha256 " + F_CAND, "%s (%d bytes)" % (h, os.path.getsize(path)))
    rep.check("sha256 matches SHA256SUMS @680e074 (the 70a page states the commit only)", h == H_CAND)
    cand = json.load(open(path))
    rep.info("c70_candidate.json", json.dumps(cand, sort_keys=True))
    n, k = cand["n"], cand["k"]
    assert isinstance(n, int) and isinstance(k, int) and not isinstance(n, bool) and not isinstance(k, bool)
    s = Fraction(cand["normalized_side"])
    claim = Fraction(cand["claim"])
    rep.check("certificate: n = 10^12, k = 367879441172, both even, 0 < k < n (as on the 70a page)",
              n == 10 ** 12 and k == 367879441172 and n % 2 == 0 and k % 2 == 0 and 0 < k < n)
    rep.check("normalized side a/kappa_n^(1/n) = 632120558828/367879441172 = (n-k)/k exactly, as on the page",
              cand["normalized_side"] == "632120558828/367879441172" and s == Fraction(n - k, k))
    rep.check("claim 1.20194336 = 3756073/3125000 (the page's bound, the paper's Theorem 1)",
              cand["claim"] == PAGE_CLAIM and claim == Fraction(3756073, 3125000))
    rep.check("limit claim 1.20194336847 and formula exp(1/(2e)) as on the page",
              cand.get("limit_claim") == LIMIT_CLAIM and cand.get("limit_formula") == "exp(1/(2*e))")
    rep.info("t = k/n", "%s = %s (1/e = 0.36787944117144233...)" % (Fraction(k, n), float(Fraction(k, n))))
    selftests(rep)
    algebra_check(rep)
    # ---- the finite certificate
    t0 = time.time()
    L = log_R(n, k)
    R = exp_iv(L)
    rep.info("log R_{n,k}", enc(L, 45))
    rep.info("R_{n,k} enclosure", enc(R, 45))
    rep.info("enclosure width", "%.2e (computed in %.2f s)" % (float(width(R)), time.time() - t0))
    lc = ln_q(claim)
    rep.check("log R > log(1.20194336): the claim holds (C70 >= R > 1.20194336)", L[0] > lc[1],
              "margin in R: %s" % dec_floor(R[0] - qiv(claim)[1], 30))
    rep.check("R > 1.20194336 also compared directly after exp", R[0] > qiv(claim)[1])
    rep.check("every point of the enclosure begins 1.20194336844315056699 and rounds to it (the 70a page's digits)",
              digits_agree(R, PAGE_DIGITS) and digits_round(R, PAGE_DIGITS))
    rep.check("every point of the enclosure begins 1.2019433684431505669948027797355 and rounds to it (the README's "
              "and the claimant's Arb digits, 31 decimals)", digits_agree(R, README_DIGITS) and
              digits_round(R, README_DIGITS))
    rep.info("R_{n,k} to 75 decimals", enc(R, 75))
    rep.check("our enclosures of log R and R lie inside the claimant's 256-bit Arb balls (c70_arb.log: radius "
              "1.65e-75 and 5.01e-75): the same value, 117-digit-narrow here", inside_ball(L, ARB_LOG) and
              inside_ball(R, ARB_VALUE))
    G2 = log_R_general(n, k, s)
    rep.check("cross-check: the general-s form with explicit pi and kappa's at s = (n-k)/k gives an overlapping "
              "enclosure", max(G2[0], L[0]) <= min(G2[1], L[1]), "width %.1e" % float(width(G2)))
    ok1 = True
    for N in (3, 5):
        LN = log_R(n, k, lf=lambda m, N=N: lnfact(m, N=N))
        ok1 = ok1 and max(LN[0], L[0]) <= min(LN[1], L[1])
    rep.check("cross-check: Stirling with 3 and 5 Bernoulli terms gives overlapping enclosures", ok1)
    # moderate n against exact factorials
    n2, k2 = 20000, 7358
    a1 = log_R(n2, k2)
    a2 = log_R(n2, k2, lf=lnfact_exact)
    rep.check("cross-check at n = 20000, k = 7358: Stirling enclosure overlaps the one from the exact integers m!",
              max(a1[0], a2[0]) <= min(a1[1], a2[1]), "R_{20000,7358} = %s..." % dec_floor(exp_iv(a2)[0], 25))
    getcontext().prec = 90
    rep.check("cross-check (non-rigorous, other code): 90-digit decimal, pi by Gauss-Legendre AGM, Stirling with 6 "
              "terms, decimal ln and exp, lies in the enclosure up to 1e-80", abs(decimal_R(n, k) -
              Decimal(R[0]) / Decimal(2) ** P) < Decimal(10) ** -80)
    lg = (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) - math.lgamma(k // 2 + 1) +
          (k / n) * math.lgamma(n // 2 + 1)) / n + (1 - k / n) * math.log(1 - k / n) + (k / n) * math.log(k / n)
    rep.check("sanity: double-precision lgamma gives the same R to 1e-11", abs(math.exp(lg) - float(Fraction(R[0], ONE)))
              < 1e-11, "%.15f" % math.exp(lg))
    # ---- the separate limit exp(1/(2e))
    elo, ehi = e_bounds()
    x = qiv(1 / (2 * ehi))[0], qiv(1 / (2 * elo))[1]
    Linf = exp_iv(x)
    rep.info("exp(1/(2e)) enclosure", enc(Linf, 45))
    rep.check("exp(1/(2e)) lies inside the claimant's Arb ball (radius 3.43e-76)", inside_ball(Linf, ARB_LIMIT))
    rep.check("exp(1/(2e)) > 1.20194336847 (the page's limit decimal; the limit itself is the paper's prose argument)",
              Linf[0] > qiv(Fraction(LIMIT_CLAIM))[1])
    rep.check("exp(1/(2e)) rounds to 1.2019433684703144671942411393813 (FINAL_LIVE_CHECK's Arb digits, rounded to "
              "nearest; truncation gives ...3812) and begins 1.20194336847031446719424113938 (the README's)",
              digits_round(Linf, LIMIT_DIGITS) and not digits_agree(Linf, LIMIT_DIGITS) and
              digits_agree(Linf, "1.20194336847031446719424113938"))
    rep.check("the finite R_{n,k} is below the limit, as expected", R[1] < Linf[0],
              "limit - R = %.4e" % float(Fraction(Linf[0] - R[1], ONE)))
    conv = []
    prev = R
    mono = True
    for ex in (14, 16, 18):
        N = 10 ** ex
        kk = 2 * round(Fraction(N) / (2 * elo))
        RR = exp_iv(log_R(N, kk))
        mono = mono and prev[1] < RR[0] < Linf[0]
        conv.append("n = 10^%d: %s" % (ex, dec_floor(RR[0], 22)))
        prev = RR
    rep.check("R_{n,k} with k the even integer nearest n/e increases toward exp(1/(2e)) at n = 10^14, 10^16, 10^18",
              mono, "; ".join(conv))
    if do_controls:
        controls(rep, n, k, s, claim, L, R, Linf)
    rep.info("what rests on mathematics", "Patankar's criterion (Gram-Schmidt / Szasz, AM-GM, zonotope intrinsic "
             "volumes, Steiner), the positivity of the dropped terms, the Stirling remainder theorem; the limit "
             "exp(1/(2e)) is a separate prose argument (its decimal is enclosed above)")
    return rep.done()


def controls(rep, n, k, s, claim, L, R, Linf):
    print("  -- red controls (each must fire)", flush=True)
    for c in ("1.2019433685", "1.20194336844316"):
        rep.check("control: a bound raised to %s (just above the enclosure) is refuted" % c,
                  R[1] < qiv(Fraction(c))[0])
    rep.check("control: a limit bound raised to 1.20194336848 is refuted", Linf[1] < qiv(Fraction("1.20194336848"))[0])
    out = []
    disjoint = True
    for nn, kk in ((n, k - 2), (n, k + 2), (n - 2, k), (n + 2, k)):
        RR = exp_iv(log_R(nn, kk))
        disjoint = disjoint and (RR[1] < R[0] or RR[0] > R[1])
        out.append("n%+d, k%+d: %s (%s)" % (nn - n, kk - k, dec_floor(RR[0], 30),
                                           "below" if RR[1] < R[0] else "above"))
    rep.check("control: k -+ 2 and n -+ 2 give enclosures disjoint from the certificate's", disjoint, "; ".join(out))
    rep.info("note", "k + 2 gives a larger term than the certificate's k (by %.2e): k is not the exact maximiser over "
             "even k at n = 10^12 (none is claimed)" % float(Fraction(exp_iv(log_R(n, k + 2))[0] - R[1], ONE)))
    Rh = exp_iv(log_R(n, n // 2))
    rep.check("control: the valid but insufficient half-dimensional term k = n/2 does not reach the claim",
              Rh[1] < qiv(claim)[0], "R = %s" % dec_floor(Rh[0], 12))
    try:
        log_R(n, k + 1)
        refused = False
    except ValueError:
        refused = True
    rep.check("control: odd k is refused by the factorial form", refused)
    sp = s * (1 + Fraction(1, 10 ** 6))
    Rs = exp_iv(log_R_general(n, k, sp))
    rep.check("control: a side 1e-6 off (n-k)/k gives a strictly smaller, disjoint enclosure (s = (n-k)/k is the "
              "maximiser)", Rs[1] < R[0], "R(s') = %s" % dec_floor(Rs[0], 30))
    Rbad = exp_iv(log_R_general(n, k, Fraction(13, 1000)))
    rep.check("control: a scale error (normalized side 13/1000) does not reach the claim", Rbad[1] < qiv(claim)[0],
              "R = %s" % dec_floor(Rbad[0], 12))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        print(__doc__)
        sys.exit(2)
    ok = run(args[0], "--controls" in sys.argv)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
