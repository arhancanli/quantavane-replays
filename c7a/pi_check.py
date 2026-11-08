"""Replay of the 7a row [B2026], mu(pi) < 7.101862832357 (Bai, arXiv 2609.11276v2, Theorem 1.1), from the paper's
LaTeX source alone: the numbers sigma = r + C and tau = -s - C and the bound 1 + sigma/tau.

Our own code, standard library only. Exact rationals where the data are rational; elsewhere interval arithmetic of our
own (Python decimal, 60 digits, outward rounding; ln, exp and sqrt are correctly rounded in decimal and every result is
widened by one unit in the last place).

1. The periodic set E = {u in [0,1): {a u + 1/2} + 2{b u} < {c u}} for (a, b, c) = (1857, 3714, 5570) (the paper's
   first display in Section 3), found directly: on each interval between consecutive multiples of 1/3714 and 1/5570 the floors are constant and the
   condition is linear in u. We compare the components with the paper's count (2784), measure and first endpoint.
2. omega = sum over the components [l, r) of psi(r) - psi(l) (Section 3), with psi(x) = psi(x + 100) -
   sum_{q<100} 1/(x + q) and the asymptotic series at x + 100 with B_2..B_22, remainder |B_24|/(24 y^24); a second
   route sums the series sum_q (1/(q + l) - 1/(q + r)) to q < K with the paper's integral tail bounds.
   C = (7430 - 4645 log 2 - omega)/5570.
3. r = (log 25 + 3716 log 2 + log S(zeta) - 5570 log zeta)/5570 (Proposition 4.2), zeta the root in (0,1) of
   m(x) = 3714 x/(1+x) + 3714 x P'(x)/P(x) + 7430 x/(1-x) = 5570, P = 2 + 6x + 9x^2 + 6x^3 + 2x^4 (m is increasing).
   Cross-check: 25 (1+zeta)^2/(1-zeta)^2 is the real root of the stationary cubic p.
4. s = F(y_c) = alpha log|y_c| + beta log|y_c^2 + 6 y_c + 25| - log|25 - y_c| (Proposition 4.3), alpha = 1857/5570,
   beta = 3714/5570, y_c the root of p(y) = 3715y^3 - 232119y^2 - 928475y - 1160625 with positive imaginary part.
5. sigma = r + C, tau = -s - C, and 1 + sigma/tau against the stated 7.1018628323563507957... and the bound
   7.101862832357.
Not checked: the arithmetic of Section 2 (the denominators and the removable primes), the prime-number-theorem
passage, the contour argument of Section 4 (that F along the arc has its unique maximum at y_c), and Hata's lemma.

Usage: python3 -I pi_check.py CERTS_ROOT [--controls]
"""
import hashlib
import os
import re
import sys
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction as Fr
from math import floor

PREC = 60
CF = Context(prec=PREC, rounding=ROUND_FLOOR)
CC = Context(prec=PREC, rounding=ROUND_CEILING)
CN = Context(prec=PREC + 5, rounding=ROUND_HALF_EVEN)


def fl(q):
    q = Fr(q)
    return CF.divide(D(q.numerator), D(q.denominator))


def ce(q):
    q = Fr(q)
    return CC.divide(D(q.numerator), D(q.denominator))


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        if isinstance(lo, (Fr, int)) and hi is None:
            self.lo, self.hi = fl(lo), ce(lo)
        else:
            self.lo, self.hi = D(lo), D(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(a, b):
        b = b if isinstance(b, I) else I(b)
        return I(CF.add(a.lo, b.lo), CC.add(a.hi, b.hi))

    __radd__ = __add__

    def __neg__(a):
        return I(a.hi.copy_negate(), a.lo.copy_negate())

    def __sub__(a, b):
        b = b if isinstance(b, I) else I(b)
        return I(CF.subtract(a.lo, b.hi), CC.subtract(a.hi, b.lo))

    def __rsub__(a, b):
        return I(b) - a

    def __mul__(a, b):
        b = b if isinstance(b, I) else I(b)
        ps = [(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
        return I(min(CF.multiply(x, y) for x, y in ps), max(CC.multiply(x, y) for x, y in ps))

    __rmul__ = __mul__

    def __truediv__(a, b):
        b = b if isinstance(b, I) else I(b)
        assert b.lo > 0 or b.hi < 0
        ps = [(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
        return I(min(CF.divide(x, y) for x, y in ps), max(CC.divide(x, y) for x, y in ps))

    def __rtruediv__(a, b):
        return I(b) / a

    def __repr__(self):
        return "[%s, %s]" % (str(self.lo)[:34], str(self.hi)[:34])


def _dn(r):
    return CF.plus(CN.next_minus(r))


def _up(r):
    return CC.plus(CN.next_plus(r))


def ln(a):
    assert a.lo > 0
    return I(_dn(CN.ln(a.lo)), _up(CN.ln(a.hi)))


def sqrt(a):
    assert a.lo >= 0
    return I(_dn(CN.sqrt(a.lo)), _up(CN.sqrt(a.hi)))


# ----------------------------------------------------------------------------------------------- 1. the set E
A, B, C = 1857, 3714, 5570


def periodic_set(a, b, c):
    """Components of {u in [0,1): {a u + 1/2} + 2{b u} < {c u}} as exact rationals."""
    pts = sorted({Fr(k, 2 * a) for k in range(2 * a + 1)} | {Fr(k, b) for k in range(b + 1)} |
                 {Fr(k, c) for k in range(c + 1)})
    d = a + 2 * b - c
    pieces = []
    for u0, u1 in zip(pts, pts[1:]):
        m = (u0 + u1) / 2
        f1 = floor(a * m + Fr(1, 2))
        f2 = floor(b * m)
        f3 = floor(c * m)
        # on (u0, u1) the floors are constant (no breakpoint inside); the condition is
        # (a u + 1/2 - f1) + 2 (b u - f2) < c u - f3, i.e. d u < f1 + 2 f2 - f3 - 1/2
        rhs = Fr(f1 + 2 * f2 - f3) - Fr(1, 2)
        if d > 0:
            top = rhs / d
            if top > u0:
                pieces.append((u0, min(u1, top)))
        else:
            raise ValueError("d <= 0 not handled")
    comps = []
    for l, r in pieces:
        if comps and comps[-1][1] == l:
            comps[-1] = (comps[-1][0], r)
        else:
            comps.append((l, r))
    return comps


# ----------------------------------------------------------------------------------------------- 2. digamma
BERN = [Fr(1, 6), Fr(-1, 30), Fr(1, 42), Fr(-1, 30), Fr(5, 66), Fr(-691, 2730), Fr(7, 6), Fr(-3617, 510),
        Fr(43867, 798), Fr(-174611, 330), Fr(854513, 138), Fr(-236364091, 2730)]   # B_2 .. B_24


def psi_minus_ln100(x):
    """psi(x) - log 100 for rational 0 < x < 1, as an interval: psi(x) = psi(y) - sum_{q<100} 1/(x+q), y = x + 100,
    psi(y) = log y - 1/(2y) - sum_{m=1}^{11} B_2m/(2m y^2m) + R, |R| < |B_24|/(24 y^24); log y - log 100 =
    log(1 + x/100) directly in interval arithmetic."""
    xi = I(Fr(x))
    y = xi + 100
    series = -(I(1) / (I(2) * y))
    yp = y * y
    for m in range(1, 12):
        series = series - I(BERN[m - 1]) / (I(2 * m) * yp)
        yp = yp * y * y
    harm = I(0)
    for q in range(100):
        harm = harm + I(1) / (xi + q)
    rem = ce(abs(BERN[11]) / (24 * (Fr(x) + 100) ** 24))
    val = ln(1 + xi / I(100)) + series - harm
    return I(CF.subtract(val.lo, rem), CC.add(val.hi, rem))


def omega_digamma(comps):
    tot = I(0)
    for l, r in comps:
        tot = tot + (psi_minus_ln100(r) - psi_minus_ln100(l))
    return tot


def omega_series(comps, K):
    """sum_k sum_q (1/(q+l) - 1/(q+r)) with the integral tail bounds: log((K+r)/(K+l)) <= tail <= that + f(K)."""
    Ctx = Context(prec=40)
    s = D(0)
    lo_tail = D(0)
    hi_extra = D(0)
    for l, r in comps:
        ld, rd = Ctx.divide(D(l.numerator), D(l.denominator)), Ctx.divide(D(r.numerator), D(r.denominator))
        for q in range(K):
            s = Ctx.add(s, Ctx.subtract(Ctx.divide(1, q + ld), Ctx.divide(1, q + rd)))
        lo_tail = Ctx.add(lo_tail, Ctx.ln(Ctx.divide(K + rd, K + ld)))
        hi_extra = Ctx.add(hi_extra, Ctx.divide(rd - ld, (K + ld) * (K + rd)))
    return s + lo_tail, s + lo_tail + hi_extra


# ----------------------------------------------------------------------------------------------- 3. r
def P_iv(x):
    return 2 + x * (6 + x * (9 + x * (6 + x * 2)))


def dP_iv(x):
    return 6 + x * (18 + x * (18 + x * 8))


def m_of(x):
    """m(x) at an exact rational point, as an interval"""
    xi = I(x)
    return I(3714) * xi / (1 + xi) + I(3714) * xi * dP_iv(xi) / P_iv(xi) + I(7430) * xi / (1 - xi)


def zeta_root(target=5570, iters=200):
    lo, hi = Fr(1, 10), Fr(9, 10)
    assert m_of(lo).hi < target < m_of(hi).lo
    for _ in range(iters):
        mid = (lo + hi) / 2
        v = m_of(mid)
        if v.hi < target:
            lo = mid
        elif v.lo > target:
            hi = mid
        else:
            break
    return lo, hi


def r_value(zlo, zhi):
    """r = (log 25 + 3716 log 2 + log S(zeta) - 5570 log zeta)/5570 on the bracket [zlo, zhi]; the function
    g(x) = log S(x) - 5570 log x has g'(x) = (m(x) - 5570)/x, so on the bracket its range is within g(zlo..zhi)
    evaluated with interval x."""
    x = I(fl(zlo), ce(zhi))
    g = I(3714) * ln(1 + x) + I(3714) * ln(P_iv(x)) - I(7430) * ln(1 - x) - I(5570) * ln(x)
    return (ln(I(25)) + I(3716) * ln(I(2)) + g) / I(5570)


# ----------------------------------------------------------------------------------------------- 4. s
def cubic(y):
    return 3715 * y ** 3 - 232119 * y ** 2 - 928475 * y - 1160625


def real_root():
    lo, hi = Fr(66), Fr(67)
    assert cubic(lo) < 0 < cubic(hi)
    for _ in range(220):
        mid = (lo + hi) / 2
        if cubic(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


class Cx:
    """a complex box re + i im with interval parts"""

    def __init__(self, re, im):
        self.re, self.im = re, im

    def __add__(a, b):
        b = b if isinstance(b, Cx) else Cx(b if isinstance(b, I) else I(b), I(0))
        return Cx(a.re + b.re, a.im + b.im)

    def __mul__(a, b):
        b = b if isinstance(b, Cx) else Cx(b if isinstance(b, I) else I(b), I(0))
        return Cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)

    def abs2(self):
        return self.re * self.re + self.im * self.im


def s_value(rho):
    """y_c from the quadratic factor of p after the real root rho: y^2 + Bq y + Cq with Bq = rho - 232119/3715,
    Cq = 1160625/(3715 rho); y_c = (-Bq + i sqrt(4 Cq - Bq^2))/2."""
    Bq = rho - I(Fr(232119, 3715))
    Cq = I(Fr(1160625, 3715)) / rho
    disc = I(4) * Cq - Bq * Bq
    assert disc.lo > 0
    yc = Cx((-Bq) / I(2), sqrt(disc) / I(2))
    t1 = ln(yc.abs2()) * I(Fr(1, 2))
    q = yc * yc + yc * I(6) + I(25)
    t2 = ln(q.abs2()) * I(Fr(1, 2))
    w = Cx(I(25) - yc.re, -yc.im)
    t3 = ln(w.abs2()) * I(Fr(1, 2))
    s = I(Fr(1857, 5570)) * t1 + I(Fr(3714, 5570)) * t2 - t3
    return s, yc


TEX = {
    "sections/03_prime_saving.tex": [
        r"\left\{1857u+\frac12\right\}+2\{3714u\}<\{5570u\}",
        r"|E|=\frac{1724689}{13797510}.",
        r"=\frac{7430-4645\log2-\omega}{5570}.",
    ],
    "sections/04_asymptotics.tex": [
        r"=2+6z+9z^2+6z^3+2z^4.",
        r"S(z)=\frac{(1+z)^{3714}P(z)^{3714}}{(1-z)^{7430}}.",
        r"r=\frac{\log25+3716\log2+\log S(\zeta)-5570\log\zeta} {5570}.",
        r"p(y)=3715y^3-232119y^2-928475y-1160625.",
        r"F(y)=\alpha\log|y|+\beta\log|y^2+6y+25|-\log|25-y|.",
        r"s=F(y_c),",
    ],
    "sections/05_measure.tex": [
        r"=1+\frac{r+\mathcal C}{-s-\mathcal C}.",
        r"=7.1018628323563507957\ldots <7.101862832357.",
    ],
    "sections/01_introduction.tex": [
        r"(a,b,c)=(1857,3714,5570).",
    ],
}


# ============================================================================ [ZZ2020]: mu(pi) <= 7.103205334137...
ZZ_TEX = [
    r"\cP_n=\biggl\{p>\max\{5,\sqrt{3n}\}:\frac12\le\Bigl\{\frac np\Bigr\}<\frac23\biggr\}",
    r"=\frac{\Gamma'(2/3)}{\Gamma(2/3)}-\frac{\Gamma'(1/2)}{\Gamma(1/2)} =\frac\pi{2\sqrt3}-\log\frac{3\sqrt3}4 =0.64527561\dotsc",
    r"108N^3-2359989N^2+138304N-2048.",
    r"\frac{g'(y)}{g(y)}=\frac{2y^3-125y^2-500y-625}{y(y^2+6y+25)(y-25)}.",
    r"\tilde R(y)=\frac{5g(y)^n}{y-25}, \quad\text{where}\; g(y)=\frac{y(y^2+6y+25)^2}{(y-25)^3},",
    r"&=\log|N_1|-\frac52\,\log2+4-\frac\pi{2\sqrt3}+\log\frac{3\sqrt3}4 =-1.90291648559998\dots",
    r"&=\log N_3-\frac52\,\log2+4-\frac\pi{2\sqrt3}+\log\frac{3\sqrt3}4 =11.613890045331\dots",
    r"=7.10320533413700172750577342281\dotsc.",
]


def arctan_inv(k, terms=60):
    """arctan(1/k) bracketed by consecutive partial sums of the alternating series (terms decrease)"""
    s = Fr(0)
    sums = []
    for j in range(terms):
        s += Fr((-1) ** j, (2 * j + 1) * k ** (2 * j + 1))
        sums.append(s)
    return min(sums[-1], sums[-2]), max(sums[-1], sums[-2])


def pi_iv():
    a5, b5 = arctan_inv(5)
    a239, b239 = arctan_inv(239)
    return I(fl(16 * a5 - 4 * b239), ce(16 * b5 - 4 * a239))


def cubic_real_root(f, lo, hi, iters=220):
    assert f(lo) < 0 < f(hi) or f(lo) > 0 > f(hi)
    up = f(lo) < 0
    for _ in range(iters):
        mid = (lo + hi) / 2
        if (f(mid) < 0) == up:
            lo = mid
        else:
            hi = mid
    return lo, hi


def main_zz(root, controls):
    def out(x):
        print(x, flush=True)

    out("7a row [ZZ2020]: mu(pi) <= 7.103205334137... (Zeilberger and Zudilin, arXiv 1912.06345v2, Part II)")
    status = True
    path = os.path.join(root, "zz2020/main.tex")
    with open(path, "rb") as fh:
        raw = fh.read()
    out("zz2020/main.tex sha256 %s" % hashlib.sha256(raw).hexdigest())
    text = " ".join(raw.decode("utf-8").split())
    for nd in ZZ_TEX:
        hit = " ".join(nd.split()) in text
        status &= hit
        out("  our transcription appears verbatim: %s   %s" % (hit, nd[:100]))
    # the prime saving: E = [1/2, 2/3), so omega = psi(2/3) - psi(1/2), two routes
    om_a = psi_minus_ln100(Fr(2, 3)) - psi_minus_ln100(Fr(1, 2))
    pi = pi_iv()
    s3 = sqrt(I(3))
    om_b = pi / (I(2) * s3) - ln(I(3) * s3 / I(4))
    agree = om_a.lo <= om_b.hi and om_b.lo <= om_a.hi
    out("pi in %s (Machin, alternating-series brackets)" % pi)
    out("omega = psi(2/3) - psi(1/2) in %s (digamma route); pi/(2 sqrt 3) - log(3 sqrt 3/4) in %s; overlap: %s; "
        "the stated 0.64527561... is a truncation: %s" % (om_a, om_b, agree,
                                                        D("0.64527561") <= om_a.lo and om_a.hi < D("0.64527562")))
    status &= agree and D("0.64527561") <= om_a.lo
    # the indicial cubic: one real root N_3, and |N_1|^2 = 2048/(108 N_3) by Vieta
    fN = lambda x: 108 * x ** 3 - 2359989 * x ** 2 + 138304 * x - 2048
    disc = 18 * 108 * (-2359989) * 138304 * (-2048) - 4 * (-2359989) ** 3 * (-2048) + (-2359989) ** 2 * 138304 ** 2 \
        - 4 * 108 * 138304 ** 3 - 27 * 108 ** 2 * 2048 ** 2
    lo, hi = cubic_real_root(fN, Fr(21851), Fr(21852))
    N3 = I(fl(lo), ce(hi))
    absN1 = sqrt(I(Fr(2048, 108)) / N3)
    out("discriminant of 108N^3 - 2359989N^2 + 138304N - 2048 is %s (< 0: one real root): %s" % (str(disc)[:12] + "...",
                                                                                              disc < 0))
    out("N_3 in %s; |N_1| = sqrt(2048/(108 N_3)) in %s (stated 21851.691396... and 0.029458495928...)" % (N3, absN1))
    status &= disc < 0 and D("21851.691396") <= N3.lo and D("0.029458495928") <= absN1.lo and absN1.hi < D("0.029458495929")
    # second route: the saddle points y of g'(y)/g(y), N_j = g(y_j)
    fy = lambda y: 2 * y ** 3 - 125 * y ** 2 - 500 * y - 625
    ylo, yhi = cubic_real_root(fy, Fr(66), Fr(67))
    y3 = I(fl(ylo), ce(yhi))
    g3 = y3 * (y3 * y3 + I(6) * y3 + I(25)) * (y3 * y3 + I(6) * y3 + I(25)) / ((y3 - I(25)) * (y3 - I(25)) * (y3 - I(25)))
    Bq = y3 - I(Fr(125, 2))
    Cq = I(Fr(625, 2)) / y3
    dq = I(4) * Cq - Bq * Bq
    y1 = Cx((-Bq) / I(2), sqrt(dq) / I(2))
    y1sq = y1 * y1
    q = y1sq + y1 * I(6) + I(25)
    num2 = y1.abs2() * q.abs2() * q.abs2()
    w = Cx(y1.re - I(25), y1.im)
    den2 = w.abs2() * w.abs2() * w.abs2()
    g1 = sqrt(num2 / den2)
    r2 = g3.lo <= N3.hi and N3.lo <= g3.hi and g1.lo <= absN1.hi and absN1.lo <= g1.hi
    out("second route: y_3 in %s, g(y_3) in %s; |g(y_1)| in %s; agree with N_3 and |N_1|: %s" % (y3, g3, g1, r2))
    status &= r2
    cst = I(4) - I(Fr(5, 2)) * ln(I(2)) - om_a
    rate_I = ln(absN1) + cst
    rate_b = ln(N3) + cst
    Bzz = 1 + rate_b / (-rate_I)
    out("log|I'_n|/n -> %s (stated -1.90291648559998...); log b'_n/n -> %s (stated 11.613890045331...)" % (rate_I, rate_b))
    out("1 + rate_b/(-rate_I) in %s (stated 7.10320533413700172750577342281...)" % Bzz)
    ok = (D("-1.90291648559999") < rate_I.lo and rate_I.hi < D("-1.90291648559998")
          and D("11.613890045331") <= rate_b.lo and rate_b.hi < D("11.613890045332")
          and D("7.10320533413700172750577342281") <= Bzz.lo and Bzz.hi < D("7.10320533413700172750577342282"))
    out("all three agree with the paper's printed digits: %s" % ok)
    status &= ok
    below = Bzz.lo > D("7.103205334137")
    out("the page prints 7.103205334137, which is %.4e below the bound (7.103205334137... in the paper's title): %s"
        % (float(CF.subtract(Bzz.lo, D("7.103205334137"))), below))
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        om_bad = psi_minus_ln100(Fr(3, 4)) - psi_minus_ln100(Fr(1, 2))
        r1 = om_bad.lo > om_a.hi + D("0.01")
        out("[ctl] the set [1/2, 3/4) instead of [1/2, 2/3) changes omega to %s: %s" % (om_bad, r1))
        ctl &= r1
        r2c = not (Bzz.hi < D("7.103205334136"))
        out("[ctl] a bound 7.103205334136 is refused: %s" % r2c)
        ctl &= r2c
        lo2, hi2 = cubic_real_root(lambda x: 108 * x ** 3 - 2359989 * x ** 2 + 138304 * x - 2049, Fr(21851), Fr(21852))
        absN1b = sqrt(I(Fr(2049, 108)) / I(fl(lo2), ce(hi2)))
        r3 = absN1b.lo > absN1.hi
        out("[ctl] the constant term 2048 -> 2049 moves |N_1| to %s: %s" % (absN1b, r3))
        ctl &= r3
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT: %s" % ("PASS (1 + 11.613890045331.../1.90291648559998... = 7.10320533413700172750577342281...; two "
                        "routes for the prime saving and for N_1, N_3)" if status else "FAIL"))
    return 0 if status else 1


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    controls = "--controls" in sys.argv
    if "--row=zz2020" in sys.argv:
        return main_zz(root, controls)

    def out(x):
        print(x, flush=True)

    out("7a row [B2026]: mu(pi) < 7.101862832357 (Bai, arXiv 2609.11276v2, Theorem 1.1)")
    status = True
    for rel, needles in TEX.items():
        path = os.path.join(root, "b2026", rel)
        with open(path, "rb") as fh:
            raw = fh.read()
        out("b2026/%s sha256 %s" % (rel, hashlib.sha256(raw).hexdigest()))
        text = " ".join(raw.decode("utf-8").split())
        for nd in needles:
            hit = " ".join(nd.split()) in text
            status &= hit
            out("  our transcription appears verbatim: %s   %s" % (hit, nd))
    comps = periodic_set(A, B, C)
    meas = sum(r - l for l, r in comps)
    estats = (len(comps), meas, comps[0][0])
    eok = estats == (2784, Fr(1724689, 13797510), Fr(1, 3714))
    out("E: %d components, measure %s, first left endpoint %s (paper: 2784, 1724689/13797510, 1/3714): %s"
        % (estats[0], estats[1], estats[2], eok))
    status &= eok
    om = omega_digamma(comps)
    out("omega (digamma route) in %s" % om)
    s_lo, s_hi = omega_series(comps, 2000)
    agree = s_lo - D("1e-30") <= om.hi and om.lo <= s_hi + D("1e-30")
    out("omega (series route, K = 2000, tail bounds) in [%s, %s]; the two routes overlap: %s"
        % (str(s_lo)[:22], str(s_hi)[:22], agree))
    status &= agree
    Cv = (I(7430) - I(4645) * ln(I(2)) - om) / I(5570)
    out("C in %s (paper: 0.539993819198880114452898006 < C < ...061)" % Cv)
    zlo, zhi = zeta_root()
    rv = r_value(zlo, zhi)
    out("zeta in [%s, %s]; r in %s (paper: 3.33194398617829933828283250382326 < r < ...327)"
        % (str(fl(zlo))[:30], str(ce(zhi))[:30], rv))
    rlo, rhi = real_root()
    rho = I(fl(rlo), ce(rhi))
    zi = I(fl(zlo), ce(zhi))
    rho_from_zeta = I(25) * (1 + zi) * (1 + zi) / ((1 - zi) * (1 - zi))
    cross = rho_from_zeta.lo <= rho.hi and rho.lo <= rho_from_zeta.hi
    out("rho (real root of p) in %s; 25 (1+zeta)^2/(1-zeta)^2 in %s; overlap: %s" % (rho, rho_from_zeta, cross))
    status &= cross
    sv, yc = s_value(rho)
    out("y_c = %s + i %s; s in %s (paper: -1.174543941310025108442583231 < s < ...192)" % (yc.re, yc.im, sv))
    sigma = rv + Cv
    tau = -sv - Cv
    Bv = 1 + sigma / tau
    out("sigma in %s; tau in %s" % (sigma, tau))
    out("1 + sigma/tau in %s" % Bv)
    stated = D("7.1018628323563507957") <= Bv.lo and Bv.hi < D("7.1018628323563507958")
    bound = Bv.hi < D("7.101862832357")
    paper_encl = D("7.101862832356350795733674505") < Bv.lo and Bv.hi < D("7.101862832356350795733675497")
    out("the stated 7.1018628323563507957... is its truncation: %s; inside the paper's enclosure (7.101862832356350795733674505,"
        " ...675497): %s; below the printed 7.101862832357: %s" % (stated, paper_encl, bound))
    status &= stated and bound and paper_encl
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        comps2 = periodic_set(1856, B, C)
        r1 = len(comps2) != 2784 or sum(r - l for l, r in comps2) != Fr(1724689, 13797510)
        out("[ctl] E with a = 1856: %d components, measure %s: differs: %s" % (len(comps2), sum(r - l for l, r in comps2),
                                                                             r1))
        ctl &= r1
        z2 = zeta_root(5571)
        r2 = abs(r_value(*z2).lo - rv.lo) > D("1e-10")
        out("[ctl] zeta from m = 5571 moves r to %s: %s" % (r_value(*z2), r2))
        ctl &= r2
        r3 = not (Bv.hi < D("7.101862832356"))
        out("[ctl] the bound 7.101862832356 is refused: %s" % r3)
        ctl &= r3
        sv_bad = I(Fr(1857, 5570)) * ln(rho) + I(Fr(3714, 5570)) * ln(rho * rho + I(6) * rho + I(25)) - \
            ln(I(0) - (I(25) - rho))
        r4 = abs(sv_bad.lo - sv.lo) > D("1e-6")
        same_r = sv_bad.lo <= rv.hi + D("1e-50") and rv.lo <= sv_bad.hi + D("1e-50")
        out("[ctl] F at the real critical point rho instead of y_c gives %s, not s: %s (it equals r, a second route to r: "
            "%s)" % (sv_bad, r4, same_r))
        r4 = r4 and same_r
        ctl &= r4
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT: %s" % ("PASS (1 + sigma/tau = 7.1018628323563507957... < 7.101862832357; E, omega, C, r and s "
                        "re-derived; the analysis behind them is not checked)" if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
