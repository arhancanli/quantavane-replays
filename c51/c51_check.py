#!/usr/bin/env python3
"""Exact replay of the C_51 lower bound 0.5850788196744* [R2026b] (the Erdos maximum-term constant, constants/51a.md;
the row of PR #184), from He and Tang's reduction (arXiv:2602.12217, Theorem 2.8 and Lemma 3.1) and the two exact
parameters alone. Roehrig's certify_maximum_term.py and audit.py were neither read nor run.

  Theorem 2.8 (He-Tang). For K > 1 and |eps| = 1, f(z) = sum_{n>=0} eps^(n(n-1)/2) K^(-n(n+1)/2) z^n is transcendental
  entire and beta(f) = 1/A(K, eps), A = max_{|z|=1} |k(z)|, k the bilateral series.
  Lemma 3.1. A = max_theta |P(theta)|, P(theta) = 2 sum_{n>=0} c_n cos((2n+1) theta), c_n = K^(-T_n) e^(i alpha T_n),
  T_n = n(n+1)/2, eps = e^(i alpha).

So C_51 >= 1/A, and a row C_51 >= b needs A < 1/b. Our method (neither Roehrig's Bernstein subdivision in x = cos^2,
nor He and Tang's mesh with a Lipschitz bound):
  - P = P_N + R_N with N = 12, |R_N| <= 2 sum_{n>N} K^(-T_n) <= 2 K^(-T_(N+1)) / (1 - K^(-(N+2))), an exact rational
    (T_(n+1) - T_n = n + 1 >= N + 2 for n > N).
  - |P_N|^2 = F = A^2 + B^2 with A = sum a_n cos((2n+1) t), B = sum b_n cos((2n+1) t), (a_n, b_n) = 2 (Re c_n, Im c_n).
    P(-t) = P(t) and P(pi - t) = -P(t), so t in [0, pi/2] suffices; we cover [0, 1.5708], which contains it.
  - On a cell [m - h, m + h]: F(m + u) <= F(m) + F'(m) u + F''(m) u^2/2 + G3 |u|^3/6, with
    G3 = 2 (g0 g3 + 3 g1 g2) >= sup |F'''|, g_k = sum_{n<=N} 2 K^(-T_n) (2n+1)^k (F''' = 2 (A A''' + 3 A' A'' + B B'''
    + 3 B' B''), and Cauchy-Schwarz on the vectors (A^(j), B^(j)), whose length is at most g_j). F, F', F'' at m in
    decimal interval arithmetic (60 digits, every operation rounded outward); the maximum of the quadratic over
    |u| <= h, and every comparison, in exact rationals. A cell above the target is bisected.
  - A lower bound A >= |P_N(t)| - |R_N| at the maxima found (located in floating point, then evaluated in intervals).
  - e^(i alpha T_n) as powers of e^(i alpha) (no reduction modulo 2 pi); pi (for He and Tang's mesh point only) from
    Machin's formula with the series tail bounded.

Rows (with --controls, red controls that must fail):
  r2026b   K = 3.5681854310698435, alpha = 3.9615431126097391 (certificates/maximum_term.json of Roehrig's package at
           the pinned commit, SHA-256 checked), claim 0.5850788196744
  s2026    the previous parameters K = 3.568182317714, alpha = 3.961543335688 (as Roehrig's note states them): beta
           enclosed, against the 0.585078819653 of [S2026] and Roehrig's beta < 0.585078819653352950854632621542
  ht2026   He and Tang's own K0 = 7137/2000, alpha0 = 198074929/50000000: beta > 0.58507 (their Theorem 1.2), with their
           mesh maximum at j = 1794536 recomputed (a test of our arithmetic against theirs)

  python3 -I c51_check.py [r2026b|s2026|ht2026|selftest|all] [CERTS] [--controls]
Python 3.9+, standard library only.
"""
import decimal
import hashlib
import json
import math
import os
import random
import sys
import time
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PREC = 60
LO = decimal.Context(prec=PREC, rounding=decimal.ROUND_FLOOR)
HI = decimal.Context(prec=PREC, rounding=decimal.ROUND_CEILING)
NEAR = decimal.Context(prec=PREC, rounding=decimal.ROUND_HALF_EVEN)
N_TRUNC = 12
COVER_END = Fraction(15708, 10000)        # > pi/2 = 1.5707963267948966...
CERT = ("r2026b/maximum_term.json", "dcf9d1bc0f6d7b295e0eff548cf73faf092efd8d038829e4d19755f07df33683")


# ---------------------------------------------------------------------------------------------- interval arithmetic
class I:
    """A closed interval [lo, hi] of decimals; every operation rounds outward."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo if isinstance(lo, Decimal) else Decimal(lo)
        self.hi = self.lo if hi is None else (hi if isinstance(hi, Decimal) else Decimal(hi))
        assert self.lo <= self.hi, (self.lo, self.hi)

    @staticmethod
    def frac(x):
        """The rational x: its quotient rounded down and up."""
        x = Fraction(x)
        n, d = Decimal(x.numerator), Decimal(x.denominator)
        return I(LO.divide(n, d), HI.divide(n, d))

    def __add__(self, o):
        return I(LO.add(self.lo, o.lo), HI.add(self.hi, o.hi))

    def __sub__(self, o):
        return I(LO.subtract(self.lo, o.hi), HI.subtract(self.hi, o.lo))

    def __neg__(self):
        return I(self.hi.copy_negate(), self.lo.copy_negate())

    def __mul__(self, o):
        a, b, c, d = self.lo, self.hi, o.lo, o.hi
        if a >= 0 and c >= 0:
            return I(LO.multiply(a, c), HI.multiply(b, d))
        # every product rounded in both directions (never in decimal's default 28-digit context)
        lo = min(LO.multiply(a, c), LO.multiply(a, d), LO.multiply(b, c), LO.multiply(b, d))
        hi = max(HI.multiply(a, c), HI.multiply(a, d), HI.multiply(b, c), HI.multiply(b, d))
        return I(lo, hi)

    def scale(self, k):
        """Multiplication by an exact rational k."""
        return self * I.frac(k)

    def muli(self, k):
        """Multiplication by an integer k >= 0 (k as a decimal is exact)."""
        k = Decimal(k)
        return I(LO.multiply(self.lo, k), HI.multiply(self.hi, k))

    def sq(self):
        a, b = self.lo, self.hi
        if a >= 0:
            return I(LO.multiply(a, a), HI.multiply(b, b))
        if b <= 0:
            return I(LO.multiply(b, b), HI.multiply(a, a))
        return I(Decimal(0), max(HI.multiply(a, a), HI.multiply(b, b)))

    def mag(self):
        return max(self.lo.copy_abs(), self.hi.copy_abs())

    def widen(self, r):
        return I(LO.subtract(self.lo, r), HI.add(self.hi, r))

    def flo(self):
        return Fraction(self.lo)

    def fhi(self):
        return Fraction(self.hi)


ZERO, ONE = I(0), I(1)


class C:
    """A complex rectangle re + i im."""
    __slots__ = ("re", "im")

    def __init__(self, re, im):
        self.re, self.im = re, im

    def __mul__(self, o):
        return C(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)


def sincos(x):
    """(sin x, cos x) for an exact rational x with |x| <= 8: Taylor series, the remainder bounded by
    |x|^(k+1)/(k+1)! / (1 - |x|/(k+2)) once k + 2 > 2|x|."""
    x = Fraction(x)
    assert abs(x) <= 8
    X = I.frac(x)
    s, c, term = ZERO, ONE, ONE
    tol = Fraction(1, 10 ** (PREC + 2))
    k = 0
    while True:
        k += 1
        term = term * X * I.frac(Fraction(1, k))
        if k % 2:
            s = s + term if k % 4 == 1 else s - term
        else:
            c = c + term if k % 4 == 0 else c - term
        if k + 2 > 2 * abs(x) + 2 and Fraction(term.mag()) * abs(x) < tol:
            break
    nxt = Fraction(term.mag()) * abs(x) / (k + 1)                 # bounds |x|^(k+1)/(k+1)!
    rem = nxt / (1 - abs(x) / (k + 2))
    r = HI.divide(Decimal(rem.numerator), Decimal(rem.denominator)) if rem else Decimal(0)
    return s.widen(r), c.widen(r)


def isqrt_lo(x):
    """A decimal y <= sqrt(x) for a rational x >= 0 (checked exactly)."""
    x = Fraction(x)
    if x <= 0:
        return Decimal(0)
    y = LO.sqrt(LO.divide(Decimal(x.numerator), Decimal(x.denominator)))
    while Fraction(y) ** 2 > x:
        y = LO.next_minus(y)
    return y


def isqrt_hi(x):
    """A decimal y >= sqrt(x) for a rational x >= 0 (checked exactly)."""
    x = Fraction(x)
    y = HI.sqrt(HI.divide(Decimal(x.numerator), Decimal(x.denominator)))
    while Fraction(y) ** 2 < x:
        y = HI.next_plus(y)
    return y


def pi_interval():
    """pi = 16 arctan(1/5) - 4 arctan(1/239) (Machin), each arctan by its alternating series with the tail bounded by
    the first omitted term."""
    def atan_inv(m):
        tot, k = ZERO, 0
        while True:
            t = I.frac(Fraction(1, (2 * k + 1) * m ** (2 * k + 1)))
            if Fraction(t.hi) < Fraction(1, 10 ** (PREC + 3)):
                return tot.widen(t.hi)
            tot = tot + t if k % 2 == 0 else tot - t
            k += 1
    return atan_inv(5).scale(16) - atan_inv(239).scale(4)


# -------------------------------------------------------------------------------------------------- the function P
class Theta:
    """P_N(t) = 2 sum_{n<=N} c_n cos((2n+1) t) for one parameter pair, with everything a cell needs."""

    def __init__(self, K, alpha, N=N_TRUNC):
        self.K, self.alpha, self.N = Fraction(K), Fraction(alpha), N
        assert self.K > 1
        T = [n * (n + 1) // 2 for n in range(N + 2)]
        s, c = sincos(self.alpha)
        eps = C(c, s)
        epow = [C(ONE, ZERO)]                                    # eps^j for j <= T_N
        for _ in range(T[N]):
            epow.append(epow[-1] * eps)
        self.mods = [self.K ** -T[n] for n in range(N + 1)]      # |c_n| = K^(-T_n), exact
        self.a, self.b = [], []
        for n in range(N + 1):
            m = I.frac(2 * self.mods[n])
            self.a.append(epow[T[n]].re * m)
            self.b.append(epow[T[n]].im * m)
        # |R_N| <= 2 K^(-T_(N+1)) / (1 - K^(-(N+2)))
        self.tail = 2 * self.K ** -T[N + 1] / (1 - self.K ** -(N + 2))
        g = [sum(2 * self.mods[n] * (2 * n + 1) ** k for n in range(N + 1)) for k in range(4)]
        self.G3 = 2 * (g[0] * g[3] + 3 * g[1] * g[2])          # >= sup |F'''|
        self.G2 = 2 * (g[1] * g[1] + g[0] * g[2])              # >= sup |F''|
        self.g = g

    def derivs(self, t):
        """(F, F', F'') at the exact rational t, F = |P_N|^2, as intervals."""
        s1, c1 = sincos(t)
        z1 = C(c1, s1)
        z2 = z1 * z1
        z = z1
        A0 = A1 = A2 = B0 = B1 = B2 = ZERO
        for n in range(self.N + 1):
            k = 2 * n + 1
            if n:
                z = z * z2
            cr, si = z.re, z.im                                   # cos(k t), sin(k t)
            ks = si.muli(k)
            kkc = cr.muli(k * k)
            A0, B0 = A0 + self.a[n] * cr, B0 + self.b[n] * cr
            A1, B1 = A1 - self.a[n] * ks, B1 - self.b[n] * ks
            A2, B2 = A2 - self.a[n] * kkc, B2 - self.b[n] * kkc
        F0 = A0.sq() + B0.sq()
        F1 = (A0 * A1 + B0 * B1).muli(2)
        F2 = (A1.sq() + A0 * A2 + B1.sq() + B0 * B2).muli(2)
        return F0, F1, F2

    def float_abs(self, t):
        """|P_N(t)| in floating point (to locate maxima only)."""
        re = im = 0.0
        for n in range(self.N + 1):
            cs = math.cos((2 * n + 1) * t)
            re += float(self.a[n].lo) * cs
            im += float(self.b[n].lo) * cs
        return math.hypot(re, im)


def quad_max(p, r, h):
    """max over 0 <= u <= h of p u + r u^2 (exact rationals)."""
    best = max(Fraction(0), p * h + r * h * h)
    if r < 0 and p > 0:
        u = -p / (2 * r)
        if u < h:
            best = max(best, -p * p / (4 * r))
    return best


def cell_bound(P, a, b):
    """An upper bound for F = |P_N|^2 on [a, b] (exact rational)."""
    m, h = (a + b) / 2, (b - a) / 2
    F0, F1, F2 = P.derivs(m)
    p = max(F1.fhi(), -F1.flo())
    r = F2.fhi() / 2
    return F0.fhi() + quad_max(p, r, h) + P.G3 * h ** 3 / 6, F0


def range_ok(end):
    """The cover [0, end] contains [0, pi/2] (pi from Machin's formula, compared at its upper end)."""
    return Fraction(end) >= pi_interval().fhi() / 2


def cover(P, target, start=Fraction(0), end=COVER_END, init=1024, max_depth=64):
    """Check F <= target on [start, end] by adaptive bisection. Returns a dict of statistics; ok is False if a cell
    could not be accepted within max_depth."""
    stack = [(start + (end - start) * j / init, start + (end - start) * (j + 1) / init, 0) for j in range(init)]
    stack.reverse()
    nodes = leaves = 0
    deepest, smallest, worst = 0, end - start, Fraction(0)
    covered = []                                             # accepted cells, in order (checked to tile the range)
    while stack:
        a, b, d = stack.pop()
        nodes += 1
        ub, _ = cell_bound(P, a, b)
        if ub <= target:
            leaves += 1
            covered.append((a, b))
            deepest, smallest, worst = max(deepest, d), min(smallest, b - a), max(worst, ub)
            continue
        if d >= max_depth:
            return dict(ok=False, nodes=nodes, fail=(a, b, ub))
        mid = (a + b) / 2
        stack.append((mid, b, d + 1))
        stack.append((a, mid, d + 1))
    tiles = covered[0][0] == start and covered[-1][1] == end and all(x[1] == y[0] for x, y in zip(covered, covered[1:]))
    return dict(ok=tiles, nodes=nodes, leaves=leaves, deepest=deepest, smallest=smallest, worst=worst, tiles=tiles)


def find_maxima(P, grid=40000):
    """Local maxima of |P_N| on [0, pi/2], located in floating point and polished by Newton steps on F' = 0 in
    interval midpoints; returned as exact rationals with a rigorous lower bound for |P_N| there."""
    hs = math.pi / 2 / grid
    vals = [P.float_abs(j * hs) for j in range(grid + 1)]
    cands = [0] if vals[0] >= vals[1] else []
    cands += [j for j in range(1, grid) if vals[j] >= vals[j - 1] and vals[j] >= vals[j + 1]]
    out = []
    for j in cands:
        t = Fraction(0) if j == 0 else Fraction(j * hs).limit_denominator(10 ** 15)
        if j:
            for _ in range(8):                                  # Newton on F'(t) = 0
                _, F1, F2 = P.derivs(t)
                f1 = (Fraction(F1.lo) + Fraction(F1.hi)) / 2
                f2 = (Fraction(F2.lo) + Fraction(F2.hi)) / 2
                if f2 >= 0:
                    break
                t = (t - f1 / f2).limit_denominator(10 ** 40)
        F0, _, _ = P.derivs(t)
        out.append((t, isqrt_lo(F0.flo()), isqrt_hi(F0.fhi())))
    out.sort(key=lambda x: -x[1])
    return out


# ------------------------------------------------------------------------------------------------------ one row
def enclose(K, alpha, out=print, label="", delta=Fraction(1, 10 ** 24), N=N_TRUNC):
    """A rigorous enclosure [A_lo, A_hi] of A(K, e^(i alpha)) = max |P|; returns a dict."""
    t0 = time.time()
    P = Theta(K, alpha, N)
    maxima = find_maxima(P)
    tbest, lo_best, hi_best = maxima[0]
    A_N_lo = Fraction(lo_best)
    target_A = A_N_lo + delta                                   # we prove |P_N| <= A_N_lo + delta everywhere
    res = cover(P, target_A ** 2)
    res["ok"] = res["ok"] and range_ok(COVER_END)
    A_lo = A_N_lo - P.tail
    A_hi = target_A + P.tail
    out(f"  {label}K = {K} = {float(K)}, alpha = {alpha} = {float(alpha)}; N = {N}, tail |R_N| <= {float(P.tail):.3e}, "
        f"G3 = {float(P.G3):.2f}")
    for t, l, h in maxima[:3]:
        out(f"    local maximum of |P_N| near theta = {float(t):.15f}: |P_N| in [{l:.32f}, {h:.32f}]")
    if res["ok"]:
        out(f"    cover of [0, {COVER_END}] (> pi/2): {res['leaves']:,} cells accepted ({res['nodes']:,} evaluated, "
            f"depth up to {res['deepest']}, smallest width {float(res['smallest']):.2e}); "
            f"largest cell bound for |P_N| {float(isqrt_hi(res['worst'])):.25f}")
    else:
        out(f"    cover FAILED: {res}")
    out(f"    A = max |P| in [{dec(A_lo, 30, 'lo')}, {dec(A_hi, 30, 'hi')}]  (width {float(A_hi - A_lo):.1e}; "
        f"{time.time() - t0:.1f} s)")
    return dict(ok=res["ok"], A_lo=A_lo, A_hi=A_hi, P=P, maxima=maxima, cover=res)


def beta_text(e):
    blo, bhi = 1 / e["A_hi"], 1 / e["A_lo"]
    return blo, bhi, f"beta = 1/A in [{dec(blo, 30, 'lo')}, {dec(bhi, 30, 'hi')}]"


def dec(x, digits, mode):
    """x rounded down (lo) or up (hi) to `digits` decimals, as a string."""
    x = Fraction(x)
    s = 10 ** digits
    v = (x.numerator * s) // x.denominator if mode == "lo" else -((-x.numerator * s) // x.denominator)
    sign = "-" if v < 0 else ""
    v = abs(v)
    return f"{sign}{v // s}.{str(v % s).zfill(digits)}"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def row_r2026b(certs, controls):
    path = os.path.join(certs, CERT[0])
    digest = sha256(path)
    print(f"Row r2026b ([R2026b], Roehrig, PR #184): {CERT[0]}, sha256 {digest}")
    if digest != CERT[1]:
        print(f"  SHA-256 differs from the pinned {CERT[1]}")
        print("RESULT: r2026b FAIL")
        return False
    d = json.load(open(path))
    K, alpha, claim = Fraction(d["K"]), Fraction(d["alpha"]), Fraction(d["lower_bound"])
    print(f"  parameters read as exact decimals: K = {d['K']}, alpha = {d['alpha']}; claim C_51 >= {d['lower_bound']}")
    e = enclose(K, alpha)
    ok = e["ok"]
    blo, bhi, txt = beta_text(e)
    print(f"  {txt}")
    imp = claim < blo
    print(f"  claim {d['lower_bound']} < 1/A_hi: {'implied' if imp else 'NOT implied'} (margin {float(blo - claim):.4e})")
    ok &= imp
    # the two local maxima of |P_N| on [0, pi/2]: theta = 0 and an interior one
    mx = e["maxima"]
    at0 = [m for m in mx if m[0] == 0]
    inner = [m for m in mx if m[0] != 0]
    if at0 and inner:
        glo = Fraction(inner[0][1]) - Fraction(at0[0][2])
        ghi = Fraction(inner[0][2]) - Fraction(at0[0][1])
        print(f"  |P_N| at the interior local maximum (theta = {dec(inner[0][0], 15, 'lo')}...) minus |P_N(0)| lies in "
              f"[{float(glo):.6e}, {float(ghi):.6e}]: the global maximum is "
              f"{'interior' if glo > 0 else ('at theta = 0' if ghi < 0 else 'undecided')}")
    # the note's printed enclosures
    their_hi = Fraction("1.709171425067816128872617735945")
    their_recip_lo = Fraction("0.585078819674464320841496314106")
    print(f"  Roehrig's upper bound 1.709171425067816128872617735945 >= our A_lo: {e['A_lo'] <= their_hi}; "
          f"it exceeds our A_hi by {float(their_hi - e['A_hi']):.3e}")
    print(f"  Roehrig's reciprocal 0.585078819674464320841496314106 lies in our beta enclosure or below it: "
          f"{their_recip_lo <= bhi}")
    ok &= e["A_lo"] <= their_hi and their_recip_lo <= bhi
    if controls:
        print("Red controls:")
        one = claim + Fraction(1, 10 ** 13)
        print(f"  (not a control) the claim one unit higher, {dec(one, 13, 'lo')}, is "
              f"{'also implied: the printed bound is not tight' if one < blo else 'not implied'}")
        up = Fraction(-((-bhi.numerator * 10 ** 13) // bhi.denominator), 10 ** 13)   # least 13-decimal number >= bhi
        refuted = up > bhi
        print(f"  C1 the least 13-decimal claim above our enclosure, {dec(up, 13, 'lo')}: "
              f"{'false (above 1/A_lo), as expected' if refuted else 'NOT REFUTED (control failed)'}")
        ok &= refuted
        Pold = Theta(Fraction("3.568182317714"), Fraction("3.961543335688"))
        F0, _, _ = Pold.derivs(Fraction(0))
        bold_hi = 1 / (Fraction(isqrt_lo(F0.flo())) - Pold.tail)
        c2 = bold_hi < claim
        print(f"  C2 the claim with the previous parameters (K = 3.568182317714, alpha = 3.961543335688): "
              f"beta <= 1/(|P(0)| - tail) = {dec(bold_hi, 22, 'hi')} < claim: {'refused, as expected' if c2 else 'ACCEPTED (control failed)'}")
        ok &= c2
        e3 = enclose(K, alpha + Fraction(1, 10 ** 6), out=lambda *a, **k: None)
        b3 = 1 / e3["A_lo"]
        c3 = b3 < claim
        print(f"  C3 alpha raised by 10^-6: beta <= {dec(b3, 16, 'hi')} < claim: {'refused, as expected' if c3 else 'ACCEPTED (control failed)'}")
        ok &= c3
        e4 = enclose(K, alpha, out=lambda *a, **k: None, N=2)
        c4 = 1 / e4["A_hi"] <= claim
        print(f"  C4 truncation at N = 2 (tail bound {float(e4['P'].tail):.3e}): the best bound it gives, "
              f"{dec(1 / e4['A_hi'], 16, 'lo')}, does not reach the claim: {'refused, as expected' if c4 else 'ACCEPTED (control failed)'}")
        ok &= c4
        c5 = not range_ok(Fraction(3, 2)) and range_ok(COVER_END)
        print(f"  C5 a cover of [0, 3/2] only: the range check (end >= pi/2) refuses it and accepts [0, {COVER_END}]: "
              f"{'as expected' if c5 else 'CONTROL FAILED'}")
        ok &= c5
        low = Fraction(e["maxima"][0][1]) - Fraction(1, 10 ** 20)
        r6 = cover(e["P"], low ** 2)
        c6 = not r6["ok"]
        where = f" (stopped at the cell [{float(r6['fail'][0]):.12f}, {float(r6['fail'][1]):.12f}])" if "fail" in r6 else ""
        print(f"  C6 the cover with |P_N| <= (our lower bound - 10^-20) as target: "
              f"{'refused, as expected' + where if c6 else 'ACCEPTED (control failed)'}")
        ok &= c6
    print(f"RESULT: r2026b {'PASS' if ok else 'FAIL'}")
    return ok


def row_s2026(controls):
    K, alpha = Fraction("3.568182317714"), Fraction("3.961543335688")
    print("Row s2026 (the previous exact parameters, as Roehrig's note and the PR state them)")
    e = enclose(K, alpha)
    blo, bhi, txt = beta_text(e)
    print(f"  {txt}")
    s_claim = Fraction("0.585078819653")
    r_bound = Fraction("0.585078819653352950854632621542")
    a = s_claim < blo
    b = bhi < r_bound
    print(f"  [S2026]'s 0.585078819653 < beta: {'yes' if a else 'NO'} (margin {float(blo - s_claim):.4e})")
    print(f"  Roehrig's beta(f_old) < 0.585078819653352950854632621542: {'yes' if b else 'NO'} "
          f"(our upper end is {float(r_bound - bhi):.3e} below it)")
    ok = e["ok"] and a and b
    print(f"RESULT: s2026 {'PASS' if ok else 'FAIL'}")
    return ok


def row_ht2026(controls):
    K, alpha = Fraction(7137, 2000), Fraction(198074929, 50000000)
    print("Row ht2026 (He and Tang's parameters K0 = 7137/2000, alpha0 = 198074929/50000000, Theorem 1.2)")
    e = enclose(K, alpha)
    blo, bhi, txt = beta_text(e)
    print(f"  {txt}")
    a = e["A_hi"] < Fraction(170919, 100000) and Fraction(58507, 100000) < blo
    print(f"  A0 < 1.70919 (their Proposition 3.5) and beta > 0.58507 (Theorem 1.2): {'yes' if a else 'NO'}")
    # their mesh: M = 2,000,000, theta_j = 2 pi j / M, P truncated at n = 5; maximum upper bound printed for j = 1794536
    pi = pi_interval()
    j, M = 1794536, 2000000
    # |P(theta_j)| = |P(2 pi - theta_j)|; 2 pi - theta_j = 2 pi (M - j) / M lies in [0, pi/2]
    th = (pi.scale(Fraction(2 * (M - j), M)))
    P5 = Theta(K, alpha, N=5)
    vals = []
    for t in (th.flo(), th.fhi()):
        F0, F1, _ = P5.derivs(t)
        vals.append((F0, F1))
    # theta_j lies in [th.lo, th.hi]; there F = |P_5|^2 is within |F'(th.lo)| w + G2 w^2 / 2 of F(th.lo) (w the width)
    w = th.fhi() - th.flo()
    d1 = max(abs(vals[0][1].flo()), abs(vals[0][1].fhi()))
    fl = vals[0][0].flo() - d1 * w - P5.G2 * w * w / 2
    fh = vals[0][0].fhi() + d1 * w + P5.G2 * w * w / 2
    mlo, mhi = isqrt_lo(fl), isqrt_hi(fh)
    their = Fraction("1.709176397414741464507875")
    agree = Fraction(mhi) <= their + Fraction(1, 10 ** 24) and their - Fraction(mlo) < Fraction(1, 10 ** 23)
    print(f"  their mesh point j = {j} (theta_j = 2 pi j / M, M = 2,000,000; N = 5): |P_5(theta_j)| in "
          f"[{mlo:.30f}, {mhi:.30f}]; their printed upper bound 1.709176397414741464507875: "
          f"{'agrees' if agree else 'DISAGREES'}")
    ok = e["ok"] and a and agree
    print(f"RESULT: ht2026 {'PASS' if ok else 'FAIL'}")
    return ok


def selftest():
    print("Self-tests")
    ok = True
    s, c = sincos(Fraction(1))
    t1 = Fraction(s.lo) < Fraction(math.sin(1)) + Fraction(1, 10 ** 15) and Fraction(c.hi) > Fraction(math.cos(1)) - Fraction(1, 10 ** 15)
    pi = pi_interval()
    t2 = Fraction(pi.lo) < Fraction("3.14159265358979323846264338327950288419716939937510582097494459231") < Fraction(pi.hi)
    print(f"  sin 1, cos 1 by our series agree with libm; pi by Machin encloses the 66-digit value: {t1 and t2}")
    ok &= t1 and t2
    # derivatives against finite differences, the symmetries, and the cell bound against dense sampling
    P = Theta(Fraction("3.5681854310698435"), Fraction("3.9615431126097391"))
    rnd = random.Random(51)
    worst_fd = 0.0
    for _ in range(20):
        t = Fraction(rnd.randrange(1, 10 ** 6), 10 ** 6) * Fraction(3, 2)
        F0, F1, F2 = P.derivs(t)
        hh = Fraction(1, 10 ** 8)
        Fp, _, _ = P.derivs(t + hh)
        Fm, _, _ = P.derivs(t - hh)
        fd1 = (Fp.flo() - Fm.flo()) / (2 * hh)
        fd2 = (Fp.flo() - 2 * F0.flo() + Fm.flo()) / hh ** 2
        worst_fd = max(worst_fd, float(abs(fd1 - F1.flo())), float(abs(fd2 - F2.flo())) * 1e-8)
    t3 = worst_fd < 1e-6
    print(f"  F' and F'' against central differences at 20 random points: largest discrepancy {worst_fd:.2e}: {t3}")
    ok &= t3
    sym = 0.0
    for _ in range(20):
        t = rnd.random() * 1.5
        sym = max(sym, abs(P.float_abs(t) - P.float_abs(math.pi - t)), abs(P.float_abs(t) - P.float_abs(-t)))
    t4 = sym < 1e-12
    print(f"  |P(t)| = |P(-t)| = |P(pi - t)| at 20 random points (floating point): largest difference {sym:.1e}: {t4}")
    ok &= t4
    viol = 0
    for _ in range(300):
        a = Fraction(rnd.randrange(0, 15000), 10000)
        b = a + Fraction(rnd.randrange(1, 2000), 10 ** 5)
        ub, _ = cell_bound(P, a, b)
        for k in range(21):
            t = float(a + (b - a) * k / 20)
            if P.float_abs(t) ** 2 > float(ub) + 1e-12:
                viol += 1
    t5 = viol == 0
    print(f"  300 random cells: the cell bound is above |P_N|^2 at 21 points of each: {t5}")
    ok &= t5
    # the tail bound against a long direct sum
    direct = sum(2 * Fraction(P.K) ** -(n * (n + 1) // 2) for n in range(N_TRUNC + 1, N_TRUNC + 40))
    t6 = direct <= P.tail and P.tail < direct * (1 + Fraction(1, 10 ** 6))
    print(f"  the tail bound {float(P.tail):.6e} against the sum of the next 39 terms {float(direct):.6e}: {t6}")
    ok &= t6
    print(f"RESULT: selftest {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    what = args[0] if args else "all"
    certs = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    controls = "--controls" in sys.argv
    t0 = time.time()
    ok = True
    if what in ("selftest", "all"):
        ok &= selftest()
    if what in ("r2026b", "all"):
        ok &= row_r2026b(certs, controls)
    if what in ("s2026", "all"):
        ok &= row_s2026(controls)
    if what in ("ht2026", "all"):
        ok &= row_ht2026(controls)
    print(f"run time {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
