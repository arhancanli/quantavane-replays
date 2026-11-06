"""Independent check of the C_42 upper bound 0.6873358063 [Kos2026] (Turan's pure power-sum constant, constants/42a.md,
row added in teorth/optimizationproblems PR #263), from the criterion as section 4 of Koshta's note states it for any
first-block cut-off 0 < tau < 1/2, and the exact inputs of the certificate alone: the checker/, verifier/ and
crosscheck/ code of DennisKoshta/c42-turan-power-sum-certificate was neither read nor run.

  ./fetch.sh && python3 -I kos_check.py [CERTS_ROOT] [--controls] [--workers=N]   (default 4 worker processes)

The criterion. With tau = t_0 < t_1 < ... < t_M = 1 - tau, alpha = a + ib (a > 0), s = 1 - alpha, middle values
eta_0, ..., eta_{M-1}, delta_0 = eta_0 - s, delta_j = eta_j - eta_{j-1}, and d* = ceil(1/tau) - 1 (the largest d with
d tau < 1, found in exact rationals),
  I_A(x) = int_0^x v^(A-1) / (1 - v) dv,   K = I_alpha(tau),   D = I_a(tau),
  T_d(u) = int over {x_i >= u_i, x_1 + ... + x_d <= 1} of (1 - sum x)^(alpha-1) / (x_1 ... x_d) dx  (0 if sum u >= 1),
  Y = 1 + s K - sum_j delta_j (I_alpha(1 - t_j) - K)
        + sum_(d=2..d*) ((-1)^d / d!) sum over ordered d-tuples of delta_j1 ... delta_jd T_d(t_j1, ..., t_jd),
and if |s| < C, every |eta_j| < C and |Y| < C D, then C_42 <= C. For tau >= 1/3 this is the criterion of Griego and
Pinho that c42_check.py replays (d* = 2, T_2 = J).

How T_d is evaluated: our own series, not the note's sections 5 and 7. Let G_d(S; u) be the integral of
1/(x_1 ... x_d) over the slice {x >= u, sum x = S} and V_m(c; a) the integral of prod 1/x_i over {x >= a, sum x <= c}.
On the slice, 1/(x_1 ... x_d) = (1/S) sum_k prod_(i != k) 1/x_i, so
  G_d(S; u) = (1/S) sum_k V_(d-1)(S - u_k; u without u_k),   V_m(c; a) = int_(sum a)^c G_m(c'; a) dc',
  V_1(c; a) = log(c / a),   T_d(u) = int_0^L v^(alpha-1) G_d(1 - v; u) dv   (L = 1 - sum u).
Every sub-tuple a of u is expanded about c = L + sum a in one variable w in [0, L]: its V and G series follow from
those of the smaller sub-tuples by adding, dividing by (L + sum a - w) and integrating, and the constant
V(L + sum a; a) is the integral of its G series over [0, L]. Then T_d(u) = sum_n g_n L^(alpha+n) / (alpha+n). All
these series converge for |w| < L + min u, which exceeds L for every tuple, so no condition such as L < min u is
needed; a majorant G (1 - w/(L + min u))^-2, carried through the same recursion, bounds every tail. Ratios for
Koshta's breaks: at most 0.742 (T_2), 0.653 (T_3), 0.468 (T_4). The T_d series are cut where their proved tail is below
about 10^-44 (the bounds actually added are printed), the I series as in c42_check.py (10^-42). The arithmetic is
c42_check.py's decimal intervals at 50 digits, rounded outward, with x^alpha = e^(a ln x) (cos(b ln x) + i sin(b ln x))
and the angle reduced by an enclosure of 2 pi from Machin's formula.

--controls: (a) the same code on the tau >= 1/3 certificates of Griego and Pinho (d* = 2) against c42_check.py;
(b) C lowered just below our |Y|/D, and eta_0 replaced by s, must FAIL; (c) T_2 against c42_check.py's J on every pair;
T_3 and T_4 against the note's Dirichlet series about the vertex x = u (rigorous, where q = L / min u < 1); T_3 with
q >= 1 against the integral int_0^L H_2(sigma) K(L - sigma) d sigma and T_4 against the scaling identity
T_4(u) = int (1-x)^(alpha-1) T_3(u'/(1-x)) dx/x, both by tanh-sinh quadrature (floating point, not rigorous).
Python 3.9+, standard library only."""
import decimal
import json
import math
import multiprocessing
import os
import re
import sys
import time
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)            # for c42_check.py next to this file (python3 -I leaves the script's folder out)
import c42_check as base            # noqa: E402
from c42_check import I, Z, q, iln, iexp, isincos, zq, zinv, I_series, ipow, LO, HI, NEAR, EPS_DIGITS  # noqa: E402

ZERO, ONE = Decimal(0), Decimal(1)
CERT = os.path.join("koshta", "c42_tau41_200_M64_0.6873358063.json")
# The claimant's printed values (NOTES.md section 6), read only after our own computation, for the comparison lines.
CLAIM_RATIO = Decimal("0.687335806186516616249961526259")      # printed with "+- 1e-31"
CLAIM_MARGIN = "4.348178358e-11"
CLAIM_S4 = "1.214e-5 - 3.45e-6i"
CLAIM_COUNTS = "1056 T_2 pairs, 2000 T_3 by the series of its section 7.1 (q <= 4/5), 443 by its integral, 609 T_4"

# Every rigorous computation below uses only LO/HI context methods, copy_negate and copy_abs on Decimals: the plain
# operators (+, -, *, /, abs, unary minus, sum) would round in decimal's default context (28 digits, half-even).


# ---- arithmetic added to c42_check.py's ---------------------------------------------------------------------------
def _pi():
    """pi = 16 arctan(1/5) - 4 arctan(1/239) (Machin). Each arctan is an alternating series with decreasing terms, so
    it lies between consecutive partial sums: exact rationals, then rounded outward."""
    def atan_inv(n, terms):
        s = sum(Fraction((-1) ** k, (2 * k + 1) * n ** (2 * k + 1)) for k in range(terms))
        nxt = Fraction(1, (2 * terms + 1) * n ** (2 * terms + 1))
        return (s, s + nxt) if terms % 2 == 0 else (s - nxt, s)
    a, b = atan_inv(5, 60), atan_inv(239, 20)
    lo, hi = 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]
    return I(q(lo).lo, q(hi).hi)


TWO_PI = _pi() * I(2)


def isincos_any(x):
    """(sin x, cos x) for any interval x: x - 2 pi k with the rigorous 2 pi, then c42_check.isincos (needs |x| < 8)."""
    k = round(float(x.lo) / (2 * math.pi))
    return isincos(x - TWO_PI * I(k) if k else x)


def cpow_any(x, ar, ai):
    """x^(ar + i ai) for a rational 0 < x < 1 (real logarithm), for any size of ai ln x."""
    lx = iln(q(x))
    mod = iexp(q(ar) * lx)
    s, c = isincos_any(q(ai) * lx)
    return Z(mod * c, mod * s)


def mulp(xl, xh, pl, ph):
    """[xl, xh] [pl, ph] for pl > 0, rounded outward."""
    return LO.multiply(xl, pl if xl >= 0 else ph), HI.multiply(xh, ph if xh >= 0 else pl)


def mulg(xl, xh, yl, yh):
    """[xl, xh] [yl, yh] for any signs, rounded outward."""
    if yl >= 0:
        return LO.multiply(xl, yl if xl >= 0 else yh), HI.multiply(xh, yh if xh >= 0 else yl)
    if yh <= 0:
        return LO.multiply(xh, yl if xh >= 0 else yh), HI.multiply(xl, yh if xl >= 0 else yl)
    return (min(LO.multiply(xl, yh), LO.multiply(xh, yl)), max(HI.multiply(xl, yl), HI.multiply(xh, yh)))


def abs2_lo(z):
    """A lower bound for |z|^2 over the rectangle (its distance from 0)."""
    def dist(v):
        return v.lo if v.lo > 0 else (v.hi.copy_abs() if v.hi < 0 else ZERO)
    x, y = dist(z.re), dist(z.im)
    return LO.add(LO.multiply(x, x), LO.multiply(y, y))


def isqrt(x_lo, x_hi):
    """[sqrt x_lo, sqrt x_hi] rounded outward (decimal's sqrt is correctly rounded: one unit either way, checked)."""
    lo = NEAR.next_minus(NEAR.sqrt(x_lo)) if x_lo > 0 else ZERO
    hi = NEAR.next_plus(NEAR.sqrt(x_hi))
    assert HI.multiply(lo, lo) <= x_lo and LO.multiply(hi, hi) >= x_hi
    return lo, hi


def overlap(x, y):
    return x.lo <= y.hi and y.lo <= x.hi


def zoverlap(x, y):
    return overlap(x.re, y.re) and overlap(x.im, y.im)


def mid(v):
    with decimal.localcontext(NEAR):
        return (v.lo + v.hi) / 2


def width(z):
    return max(HI.subtract(z.re.hi, z.re.lo), HI.subtract(z.im.hi, z.im.lo))


_INV = {}


def inv_alpha(ar, ai, N):
    """[1/(alpha + n) for n < N] as rectangles, kept per process."""
    lst = _INV.setdefault((ar, ai), [])
    while len(lst) < N:
        lst.append(zinv(ar + len(lst), ai))
    return lst


# ---- T_d by the series about the face sum x = 1 ----------------------------------------------------------------------
def n_terms(L, umin):
    """Series length for one L: the proved tail G L^a kappa r^N / (1 - r), r = L / (L + min u), falls below
    10^-(EPS_DIGITS + 2) if G L^a kappa <= 1000 (the bound actually added uses the true G and is reported)."""
    r = float(L / (L + umin))
    return max(8, int(math.ceil(((EPS_DIGITS + 2) * math.log(10) + math.log(1e3 / (1 - r))) / -math.log(r))) + 1)


class Face:
    """T_d(u) for every tuple u with one L = 1 - sum u. For a sub-tuple a of u (a sorted tuple of rationals), with
    c_a = L + sum a and w in [0, L]:
      V_a(w) = V_|a|(c_a - w; a) = sum_n V_a,n w^n,   g_a(w) = G_|a|(c_a - w; a) = sum_n g_a,n w^n,
      g_a = (sum over k of V_(a - a_k)) / (c_a - w),   V_a,n = -g_a,(n-1) / n (n >= 1),   V_a,0 = int_0^L g_a,
      V_(x) = log((L + x - w) / x) for a single x,     and T_d(u) = int_0^L v^(alpha-1) g_u(v) dv  (c_u = 1).
    Majorants, with R_a = L + min a: V_a << MV_a (1 - w/R_a)^-1 and g_a << G_a (1 - w/R_a)^-2, where
    MV_(x) = max(log((L+x)/x), 1), G_a = (sum_k MV_(a - a_k)) / c_a (as 1/(c_a - w) << (1/c_a)(1 - w/R_a)^-1 and
    R_a <= R_(a - a_k)), and MV_a = max(|V_a,0|, G_a R_a). So |g_a,n| <= G_a (n+1) R_a^-n, and the tails are
    |sum_(n>=N) g_a,n L^(n+1)/(n+1)| <= G_a L r^N/(1-r) and
    |sum_(n>=N) g_u,n L^(alpha+n)/(alpha+n)| <= G_u L^a kappa r^N/(1-r), r = L/R, kappa = max(1, (N+1)/(N+a))."""

    def __init__(self, L, ar, ai, N):
        self.L, self.N, self.ar, self.ai = Fraction(L), N, Fraction(ar), Fraction(ai)
        assert 0 < self.L < 1
        Lq = q(self.L)
        self.Plo, self.Phi = [], []                       # L^(n+1) / (n+1)
        pw = Lq
        for n in range(N):
            self.Plo.append(LO.divide(pw.lo, n + 1))
            self.Phi.append(HI.divide(pw.hi, n + 1))
            pw = pw * Lq
        La = cpow_any(self.L, self.ar, self.ai)
        inv = inv_alpha(self.ar, self.ai, N)
        self.W = []                                       # L^(alpha+n) / (alpha+n)
        Lm = I(1)
        for n in range(N):
            w = La * Lm * inv[n]
            self.W.append((w.re.lo, w.re.hi, w.im.lo, w.im.hi))
            Lm = Lm * Lq
        self.La_hi = iexp(q(self.ar) * iln(Lq)).hi        # L^a from above
        self.memo = {}
        self.max_tail = ZERO                              # the largest tail added to a constant V_a,0

    def V(self, a):
        hit = self.memo.get(a)
        if hit:
            return hit
        L, N = self.L, self.N
        if len(a) == 1:
            C = L + a[0]
            const = iln(q(C / a[0]))
            iC = q(1 / C)
            lo, hi = [const.lo], [const.hi]
            pl = ph = ONE
            for n in range(1, N):                         # -1 / (n C^n)
                pl, ph = LO.multiply(pl, iC.lo), HI.multiply(ph, iC.hi)
                lo.append(LO.divide(ph.copy_negate(), n))
                hi.append(HI.divide(pl.copy_negate(), n))
            MV = max(const.hi, ONE)
        else:
            glo, ghi, G, R = self.g(a)
            sl = sh = ZERO
            for n in range(N):
                x, y = mulp(glo[n], ghi[n], self.Plo[n], self.Phi[n])
                sl, sh = LO.add(sl, x), HI.add(sh, y)
            r = L / R
            tail = HI.divide(HI.multiply(HI.multiply(G, q(L).hi), ipow(q(r), N).hi), q(1 - r).lo)
            self.max_tail = max(self.max_tail, tail)
            sl, sh = LO.subtract(sl, tail), HI.add(sh, tail)
            lo, hi = [sl], [sh]
            for n in range(1, N):
                lo.append(LO.divide(ghi[n - 1].copy_negate(), n))
                hi.append(HI.divide(glo[n - 1].copy_negate(), n))
            MV = max(sl.copy_abs(), sh.copy_abs(), HI.multiply(G, q(R).hi))
        self.memo[a] = (lo, hi, MV)
        return lo, hi, MV

    def g(self, a):
        L, N = self.L, self.N
        c = L + sum(a)
        R = L + a[0]
        slo, shi = [ZERO] * N, [ZERO] * N
        GM = ZERO
        for v, m in Counter(a).items():
            sub = list(a)
            sub.remove(v)
            vlo, vhi, MV = self.V(tuple(sub))
            if m == 1:
                slo = [LO.add(x, y) for x, y in zip(slo, vlo)]
                shi = [HI.add(x, y) for x, y in zip(shi, vhi)]
            else:
                mm = Decimal(m)
                slo = [LO.add(x, LO.multiply(y, mm)) for x, y in zip(slo, vlo)]
                shi = [HI.add(x, HI.multiply(y, mm)) for x, y in zip(shi, vhi)]
            GM = HI.add(GM, HI.multiply(MV, Decimal(m)))
        glo, ghi = [], []
        gl = gh = ZERO
        if c == 1:
            for n in range(N):
                gl, gh = LO.add(gl, slo[n]), HI.add(gh, shi[n])
                glo.append(gl)
                ghi.append(gh)
        else:
            ic = q(1 / c)
            for n in range(N):
                gl, gh = mulp(LO.add(slo[n], gl), HI.add(shi[n], gh), ic.lo, ic.hi)
                glo.append(gl)
                ghi.append(gh)
        return glo, ghi, HI.divide(GM, q(c).lo), R

    def T(self, u):
        """(T_d(u) as a rectangle, the tail bound added) for a sorted tuple u with sum u = 1 - L."""
        assert self.L + sum(u) == 1
        glo, ghi, G, R = self.g(u)
        rl = rh = il = ih = ZERO
        for n in range(self.N):
            wrl, wrh, wil, wih = self.W[n]
            a1, b1 = mulg(glo[n], ghi[n], wrl, wrh)
            a2, b2 = mulg(glo[n], ghi[n], wil, wih)
            rl, rh, il, ih = LO.add(rl, a1), HI.add(rh, b1), LO.add(il, a2), HI.add(ih, b2)
        r, N = self.L / R, self.N
        kappa = max(ONE, HI.divide(Decimal(N + 1), LO.add(Decimal(N), q(self.ar).lo)))
        tail = HI.multiply(HI.multiply(HI.multiply(G, self.La_hi), kappa), ipow(q(r), N).hi)
        tail = HI.divide(tail, q(1 - r).lo)
        return Z(I(rl, rh).widen(tail), I(il, ih).widen(tail)), tail


# ---- T_d by the Dirichlet series about the vertex x = u (the note's sections 5 and 7.1; controls only) ---------------
def T_vertex(u, ar, ai):
    """T_d(u) = L^(alpha+d-1) sum_N (-L)^N h_N / (alpha)_(N+d), h_N = sum_(|k|=N) prod k_i! / u_i^(k_i+1), for
    q = L / min u < 1 (checked exactly). The tail from N0 on is at most B_N0 / (1 - q), with
    B_N = L^(a+d-1) q^N C(N+d-1, d-1) N! / (prod u (a)_(N+d)) (h_N <= C(N+d-1, d-1) N! min(u)^-N / prod u and
    |alpha + j| >= a + j; B_(N+1)/B_N <= q). Returns (rectangle, N0, tail)."""
    d, L = len(u), 1 - sum(u)
    qq, a = L / min(u), Fraction(ar)
    assert 0 < qq < 1
    pu = Fraction(1)
    for x in u:
        pu *= x
    # N0: B_N0 / (1 - q) below 10^-(EPS_DIGITS + 2), found in floating point and then bounded rigorously below
    lb = (a + d - 1) * math.log(L) - math.log(pu) - math.log(1 - qq) + math.lgamma(float(a))
    N0 = 10
    while lb + N0 * math.log(qq) + math.log(comb(N0 + d - 1, d - 1)) + math.lgamma(N0 + 1) \
            - math.lgamma(float(a) + N0 + d) > -(EPS_DIGITS + 2) * math.log(10):
        N0 += 10
    e = None                                              # e_N = L^N h_N: convolution of k! (L/u_i)^k / u_i
    for x in u:
        r = q(L / x)
        cur = q(1 / x)
        bl, bh = [cur.lo], [cur.hi]
        for k in range(1, N0):
            cur = cur * r * I(k)
            bl.append(cur.lo)
            bh.append(cur.hi)
        if e is None:
            e = (bl, bh)
            continue
        el, eh = e
        nl, nh = [], []
        for n in range(N0):
            sl = sh = ZERO
            for k in range(n + 1):
                sl = LO.add(sl, LO.multiply(el[k], bl[n - k]))
                sh = HI.add(sh, HI.multiply(eh[k], bh[n - k]))
            nl.append(sl)
            nh.append(sh)
        e = (nl, nh)
    inv = Z(I(1))
    for j in range(d):
        inv = inv * zinv(a + j, ai)
    tot = Z(I(0))
    for n in range(N0):
        x = I(e[0][n], e[1][n])
        tot = tot + inv * (-x if n % 2 else x)
        inv = inv * zinv(a + n + d, ai)
    poch = ONE
    for j in range(N0 + d):
        poch = LO.multiply(poch, q(a + j).lo)
    B = HI.multiply(iexp(q(a + d - 1) * iln(q(L))).hi, ipow(q(qq), N0).hi)
    B = HI.multiply(HI.multiply(B, Decimal(comb(N0 + d - 1, d - 1))), Decimal(factorial(N0)))
    B = HI.divide(HI.divide(B, LO.multiply(poch, q(pu).lo)), q(1 - qq).lo)
    res = cpow_any(L, a + d - 1, ai) * tot
    return Z(res.re.widen(B), res.im.widen(B)), N0, B


# ---- tanh-sinh quadrature (controls only; floating point at 50 digits, not rigorous) ---------------------------------
def tanh_sinh(f, A, B, h=Fraction(1, 32), tmax=Fraction(43, 10)):
    """int_A^B f(x) dx for rational A < B, with x = A + (B - A) y, y = 1 / (1 + e^(-2s)), s = (pi/2) sinh t, on the
    nodes t = k h, |t| <= tmax. f takes a rational x and returns (re, im) Decimals. Returns the sums at steps h and 2h
    (the h-nodes contain the 2h-nodes), each as (re, im)."""
    with decimal.localcontext(NEAR):
        pi = mid(_pi())
        K = int(tmax / h)
        acc = {1: [ZERO, ZERO], 2: [ZERO, ZERO]}
        hd = Decimal(h.numerator) / Decimal(h.denominator)
        for k in range(-K, K + 1):
            t = k * hd
            et = t.exp()
            sh, ch = (et - 1 / et) / 2, (et + 1 / et) / 2
            s = pi / 2 * sh
            es = s.exp()
            e2 = es * es
            w = pi / 4 * ch / ((es + 1 / es) / 2) ** 2
            if k <= 0:
                x = A + (B - A) * Fraction(e2 / (1 + e2))
            else:
                x = B - (B - A) * Fraction(1 / (1 + e2))
            fr, fi = f(x)
            for lvl in (1, 2):
                if k % lvl == 0:
                    acc[lvl][0] += w * fr
                    acc[lvl][1] += w * fi
        span = Decimal((B - A).numerator) / Decimal((B - A).denominator)
        return {lvl: (acc[lvl][0] * span * hd * lvl, acc[lvl][1] * span * hd * lvl) for lvl in (1, 2)}


def h2k_integrand(u, ar, ai):
    """sigma -> H_2(sigma) K(L - sigma), whose integral over [0, L] is T_3(u):
    H_2(sigma) = log((u_1+sigma)(u_2+sigma)/(u_1 u_2)) / (u_1+u_2+sigma) (the inner double integral over p_1 + p_2 =
    sigma) and K(lam) = int_0^lam (lam-p)^(alpha-1)/(u_3+p) dp = (lam^alpha/(u_3+lam)) sum_n y^n/(alpha+n),
    y = lam/(u_3+lam) < 1, tail at most y^R/((a+R)(1-y)) times the prefactor."""
    u1, u2, u3 = u
    L = 1 - u1 - u2 - u3

    def f(sig):
        H = iln(q((u1 + sig) * (u2 + sig) / (u1 * u2))) / q(u1 + u2 + sig)
        lam = L - sig
        y = lam / (u3 + lam)
        R = base.terms_for(y)
        inv = inv_alpha(Fraction(ar), Fraction(ai), R)
        tot, yr = Z(I(0)), I(1)
        for n in range(R):
            tot = tot + inv[n] * yr
            yr = yr * q(y)
        tot = tot.widen(HI.divide(yr.hi, LO.multiply(q(Fraction(ar) + R).lo, q(1 - y).lo)))
        val = cpow_any(lam, ar, ai) * tot * (H * q(1 / (u3 + lam)))
        return mid(val.re), mid(val.im)
    return f, L


def scaling_integrand(u, ar, ai):
    """x -> (1-x)^(alpha-1) T_3(u_1/(1-x), u_2/(1-x), u_3/(1-x)) / x on [u_4, 1 - u_1 - u_2 - u_3], whose integral is
    T_4(u) (the note's scaling identity: prod dx_i/x_i is scale-invariant); T_3 by the face series."""
    ell = 1 - u[0] - u[1] - u[2]

    def f(x):
        w = tuple(sorted(ui / (1 - x) for ui in u[:3]))
        L3 = 1 - sum(w)
        T3, _ = Face(L3, ar, ai, n_terms(L3, w[0])).T(w)
        val = cpow_any(1 - x, Fraction(ar) - 1, ai) * T3 * q(1 / x)
        return mid(val.re), mid(val.im)
    return f, u[3], ell


# ---- reading a certificate -------------------------------------------------------------------------------------------
_DEC = re.compile(r"-?[0-9]+(\.[0-9]+)?\Z")
_FRAC = re.compile(r"([0-9]+)/([0-9]+)\Z")


def rational(x):
    """An exact rational from a JSON string: a terminating decimal, or p/q with integers p, q > 0 (FORMAT2.md, 1)."""
    if not isinstance(x, str):
        raise ValueError(f"not a JSON string: {x!r}")
    if _DEC.match(x):
        return Fraction(x)
    m = _FRAC.match(x)
    if m and int(m.group(1)) > 0 and int(m.group(2)) > 0:
        return Fraction(int(m.group(1)), int(m.group(2)))
    raise ValueError(f"not an exact rational: {x!r}")


def load(path):
    """(name, C, alpha, breaks, etas, raw) from a certificate laid out as FORMAT2.md section 1 states; other keys (such
    as provenance) are not inputs and are ignored."""
    d = json.load(open(path))
    if d.get("constant") != "C_42" or len(d["alpha"]) != 2 or any(len(e) != 2 for e in d["eta"]):
        raise ValueError("not a C_42 certificate in the stated layout")
    return ("Kos2026", rational(d["C"]), (rational(d["alpha"][0]), rational(d["alpha"][1])),
            [rational(b) for b in d["breaks"]], [(rational(e[0]), rational(e[1])) for e in d["eta"]], d)


# ---- the criterion ---------------------------------------------------------------------------------------------------
def dstar_of(tau):
    """The largest d with d tau < 1, i.e. ceil(1/tau) - 1, in exact rationals."""
    d = 0
    while (d + 1) * tau < 1:
        d += 1
    return d


def multisets(breaks, M, d):
    """Index tuples i_1 <= ... <= i_d in 0..M-1 with t_i1 + ... + t_id < 1, decided exactly (the breaks increase)."""
    out = []

    def rec(prefix, start, partial):
        k = len(prefix)
        if k == d:
            out.append(tuple(prefix))
            return
        for i in range(start, M):
            if partial + (d - k) * breaks[i] >= 1:
                break
            rec(prefix + [i], i, partial + breaks[i])
    rec([], 0, Fraction(0))
    return out


def _job(job):
    if job[0] == "I":
        _, key, x, ar, ai = job
        return job[0], key, I_series(x, ar, ai)
    _, d, L, idxs, vals, ar, ai = job
    N = max(n_terms(L, v[0]) for v in vals)
    F = Face(L, ar, ai, N)
    res = []
    for ix, v in zip(idxs, vals):
        T, tail = F.T(v)
        res.append((ix, T, tail))
    return job[0], (d, L, N, F.max_tail), res


def run_jobs(jobs, workers):
    if workers <= 1:
        return [_job(j) for j in jobs]
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        return list(pool.imap_unordered(_job, jobs))


def analytic(alpha, breaks, M, dstar, workers=4):
    """K, D, every I_alpha(1 - t_j), and T_d on every multiset of breaks with sum < 1 (2 <= d <= d*)."""
    ar, ai = alpha
    tau = breaks[0]
    jobs = [("I", "K", tau, ar, ai), ("I", "D", tau, ar, 0)] + [("I", j, 1 - breaks[j], ar, ai) for j in range(M)]
    tuples = {d: multisets(breaks, M, d) for d in range(2, dstar + 1)}
    groups = {}
    for d, lst in tuples.items():
        for ix in lst:
            vals = tuple(breaks[i] for i in ix)
            groups.setdefault((d, 1 - sum(vals)), []).append((ix, vals))
    tjobs = [("T", d, L, [x[0] for x in g], [x[1] for x in g], ar, ai) for (d, L), g in groups.items()]
    tjobs.sort(key=lambda j: -len(j[3]) * (2 ** j[1]) * n_terms(j[2], j[4][0][0]))
    out = {"I": {}, "T": {d: {} for d in tuples}, "tuples": tuples, "N": [], "tail_T": ZERO, "tail_V": ZERO}
    for kind, key, res in run_jobs(jobs + tjobs, workers):
        if kind == "I":
            out["I"][key] = res
            continue
        d, L, N, tv = key
        out["N"].append(N)
        out["tail_V"] = max(out["tail_V"], tv)
        for ix, T, tail in res:
            out["T"][d][ix] = T
            out["tail_T"] = max(out["tail_T"], tail)
    out["K"], out["D"] = out["I"]["K"], out["I"]["D"].re
    return out


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def preconditions(C, alpha, breaks, etas):
    """The exact checks, as (what, holds, detail)."""
    ar, ai = alpha
    tau, M = breaks[0], len(etas)
    s = (1 - ar, -ai)
    m_s = C * C - (s[0] ** 2 + s[1] ** 2)
    m_eta = min(C * C - (e[0] ** 2 + e[1] ** 2) for e in etas) if etas else Fraction(-1)
    return [
        ("the cut-off 0 < tau < 1/2", 0 < tau < Fraction(1, 2),
         f"tau = {tau}, d* = {dstar_of(tau) if tau > 0 else '-'}"),
        ("breaks increase from tau to 1 - tau, one more than the middle values",
         M >= 1 and len(breaks) == M + 1 and breaks[-1] == 1 - tau and all(x < y for x, y in zip(breaks, breaks[1:])),
         ""),
        ("a = Re alpha > 0 and 0 < C < 1", ar > 0 and 0 < C < 1, ""),
        ("|s| < C and every |eta_j| < C (exact)", m_s > 0 and m_eta > 0,
         f"least C^2 - |.|^2 = {float(min(m_s, m_eta)):.6e} = {min(m_s, m_eta)}"),
    ]


def combine(C, alpha, breaks, etas, A, dmax=None):
    """Y from the analytic values A and the exact deltas, and the decision |Y|^2 < C^2 D^2; dmax cuts the d-fold
    terms (a control only)."""
    ar, ai = alpha
    M = len(etas)
    s = (1 - ar, -ai)
    deltas = [(etas[0][0] - s[0], etas[0][1] - s[1])] + \
             [(etas[j][0] - etas[j - 1][0], etas[j][1] - etas[j - 1][1]) for j in range(1, M)]
    K, D = A["K"], A["D"]
    Y = Z(I(1)) + zq(*s) * K
    single = Z(I(0))
    for j in range(M):
        single = single + zq(*deltas[j]) * (A["I"][j] - K)
    Y = Y - single
    S = {}
    for d, Td in sorted(A["T"].items()):
        if dmax is not None and d > dmax:
            continue
        tot = Z(I(0))
        for ix, T in sorted(Td.items()):             # a fixed order, so the printed enclosures repeat exactly
            c = (Fraction(factorial(d)), Fraction(0))
            for v in Counter(ix).values():
                c = (c[0] / factorial(v), c[1])
            for i in ix:
                c = cmul(c, deltas[i])
            tot = tot + zq(*c) * T
        S[d] = tot
        Y = Y + zq(Fraction((-1) ** d, factorial(d)), 0) * tot
    y2_hi, y2_lo = Y.abs2_hi(), abs2_lo(Y)
    c2d2 = q(C * C) * D * D
    r_lo, r_hi = isqrt(y2_lo, y2_hi)
    return {"Y": Y, "S": S, "D": D, "K": K, "single": single, "y2_hi": y2_hi, "c2d2": c2d2,
            "ratio": (LO.divide(r_lo, D.hi), HI.divide(r_hi, D.lo)),
            "margin": (LO.subtract(c2d2.lo, y2_hi), HI.subtract(c2d2.hi, y2_lo)),
            "pass": y2_hi < c2d2.lo}


def fmt(z, digits=40):
    return f"Re [{z.re.lo:.{digits}f}, {z.re.hi:.{digits}f}]\n       Im [{z.im.lo:.{digits}f}, {z.im.hi:.{digits}f}]"


def check_certificate(name, C, alpha, breaks, etas, workers, out=print):
    """All checks on one certificate; returns (ok, A, res)."""
    ok = True
    for what, cond, detail in preconditions(C, alpha, breaks, etas):
        ok &= bool(cond)
        out(f"  check {what}: {'PASS' if cond else 'FAIL'}{' (' + detail + ')' if detail else ''}")
    out(f"  (C < {float(base.HARCOS)}, the page's best unstarred upper bound: {'yes' if C < base.HARCOS else 'no'})")
    if not ok:
        return False, None, None
    tau, M = breaks[0], len(etas)
    dstar = dstar_of(tau)
    t0 = time.time()
    A = analytic(alpha, breaks, M, dstar, workers)
    names = {2: "pairs", 3: "triples", 4: "quadruples"}
    parts = []
    for d, lst in sorted(A["tuples"].items()):
        qs = [(1 - sum(breaks[i] for i in ix)) / breaks[ix[0]] for ix in lst]
        extra = ""
        if d >= 3:
            extra = (f" ({sum(1 for x in qs if x <= Fraction(4, 5))} with q = L/min u <= 4/5, "
                     f"{sum(1 for x in qs if x >= 1)} with q >= 1; largest q = {float(max(qs)):.4f})")
        parts.append(f"T_{d} on {len(lst)} {names.get(d, str(d) + '-tuples')}{extra}")
    out(f"  d* = {dstar}: " + "; ".join(parts) + " (multisets of breaks with sum < 1, decided exactly)")
    out(f"  series lengths {min(A['N'])} to {max(A['N'])} terms; largest tail bound added: {A['tail_T']:.2e} (to a T), "
        f"{A['tail_V']:.2e} (to a constant V); {time.time() - t0:.0f} s with {workers} worker(s)")
    res = combine(C, alpha, breaks, etas, A)
    out(f"  K = I_alpha(tau): {fmt(A['K'])}")
    out(f"  D = I_a(tau) in [{res['D'].lo:.40f}, {res['D'].hi:.40f}]")
    out(f"  sum_j delta_j (I_alpha(1 - t_j) - K): {fmt(res['single'])}")
    for d, S in sorted(res["S"].items()):
        out(f"  sum over ordered {d}-tuples of delta...delta T_{d}: {fmt(S)}")
    out(f"  Y in {fmt(res['Y'])}   (width {width(res['Y']):.1e})")
    out(f"  |Y|/D in [{res['ratio'][0]:.40f}, {res['ratio'][1]:.40f}]   "
        f"(width {HI.subtract(res['ratio'][1], res['ratio'][0]):.1e})")
    out(f"  C^2 D^2 - |Y|^2 in [{res['margin'][0]:.15e}, {res['margin'][1]:.15e}]   (C = {float(C):.10f})")
    out(f"  check |Y|^2 < C^2 D^2, so C_42 <= C: {'PASS' if res['pass'] else 'FAIL'}")
    return ok and res["pass"], A, res


# ---- controls --------------------------------------------------------------------------------------------------------
def nearest(items, key, targets):
    """For each target, the item whose key is closest to it (distinct items)."""
    chosen = []
    for t in targets:
        best = min((x for x in items if x not in chosen), key=lambda x: abs(key(x) - t))
        chosen.append(best)
    return chosen


def controls(root, C, alpha, breaks, etas, A, res, workers, out=print):
    ok = True
    ar, ai = alpha
    quiet = lambda *a, **k: None   # noqa: E731

    out("Control (a): the same code at tau >= 1/3 (d* = 2) against c42_check.py on Griego's and Pinho's certificates")
    for sub, fname in (("griego", "turan42_certificate.json"), ("pinho", "c42_256.json")):
        path = os.path.join(root, sub, fname)
        if not os.path.exists(path):
            out(f"  {path} missing (run ./fetch.sh): control not run")
            ok = False
            continue
        nm, C1, al1, br1, et1, _ = base.load(path)
        t0 = time.time()
        good, A1, r1 = check_certificate(nm, C1, al1, br1, et1, workers, out=quiet)
        t1 = time.time()
        old_ok, old = base.evaluate(C1, al1, br1, et1, out=quiet)
        t2 = time.time()
        same = good and old_ok and zoverlap(r1["Y"], old["Y"]) and overlap(r1["D"], old["D"])
        dif = max(HI.subtract(mid(r1["Y"].re), mid(old["Y"].re)).copy_abs(),
                  HI.subtract(mid(r1["Y"].im), mid(old["Y"].im)).copy_abs())
        out(f"  {nm} (tau = {br1[0]}, {len(et1)} block(s), {len(A1['tuples'][2])} pairs): new PASS = {good}, "
            f"|Y|/D in [{r1['ratio'][0]:.28f}, {r1['ratio'][1]:.28f}], C^2D^2 - |Y|^2 >= {r1['margin'][0]:.6e} "
            f"({t1 - t0:.0f} s)")
        out(f"    c42_check.py: PASS = {old_ok}, |Y|/D <= {old['ratio']:.25f} ({t2 - t1:.0f} s); Y and D overlap: "
            f"{'yes' if same else 'NO'}; Y midpoints differ by {dif:.1e}")
        ok &= same

    out("Control (b): altered data must be refused (the T_d values are reused; only C or the deltas change)")
    C_low = Fraction(math.floor(Fraction(res["ratio"][0]) * 10 ** 15), 10 ** 15) - Fraction(1, 10 ** 15)
    r = combine(C_low, alpha, breaks, etas, A)
    pre = [w for w, c, _ in preconditions(C_low, alpha, breaks, etas) if not c]
    out(f"  C lowered to {NEAR.divide(C_low.numerator, C_low.denominator)} (below our |Y|/D): |Y|^2 < C^2 D^2 "
        f"{'holds (control FAILED)' if r['pass'] else 'fails'}"
        f", C^2D^2 - |Y|^2 <= {r['margin'][1]:.3e}; exact checks failing too: {'; '.join(pre) if pre else 'none'}")
    ok &= not r["pass"]
    et2 = [(1 - ar, -ai)] + etas[1:]
    r = combine(C, alpha, breaks, et2, A)
    out(f"  eta_0 replaced by s (the first step removed): |Y|/D >= {r['ratio'][0]:.10f}: "
        f"{'accepted (control FAILED)' if r['pass'] else 'refused, as expected'}")
    ok &= not r["pass"]
    r = combine(C, alpha, breaks, etas, A, dmax=dstar_of(breaks[0]) - 1)
    out(f"  (for scale: without the d* = {dstar_of(breaks[0])} term, |Y|/D would be in [{r['ratio'][0]:.12f}, "
        f"{r['ratio'][1]:.12f}], {'PASS' if r['pass'] else 'FAIL'}; that term moves |Y|/D by "
        f"{HI.subtract(res['ratio'][0], r['ratio'][0]):.3e})")

    out("Control (c): T_d by other methods")
    P = base.Pairs(ar, ai)
    worst, bad = ZERO, 0
    for ix, T in A["T"][2].items():
        J = P.J(breaks[ix[0]], breaks[ix[1]])
        bad += not zoverlap(T, J)
        worst = max(worst, HI.subtract(mid(T.re), mid(J.re)).copy_abs(), HI.subtract(mid(T.im), mid(J.im)).copy_abs())
    verdict = "all overlap" if not bad else f"{bad} DO NOT overlap"
    out(f"  T_2 = J on all {len(A['T'][2])} pairs against c42_check.py's series: {verdict}; largest midpoint "
        f"difference {worst:.1e}")
    ok &= bad == 0

    def qv(ix):
        return float((1 - sum(breaks[i] for i in ix)) / breaks[ix[0]])
    tri = A["tuples"].get(3, [])
    quad = A["tuples"].get(4, [])
    for d, lst, targets in ((3, tri, (0.03, 0.25, 0.5, 0.65, 0.75, 0.85)), (4, quad, (0.05, 0.4, 0.7, 0.9))):
        if not lst:
            continue
        for ix in nearest([x for x in lst if qv(x) < 0.9], qv, targets) if d == 3 else \
                nearest(lst, qv, targets):
            t0 = time.time()
            Tv, N0, B = T_vertex(tuple(breaks[i] for i in ix), ar, ai)
            T = A["T"][d][ix]
            agree = zoverlap(T, Tv)
            ok &= agree
            dif = max(HI.subtract(mid(T.re), mid(Tv.re)).copy_abs(), HI.subtract(mid(T.im), mid(Tv.im)).copy_abs())
            out(f"  T_{d}{ix} (q = {qv(ix):.3f}): vertex series ({N0} terms, tail <= {B:.1e}) and face series overlap: "
                f"{'yes' if agree else 'NO'}; midpoints differ by {dif:.1e}, widths {width(T):.1e} and "
                f"{width(Tv):.1e} ({time.time() - t0:.0f} s)")
    big = [x for x in tri if qv(x) >= 1]
    for ix in (nearest(big, qv, (1.0, 1.4, 1.9)) if big else []):
        t0 = time.time()
        f, L = h2k_integrand(tuple(breaks[i] for i in ix), ar, ai)
        Q = tanh_sinh(f, Fraction(0), L)
        ok &= quad_line(out, f"T_3{ix} (q = {qv(ix):.3f}) by int_0^L H_2 K", A["T"][3][ix], Q, t0)
    for ix in (nearest(quad, qv, (0.9, 0.5)) if quad else []):
        t0 = time.time()
        f, x0, x1 = scaling_integrand(tuple(breaks[i] for i in ix), ar, ai)
        Q = tanh_sinh(f, x0, x1)
        ok &= quad_line(out, f"T_4{ix} (q = {qv(ix):.3f}) by the scaling identity", A["T"][4][ix], Q, t0)
    out(f"  controls: {'all as expected' if ok else 'NOT all as expected'}")
    return ok


def quad_line(out, what, T, Q, t0, tol=Decimal("1e-20")):
    with decimal.localcontext(NEAR):
        d1 = max(abs(Q[1][0] - mid(T.re)), abs(Q[1][1] - mid(T.im)))
        d2 = max(abs(Q[2][0] - Q[1][0]), abs(Q[2][1] - Q[1][1]))
    good = d1 < tol
    out(f"  {what}, tanh-sinh (not rigorous): {Q[1][0]:.6e} {Q[1][1]:+.6e}i differs from the face series by {d1:.1e} "
        f"(step h against 2h: {d2:.1e}); within {tol:.0e}: {'yes' if good else 'NO'} ({time.time() - t0:.0f} s)")
    return good


def compare_claim(A, res, out=print):
    """Our enclosures against the numbers the note prints (section 6)."""
    lo, hi = res["ratio"]
    places = -CLAIM_RATIO.as_tuple().exponent
    unit = Decimal(1).scaleb(-places)
    same = NEAR.quantize(lo, unit) == CLAIM_RATIO == NEAR.quantize(hi, unit)
    inside = lo <= HI.add(CLAIM_RATIO, Decimal("1e-31")) and LO.subtract(CLAIM_RATIO, Decimal("1e-31")) <= hi
    out(f"  the note prints |Y|/D = {CLAIM_RATIO} +- 1e-31: both ends of our enclosure round to it at {places} places: "
        f"{'yes' if same else 'NO'}; within its +- 1e-31: {'yes' if inside else 'NO'}")
    out(f"  the note prints C^2D^2 - |Y|^2 = {CLAIM_MARGIN}: ours is {mid(I(*res['margin'])):.12e}")
    S4 = res["S"].get(4)
    if S4 is not None:
        out(f"  the note prints sum delta^4 T_4 = {CLAIM_S4}: ours is {mid(S4.re):.6e} {mid(S4.im):+.6e}i "
            f"(the four-fold term of Y is this over 24)")
    counts = {d: len(v) for d, v in A["tuples"].items()}
    out(f"  the note counts {CLAIM_COUNTS}; ours: " + ", ".join(f"{n} T_{d}" for d, n in sorted(counts.items())))
    return inside


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    workers = 4
    for a in sys.argv[1:]:
        if a.startswith("--workers="):
            workers = int(a.split("=", 1)[1])
    root = args[0] if args else os.path.join(HERE, "certs")
    path = os.path.join(root, CERT)
    t0 = time.time()
    name, C, alpha, breaks, etas, _ = load(path)
    out = print
    out(f"{name}: {os.path.relpath(path, HERE)}; C = {C} = {float(C):.10f}; alpha = {float(alpha[0])!r} + "
        f"{float(alpha[1])!r}i (exact decimals); tau = {breaks[0]}; {len(etas)} middle blocks")
    ok, A, res = check_certificate(name, C, alpha, breaks, etas, workers)
    if A is not None:
        out("Comparison with the claimant's printed numbers (after our computation; not part of the check):")
        compare_claim(A, res)
    if "--controls" in sys.argv and A is not None:
        ok &= controls(root, C, alpha, breaks, etas, A, res, workers)
    out(f"time {time.time() - t0:.0f} s")
    out(f"RESULT C42 <= {C.numerator}/{C.denominator}: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
