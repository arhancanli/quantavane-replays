#!/usr/bin/env python3
"""83a, program B (written separately from program A): rigorous bounds for C_83 = |lambda_2| (Wirsing's constant).

Same statement as program A: for the polynomial psi of psi.json, prove psi > 0 and lo psi <= V psi <= hi psi on [0, 1],
where V is Wirsing's positive operator
    V psi(x) = sum_{k>=1} [ k/(x+k+1)^2 int_{1/(x+k+1)}^{1/(x+k)} psi + (1+x)/((x+k)^3 (x+k+1)) psi(1/(x+k)) ],
so that lo <= rho(V) = |lambda_2| <= hi (Collatz-Wielandt; rho(V) = |lambda_2| is Wirsing's theorem, not re-proved).
A value from its definition. The self-test checks (T phi)' = -V(phi') against the Gauss-Kuzmin operator itself.

Route B (no partial fractions, no Hurwitz zeta at non-integer points):
  * w = 1/(x + k + 1/2): then 1/(x+k) = 2w/(2-w), 1/(x+k+1) = 2w/(2+w), k = 1/w - 1/2 - x. The k-th summand is a
    power series in w whose coefficients are linear in x, built by series arithmetic (psi(1/(x+k)) by Horner on
    series), exact rationals, to w^S. Tail: the summand is analytic in |w| < 2; on |w| = rho a crude bound M gives
    |coefficient of w^s| <= M rho^-s (Cauchy), and |w| <= 2/3.
  * Summing over k: sum_k w^s = zeta(s, x + 3/2). Around x = 1/2 (t = x - 1/2): zeta(s, 2 + t) =
    sum_m (-1)^m C(s+m-1, m) zeta(s+m, 2) t^m, Lagrange remainder with zeta(., 3/2) <= (2/3)^n (1 + 3/(2(n-1))).
    So V psi(1/2 + t) is ONE polynomial of degree Mt + 1 in t plus a uniform remainder, on |t| <= 1/2.
  * zeta(n, 2) = zeta(n) - 1, zeta(n) = eta(n)/(1 - 2^(1-n)) with eta(n) by the Cohen-Rodriguez Villegas-Zagier
    acceleration of the alternating series (eta is a moment sequence of a positive measure: error <= eta(n)/T_N(3));
    for n > 40 by direct summation with the integral tail.
  * [-1/2, 1/2] cut into cells; the polynomial shifted to each cell centre; a cell passes if
    h_0 - sum |h_j| r^j - remainder > 0.
Arithmetic: decimal intervals (70 digits, ends rounded outward); exact rationals where stated. Standard library only.
Usage: c83a_b.py [--psi psi.json] [--claim-lower Q] [--claim-upper Q] [--perturb J:Q] [--variant wrong-weight]
"""
import sys, json, os
from fractions import Fraction as Q
from decimal import Decimal as D, localcontext, Context, ROUND_FLOOR, ROUND_CEILING
from math import comb

DPREC = 70

def _ctx(rnd):
    return localcontext(Context(prec=DPREC, rounding=rnd))

def qdown(q):
    with _ctx(ROUND_FLOOR):
        return D(q.numerator) / D(q.denominator)
def qup(q):
    with _ctx(ROUND_CEILING):
        return D(q.numerator) / D(q.denominator)

class DI:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        self.lo = lo; self.hi = lo if hi is None else hi
    @staticmethod
    def q(v):
        v = Q(v)
        return DI(qdown(v), qup(v))
    def __add__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        with _ctx(ROUND_FLOOR): lo = self.lo + o.lo
        with _ctx(ROUND_CEILING): hi = self.hi + o.hi
        return DI(lo, hi)
    __radd__ = __add__
    def __neg__(self):                   # copy_negate: exact (unary minus would round to the context)
        return DI(self.hi.copy_negate(), self.lo.copy_negate())
    def __sub__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        return self + (-o)
    def __mul__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        with _ctx(ROUND_FLOOR):
            lo = min(self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        with _ctx(ROUND_CEILING):
            hi = max(self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return DI(lo, hi)
    __rmul__ = __mul__
    def __truediv__(self, o):            # o a positive interval
        if not isinstance(o, DI): o = DI.q(o)
        assert o.lo > 0
        with _ctx(ROUND_FLOOR):
            lo = min(self.lo / o.lo, self.lo / o.hi)
        with _ctx(ROUND_CEILING):
            hi = max(self.hi / o.lo, self.hi / o.hi)
        return DI(lo, hi)
    def mag(self):
        return max(self.lo.copy_abs(), self.hi.copy_abs())

# ---------- series in w with coefficients a + b x (pairs of rationals) ----------
def smul(A, B, S):
    """product of two series, at most one of them depending on x"""
    out = [[Q(0), Q(0)] for _ in range(S + 1)]
    for i, (a0, a1) in enumerate(A):
        if a0 == 0 and a1 == 0: continue
        for j in range(S + 1 - i):
            b0, b1 = B[j]
            if b0 == 0 and b1 == 0: continue
            assert a1 == 0 or b1 == 0
            out[i + j][0] += a0 * b0
            out[i + j][1] += a0 * b1 + a1 * b0
    return out

def summand_series(c, S, variant=None):
    """the k-th summand of V psi as sum_s (alpha_s + beta_s x) w^s, s <= S (exact)"""
    zero = lambda: [[Q(0), Q(0)] for _ in range(S + 1)]
    u = zero(); v = zero()                       # 1/(x+k) = 2w/(2-w) = sum w^(l+1) 2^-l ; 1/(x+k+1) = 2w/(2+w)
    for l in range(S):
        u[l + 1][0] = Q(1, 2 ** l); v[l + 1][0] = Q((-1) ** l, 2 ** l)
    def horner(coefs, base):                     # sum coefs[j] base^j
        acc = zero()
        for cj in reversed(coefs):
            acc = smul(acc, base, S)
            acc[0][0] += cj
        return acc
    d = [cj / (j + 1) for j, cj in enumerate(c)]
    Psi_u = smul(horner(d, u), u, S)             # Psi(t) = sum d_j t^(j+1)
    Psi_v = smul(horner(d, v), v, S)
    integral = [[Psi_u[s][0] - Psi_v[s][0], 0] for s in range(S + 1)]
    # k/(x+k+1)^2 = (1/w - 1/2 - x) v^2 ; v^2 starts at w^2, so (1/w) v^2 is a series
    v2 = smul(v, v, S)
    kfac = zero()
    for s in range(S + 1):
        if s + 1 <= S:
            kfac[s][0] += v2[s + 1][0]            # (1/w) v^2
        kfac[s][0] -= v2[s][0] / 2
        kfac[s][1] -= v2[s][0]
    term1 = smul(kfac, integral, S)
    # (1 + x) u^3 v psi(u)   (variant: u^2 in place of u^3)
    psi_u = horner(c, u)
    p = 2 if variant == "wrong-weight" else 3
    up_ = u
    for _ in range(p - 1):
        up_ = smul(up_, u, S)
    t2 = smul(smul(up_, v, S), psi_u, S)
    term2 = [[a0, a0 + a1] for a0, a1 in t2]      # times (1 + x)
    return [[term1[s][0] + term2[s][0], term1[s][1] + term2[s][1]] for s in range(S + 1)]

def cauchy_M(c, rho):
    """an upper bound for |summand| on |w| = rho (rho < 2), any x in [0, 1]"""
    rho = Q(rho)
    Y = 2 * rho / (2 - rho)                      # |1/(x+k)|, |1/(x+k+1)| <= Y
    Psib = sum(abs(cj) / (j + 1) * Y ** (j + 1) for j, cj in enumerate(c))
    psib = sum(abs(cj) * Y ** j for j, cj in enumerate(c))
    kb = 1 / rho + Q(3, 2)
    return kb * Y * Y * 2 * Psib + 2 * Y ** 3 * Y * psib

# ---------- zeta(n, 2) ----------
def crz_weights(N):
    """Cohen-Rodriguez Villegas-Zagier: sum (-1)^k a_k ~ (sum_k w_k a_k)/d, d = T_N(3)"""
    t0, t1 = 1, 3
    for _ in range(N - 1):
        t0, t1 = t1, 6 * t1 - t0
    d = t1 if N >= 1 else 1
    b = Q(-1); cc = Q(-d); w = []
    for k in range(N):
        cc = b - cc
        w.append(cc)
        b = Q((k + N) * (k - N)) * b / (Q(2 * k + 1, 2) * (k + 1))
    return w, d

CRZ_N = 100
CRZ_W, CRZ_D = crz_weights(CRZ_N)

def eta(n):
    """eta(n) = sum_{k>=0} (-1)^k (k+1)^-n, n >= 1, as an interval"""
    s = DI(D(0))
    for k, wk in enumerate(CRZ_W):
        s = s + DI.q(wk) * DI.q(Q(1, (k + 1) ** n))
    s = s / DI.q(CRZ_D)
    e = Q(1, CRZ_D)                              # |error| <= eta(n)/d <= 1/d
    return s + DI(qdown(-e), qup(e))

def zeta2(n):
    """zeta(n, 2) = sum_{k>=2} k^-n"""
    if n <= 40:
        z = eta(n) / DI.q(1 - Q(1, 2 ** (n - 1)))
        return z - 1
    s = DI(D(0))
    for k in range(2, 101):
        s = s + DI.q(Q(1, k ** n))
    return s + DI(D(0), qup(Q(1, (n - 1) * 100 ** (n - 1))))

def pi_machin():
    def arctan_inv(m):
        s = DI(D(0)); k = 0
        while True:
            term = Q(1, (2 * k + 1) * m ** (2 * k + 1))
            if term < Q(1, 10 ** (DPREC + 5)):
                return s + DI(qdown(-term), qup(term))
            s = s + DI.q(term if k % 2 == 0 else -term)
            k += 1
    return arctan_inv(5) * 16 - arctan_inv(239) * 4

def selftests(c):
    ok = True
    pi = pi_machin()
    z2 = zeta2(2) + 1; z4 = zeta2(4) + 1
    ref2 = pi * pi / DI.q(6); ref4 = pi * pi * pi * pi / DI.q(90)
    tol = D(10) ** -60
    t1 = abs(z2.lo - ref2.lo) < tol and abs(z4.lo - ref4.lo) < tol
    with localcontext() as ctx:
        ctx.prec = DPREC
        ln2 = D(2).ln()
    t2 = abs(eta(1).lo - ln2) < tol
    print("self-test zeta(2) = pi^2/6, zeta(4) = pi^4/90 (pi by Machin): %s; eta(1) = ln 2: %s" % (t1, t2))
    ok &= t1 and t2
    # (T phi)' = -V psi with phi = int_0^x psi, against the Gauss-Kuzmin operator summed directly (floating point)
    from math import fsum
    def phi(t):
        r = 0.0
        for j in range(len(c) - 1, -1, -1):
            r = r * t + float(c[j]) / (j + 1)
        return r * t
    def Tphi(x, K=400000):
        return fsum((1 + x) / ((x + k) * (x + k + 1)) * phi(1 / (x + k)) for k in range(1, K + 1))
    x, h = 0.3, 1e-4
    deriv = (Tphi(x + h) - Tphi(x - h)) / (2 * h)
    Vx = evaluate_V(c, Q(3, 10))
    t3 = abs(deriv + float(Vx.lo)) < 1e-7
    print("self-test (T phi)'(0.3) = -V psi(0.3), T summed directly: %s (%.10f vs %.10f)" % (t3, -deriv, float(Vx.lo)))
    ok &= t3
    return ok

_cache = {}
def build(c, variant, S=150, Mt=200, rho=Q(6, 5)):
    """V psi(1/2 + t) = sum_m g_m t^m + E, |E| <= rem, for |t| <= 1/2"""
    key = (tuple(c), variant)
    if key in _cache:
        return _cache[key]
    ser = summand_series(c, S, variant)
    for s in range(3):
        if ser[s] != [0, 0]:
            raise SystemExit("the summand is not O(w^3): coefficient %d = %s" % (s, ser[s]))
    M = cauchy_M(c, rho)
    q = Q(2, 3) / rho
    tail = M * 2 * q ** (S + 1) / (1 - q)       # s > S: |coef| (2/3)^s (1 + 3/(2(s-1))) <= M rho^-s (2/3)^s * 2
    Z = {n: zeta2(n) for n in range(3, S + Mt + 1)}
    g = [DI(D(0)) for _ in range(Mt + 2)]
    rem = tail
    for s in range(3, S + 1):
        a0, a1 = ser[s]
        if a0 == 0 and a1 == 0: continue
        lin0 = a0 + a1 / 2                       # alpha + beta x at x = 1/2 ; slope a1 in t
        for m in range(Mt + 1):
            zt = Z[s + m] * ((-1) ** m * comb(s + m - 1, m))
            g[m] = g[m] + zt * lin0
            if a1:
                g[m + 1] = g[m + 1] + zt * a1
        n = s + Mt + 1
        zb = Q(2, 3) ** n * (1 + Q(3, 2 * (n - 1)))
        rem += (abs(lin0) + abs(a1) / 2) * comb(s + Mt, Mt + 1) * zb / 2 ** (Mt + 1)
    _cache[key] = (g, rem, M)
    return g, rem, M

def evaluate_V(c, x, variant=None):
    g, rem, _ = build(c, variant)
    t = DI.q(Q(x) - Q(1, 2))
    acc = DI(D(0))
    for gm in reversed(g):
        acc = acc * t + gm
    return acc + DI(qdown(-rem), qup(rem))

def shift(coefs, t0):
    """coefficients of p(t0 + u) in u, intervals"""
    out = [DI(D(0)) for _ in coefs]
    for cj in reversed(coefs):                  # out = out*(t0+u) + cj
        new = [DI(D(0)) for _ in coefs]
        for i in range(len(coefs) - 1):
            if i + 1 < len(new):
                new[i + 1] = new[i + 1] + out[i]
        for i in range(len(coefs)):
            new[i] = new[i] + out[i] * t0
        new[0] = new[0] + cj
        out = new
    return out

def cell_lower(h, r, rem):
    with _ctx(ROUND_CEILING):
        s = D(0)
        rp = D(1)
        for j in range(1, len(h)):
            rp = rp * r
            s = s + h[j].mag() * rp
        s = s + rem
    with _ctx(ROUND_FLOOR):
        return h[0].lo - s

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    a = sys.argv[1:]
    opt = {"--psi": os.path.join(here, "psi.json"), "--cells": "32"}
    i = 0
    while i < len(a):
        opt[a[i]] = a[i + 1]; i += 2
    data = json.load(open(opt["--psi"]))
    c = [Q(v) for v in data["psi_monomial_coefficients"]]
    lo = Q(opt.get("--claim-lower", data["claimed_lower"]))
    hi = Q(opt.get("--claim-upper", data["claimed_upper"]))
    if "--perturb" in opt:
        j, dd = opt["--perturb"].split(":")
        c[int(j)] += Q(dd)
        print("CONTROL: coefficient %s of psi changed by %s" % (j, dd))
    variant = opt.get("--variant")
    if variant:
        print("CONTROL: operator variant", variant, "(the point term with (x+k)^2 in place of (x+k)^3)")
    print("83a program B: psi of degree %d, claimed %s <= |lambda_2| <= %s" % (len(c) - 1, qdown(lo), qup(hi)))
    if not selftests(c):
        print("VERDICT FAILED (self-test)"); return
    g, rem, M = build(c, variant)
    print("series in w to w^150, Cauchy bound M = %.3e on |w| = 6/5; global Taylor polynomial in t of degree %d; "
          "uniform remainder %.3e" % (float(M), len(g) - 1, float(rem)))
    # psi(1/2 + t) exactly
    n = len(c)
    pt = [Q(0)] * n
    for cj in reversed(c):
        new = [Q(0)] * n
        for i2 in range(n - 1):
            new[i2 + 1] += pt[i2]
        for i2 in range(n):
            new[i2] += pt[i2] * Q(1, 2)
        new[0] += cj
        pt = new
    L = len(g)
    psit = [DI.q(v) for v in pt] + [DI(D(0)) for _ in range(L - n)]
    glo = [g[m] - psit[m] * lo for m in range(L)]
    ghi = [psit[m] * hi - g[m] for m in range(L)]
    NC = int(opt["--cells"])
    r = qup(Q(1, 2 * NC))
    res = {}
    for name, poly, rr in (("psi", psit, Q(0)), ("lo", glo, rem), ("hi", ghi, rem)):
        worst = None; bad = []
        for k in range(NC):
            t0 = Q(2 * k + 1, 2 * NC) - Q(1, 2)
            h = shift(poly, DI.q(t0))
            Lb = cell_lower(h, r, qup(rr))
            if worst is None or Lb < worst[0]:
                worst = (Lb, t0 + Q(1, 2))
            if Lb <= 0:
                bad.append((t0 + Q(1, 2), h[0]))
        res[name] = (worst, bad)
    def show(label, name):
        (w, x0), bad = res[name]
        print("%s: %s (least cell bound %.6e on the cell at x = %s)" % (label, not bad, float(w), x0))
        for x0b, v in bad[:3]:
            if v.hi < 0:
                print("  the inequality FAILS at x = %s: value in [%.6e, %.6e]" % (x0b, float(v.lo), float(v.hi)))
            else:
                print("  not proved on the cell at x = %s" % x0b)
    print("cells: %d of width 1/%d" % (NC, NC))
    show("psi > 0 on [0,1]", "psi")
    show("V psi - %s psi >= 0 on [0,1]" % qdown(lo), "lo")
    show("%s psi - V psi >= 0 on [0,1]" % qup(hi), "hi")
    if res["psi"][1] or res["lo"][1] or res["hi"][1]:
        print("VERDICT REFUSED: the claimed bounds are not proved")
        return
    print("PROVED: %s <= |lambda_2| <= %s" % (qdown(lo), qup(hi)))
    # 2/(3 + sqrt 5) > hi  <=>  (2/hi - 3)^2 > 5 with 2/hi - 3 > 0
    s5 = 2 / hi - 3
    row_s = s5 > 0 and s5 * s5 > 5
    row_up = hi <= Q("0.30366327")
    row_lo = lo >= Q("0.30366299")
    print("row [Schweiger1980] C_83 <= 2/(3 + sqrt 5) = 0.381966...:", row_s)
    print("row [MR1987] C_83 <= 0.30366327:", row_up)
    print("row [MR1987] C_83 >= 0.30366299:", row_lo)
    centre = Q("0.30366300289873265859744812190155623311087735225365")
    print("information: the 175-digit centre of [GKW2026] (first 50 digits) lies in [lo, hi]:", lo <= centre <= hi)
    print("VERDICT PASS" if (row_s and row_up and row_lo) else "VERDICT REFUSED: a row does not follow")

if __name__ == "__main__":
    main()
