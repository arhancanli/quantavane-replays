#!/usr/bin/env python3
"""80a, program A: the Krauth-Mezard threshold alpha_star of the Ising perceptron (kappa = 0), from its definition.

Rows: [H2024] (upper) and [DS2025] (lower), both alpha_star ~ 0.833, both conditional on numerical hypotheses that
are NOT checked here. A value from its definition, not a certificate.

Definition (Huang, arXiv 2404.18902, Sec. 3; Ding and Sun, arXiv 1809.07742, eqs. (q.psi.recursion), (GG) and
Proposition 1.3): with phi the standard normal density, Phi its distribution function and Z ~ N(0, 1),
    P(psi)      = E th(psi^(1/2) Z)^2,
    R_alpha(q)  = alpha E F_q(q^(1/2) Z)^2,  F_q(x) = phi(x/s)/(s Phi(x/s)), s = (1-q)^(1/2)
                = alpha/(1-q) E (phi(cZ)/Phi(cZ))^2,  c = (q/(1-q))^(1/2),
    G(alpha, q, psi) = -(1-q) psi/2 + E log(2 ch(psi^(1/2) Z)) + alpha E log Phi(cZ).
q_star(alpha) is the fixed point of q = P(R_alpha(q)) in (q_lb, q_ub) = (0.56394907949, 0.56394908030),
psi_star = R_alpha(q_star), G_star(alpha) = G(alpha, q_star, psi_star), and alpha_star is the root of G_star on
(0.833078599, 0.833078600); Ding and Sun's Proposition 1.3 states that G_star is well defined and strictly decreasing
there (this program does not re-check that; it checks the values).

This program (route A): for alpha = a_lo and a_hi, a floating-point iteration gives q~; then, rigorously,
P(R(q~ - d)) > q~ - d and P(R(q~ + d)) < q~ + d (d = 1e-13), so q_star is in Q = [q~ - d, q~ + d]; then
G(alpha, Q, R_alpha(Q)) is enclosed with Q as an interval parameter. Claimed: G_star(a_lo) > 0 > G_star(a_hi).
Integrals: Taylor models. On cells of radius r the integrand's Taylor polynomial at the cell centre (interval
Taylor-mode arithmetic, order p) is integrated exactly, plus the Lagrange remainder from the (p+1)-th coefficient
enclosed over the whole cell; tails |z| > 10 by explicit bounds (Gordon's inequality for the Mills ratio).
Arithmetic: decimal intervals, 60 digits; exp, ln and sqrt are correctly rounded in the decimal module and widened by
one unit in the last place; Phi by its everywhere-convergent series at 110 digits. Standard library only.
Also [AT2024] 0.847 (Altschuler and Tikhomirov, arXiv 2401.15092, Proposition and its proof): the capacity is at
most alpha when the Gardner-Derrida value GD(alpha) = min_q GD(alpha, q) is below -log 2, and
GD(alpha, 1/2) = alpha E log H(u) + (1 + log(1/2))/2 with H(s) = P(Z >= s), E log H(u) = -1 (H(u) is uniform).
This program computes E log H(u) = E log Phi(Z) by the same quadrature (c = 1) and GD(0.847, 1/2) in intervals.
Usage: c80a_a.py [--alpha-lo A] [--alpha-hi A] [--variant no-prefactor|sign-reversed] [--at-alpha A]
"""
import sys
from fractions import Fraction as Q
from decimal import Decimal as D, localcontext, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from math import tanh as ftanh, cosh as fcosh, log as flog, sqrt as fsqrt, exp as fexp, erfc as ferfc, pi as fpi

DPREC = 60
def _ctx(rnd, prec=DPREC):
    return localcontext(Context(prec=prec, rounding=rnd))

class DI:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        self.lo = lo; self.hi = lo if hi is None else hi
    @staticmethod
    def q(v):
        v = Q(v)
        with _ctx(ROUND_FLOOR): lo = D(v.numerator) / D(v.denominator)
        with _ctx(ROUND_CEILING): hi = D(v.numerator) / D(v.denominator)
        return DI(lo, hi)
    def __add__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        with _ctx(ROUND_FLOOR): lo = self.lo + o.lo
        with _ctx(ROUND_CEILING): hi = self.hi + o.hi
        return DI(lo, hi)
    __radd__ = __add__
    def __neg__(self):
        return DI(self.hi.copy_negate(), self.lo.copy_negate())
    def __sub__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        return self + (-o)
    def __rsub__(self, o):
        return DI.q(o) - self
    def __mul__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        a, b, c, d = self.lo, self.hi, o.lo, o.hi
        with _ctx(ROUND_FLOOR): lo = min(a * c, a * d, b * c, b * d)
        with _ctx(ROUND_CEILING): hi = max(a * c, a * d, b * c, b * d)
        return DI(lo, hi)
    __rmul__ = __mul__
    def recip(self):
        if not (self.lo > 0 or self.hi < 0):
            raise ArithmeticError("division by an interval containing 0")
        with _ctx(ROUND_FLOOR): lo = D(1) / self.hi
        with _ctx(ROUND_CEILING): hi = D(1) / self.lo
        return DI(lo, hi)
    def __truediv__(self, o):
        if not isinstance(o, DI): o = DI.q(o)
        return self * o.recip()
    def mag(self):
        return max(self.lo.copy_abs(), self.hi.copy_abs())
    def hull(self, o):
        return DI(min(self.lo, o.lo), max(self.hi, o.hi))
    def mid(self):
        return (self.lo + self.hi) / 2

def _cr(f, x):
    """a correctly rounded (half-even) decimal function, widened by one ulp each way"""
    with _ctx(ROUND_HALF_EVEN) as c:
        v = f(x, c)
        return v.next_minus(c), v.next_plus(c)

def exp_di(x):
    lo, _ = _cr(lambda v, c: v.exp(c), x.lo)
    _, hi = _cr(lambda v, c: v.exp(c), x.hi)
    return DI(lo, hi)
def ln_di(x):
    assert x.lo > 0
    lo, _ = _cr(lambda v, c: v.ln(c), x.lo)
    _, hi = _cr(lambda v, c: v.ln(c), x.hi)
    return DI(lo, hi)
def sqrt_di(x):
    assert x.lo >= 0
    lo, _ = _cr(lambda v, c: v.sqrt(c), x.lo)
    _, hi = _cr(lambda v, c: v.sqrt(c), x.hi)
    return DI(max(lo, D(0)), hi)

def machin_pi_exact(digits):
    """a rational within 10^-digits of pi (Machin, alternating series cut where the term is below 10^-(digits+5))"""
    def at(m):
        s = Q(0); k = 0
        while True:
            t = Q(1, (2 * k + 1) * m ** (2 * k + 1))
            if t < Q(1, 10 ** (digits + 5)):
                return s
            s += t if k % 2 == 0 else -t; k += 1
    return 16 * at(5) - 4 * at(239)
_pi = machin_pi_exact(DPREC + 5)
PI = DI.q(_pi - Q(1, 10 ** (DPREC + 3))).hull(DI.q(_pi + Q(1, 10 ** (DPREC + 3))))
_pi_hp = machin_pi_exact(130)
INV_SQRT_2PI = sqrt_di(PI * 2).recip()

def Phi_point(x):
    """Phi(x) for a decimal x: 1/2 + phi(x) sum_n x^(2n+1)/(2n+1)!!, at 110 digits; tail <= 2 |last term|"""
    P = 110
    with _ctx(ROUND_HALF_EVEN, P) as c:
        x2 = x * x
        term = x; s = x; n = 0
        while True:
            n += 1
            term = term * x2 / (2 * n + 1)
            s += term
            if 2 * n + 3 > 2 * x2 and abs(term) < D(10) ** -(P - 5):
                break
        pi = D(_pi_hp.numerator) / D(_pi_hp.denominator)       # within 10^-125 of pi
        ph = (-x2 / 2).exp() / (2 * pi).sqrt()
        v = D(1) / 2 + ph * s
        err = ph * abs(term) * 2 + abs(v) * D(10) ** -(P - 10) + D(10) ** -(P - 10)
    with _ctx(ROUND_FLOOR): lo = v - err
    with _ctx(ROUND_CEILING): hi = v + err
    return DI(lo, hi)

def Phi_di(x):
    return DI(Phi_point(x.lo).lo, Phi_point(x.hi).hi)

# ---------------- truncated Taylor series with interval coefficients ----------------
class TS:
    __slots__ = ("c",)
    def __init__(self, c):
        self.c = c
    @staticmethod
    def var(x0, n):
        return TS([x0, DI(D(1))] + [DI(D(0)) for _ in range(n - 2)])
    @staticmethod
    def const(v, n):
        v = v if isinstance(v, DI) else DI.q(v)
        return TS([v] + [DI(D(0)) for _ in range(n - 1)])
    def __len__(self):
        return len(self.c)
    def __add__(self, o):
        if isinstance(o, TS):
            return TS([a + b for a, b in zip(self.c, o.c)])
        return TS([self.c[0] + o] + self.c[1:])
    __radd__ = __add__
    def __neg__(self):
        return TS([-a for a in self.c])
    def __sub__(self, o):
        return self + (-o if isinstance(o, TS) else -(o if isinstance(o, DI) else DI.q(o)))
    def __mul__(self, o):
        if not isinstance(o, TS):
            return TS([a * o for a in self.c])
        n = len(self.c)
        return TS([sum((self.c[j] * o.c[k - j] for j in range(k + 1)), DI(D(0))) for k in range(n)])
    __rmul__ = __mul__
    def recip(self):
        n = len(self.c); f = self.c; r0 = f[0].recip(); r = [r0]
        for k in range(1, n):
            s = sum((f[j] * r[k - j] for j in range(1, k + 1)), DI(D(0)))
            r.append(-(s * r0))
        return TS(r)
    def __truediv__(self, o):
        return self * o.recip()

def sq(a):
    """series square"""
    n = len(a)
    return [sum((a[j] * a[k - j] for j in range(k + 1)), DI(D(0))) for k in range(n)]

def phi_point(x):
    """phi at a decimal interval x"""
    x2 = x * x
    return exp_di(DI(x2.hi.copy_negate(), x2.lo.copy_negate()) * DI.q(Q(1, 2))) * INV_SQRT_2PI if x.lo >= 0 or x.hi <= 0 \
        else exp_di(x2 * DI.q(Q(-1, 2))) * INV_SQRT_2PI

def phi_series(z0, n):
    """phi(z0 + u): (k+1) phi_(k+1) = -(z0 phi_k + phi_(k-1))"""
    f = [phi_point(z0)]
    for k in range(n - 1):
        s = z0 * f[k] + (f[k - 1] if k >= 1 else DI(D(0)))
        f.append(-s / (k + 1))
    return f

def th_point(x):                          # th = 1 - 2/(e^(2x) + 1), increasing; x used once
    return DI.q(1) - (exp_di(x * 2) + 1).recip() * 2

def th_series(sp, z0, n):
    """y = th(sp (z0 + u)): (k+1) y_(k+1) = sp ([k = 0] - (y^2)_k)"""
    y = [th_point(sp * z0)]
    for k in range(n - 1):
        y2k = sum((y[j] * y[k - j] for j in range(k + 1)), DI(D(0)))
        y.append(sp * ((DI.q(1) if k == 0 else DI(D(0))) - y2k) / (k + 1))
    return y

def mills_point(t):
    """m(t) = phi(t)/Phi(t) at a decimal interval t; m is decreasing (log-concavity of Phi)"""
    def at(x):
        return phi_point(DI(x)) / Phi_point(x)
    return DI(at(t.hi).lo, at(t.lo).hi)

def mills_series(c, z0, n):
    """m(c (z0 + u)): dm/du = -c m (t + m), t = c z0 + c u"""
    m = [mills_point(c * z0)]
    t0 = c * z0
    for k in range(n - 1):
        mt = t0 * m[k] + (c * m[k - 1] if k >= 1 else DI(D(0)))
        mm = sum((m[j] * m[k - j] for j in range(k + 1)), DI(D(0)))
        m.append(-(c * (mt + mm)) / (k + 1))
    return m

def logPhi_point(t):                     # increasing
    return DI(ln_di(Phi_point(t.lo)).lo, ln_di(Phi_point(t.hi)).hi)

def mul_series(a, b):
    n = len(a)
    return [sum((a[j] * b[k - j] for j in range(k + 1)), DI(D(0))) for k in range(n)]

def integrand(kind, z0, n, par):
    """Taylor coefficients (orders 0..n-1) at z0 (a point or a cell) of the integrand times phi(z)"""
    w = phi_series(z0, n)
    if kind == "P":
        y = th_series(par["sp"], z0, n)
        return mul_series(sq(y), w)
    if kind == "R":
        m = mills_series(par["c"], z0, n)
        return mul_series(sq(m), w)
    if kind == "LC":                      # (log 2ch)' = sp th
        sp = par["sp"]; y = th_series(sp, z0, n)
        x0 = sp * z0
        l0 = ln_di(exp_di(x0) + exp_di(-x0))
        l = [l0] + [sp * y[k] / (k + 1) for k in range(n - 1)]
        return mul_series(l, w)
    if kind == "LP":                      # (log Phi(cz))' = c m
        c = par["c"]; m = mills_series(c, z0, n)
        l = [logPhi_point(c * z0)] + [c * m[k] / (k + 1) for k in range(n - 1)]
        return mul_series(l, w)
    raise ValueError(kind)

L = 10
ORDER = 14
NCELL_PER_UNIT = 16                       # cell radius 1/32

def integrate(kind, par, even):
    """int over [-L, L] (or 2 int over [0, L] if even) of the integrand, as an interval"""
    n = ORDER + 1
    lo_end = 0 if even else -L
    ncell = (L - lo_end) * NCELL_PER_UNIT
    r = Q(1, 2 * NCELL_PER_UNIT)
    rD = DI.q(r)
    tot = DI(D(0))
    for i in range(ncell):
        m = Q(lo_end) + r * (2 * i + 1)
        mD = DI.q(m)
        g = integrand(kind, mD, n, par)
        for j in range(0, n, 2):
            tot = tot + g[j] * DI.q(2 * r ** (j + 1) / (j + 1))
        cellI = DI.q(m - r).hull(DI.q(m + r))
        gi = integrand(kind, cellI, n + 1, par)
        b = gi[n].mag()
        with _ctx(ROUND_CEILING):
            e = b * D(2) * (rD.hi ** (n + 1)) / (n + 1)
        tot = tot + DI(e.copy_negate(), e)
    return tot * 2 if even else tot

def tails(kind, par):
    """bounds for the integral over |z| > L"""
    Lq = Q(L)
    phL = phi_point(DI.q(Lq))
    PsiL = phL / L                          # Psi(L) <= phi(L)/L
    if kind == "P":
        return DI(D(0), (PsiL * 2).hi)
    if kind == "LC":                        # log 2 <= log 2ch y <= log 2 + |y|
        return DI(D(0), (PsiL * 2 * DI.q(Q(7, 10)) + phL * 2 * par["sp"]).hi)
    c = par["c"]
    if kind == "R":                         # z < -L: (phi/Phi)(t) <= |t| + 1 (Gordon); z > L: <= phi(0)/Phi(0) < 1
        up = c * c * (phL * L + PsiL) + c * phL * 2 + PsiL + PsiL
        return DI(D(0), up.hi)
    if kind == "LP":                        # z < -L: |log Phi(t)| <= t^2/2 + |t| + 2 ; z > L: |log Phi| <= 1
        lo = c * c * (phL * L + PsiL) / 2 + c * phL + PsiL * 2 + PsiL
        return DI(lo.hi.copy_negate(), D(0))
    raise ValueError(kind)

def E(kind, par, even=False):
    return integrate(kind, par, even) + tails(kind, par)

def P_of(psi):
    return E("P", {"sp": sqrt_di(psi)}, even=True)

def c_of(q):
    return sqrt_di(q / (DI.q(1) - q))

def R_of(q, alpha, variant=None):
    v = E("R", {"c": c_of(q)})
    if variant == "no-prefactor":
        return v * alpha                     # F_q without its 1/(1-q)^(1/2)
    return v * alpha / (DI.q(1) - q)

def G_of(alpha, q, psi, variant=None):
    lc = E("LC", {"sp": sqrt_di(psi)}, even=True)
    lp = E("LP", {"c": c_of(q)})
    sign = -1 if variant == "sign-reversed" else 1
    return -((DI.q(1) - q) * psi) / 2 + lc + alpha * lp * sign

# ---------------- floating-point search for q~ ----------------
def float_fixed_point(alpha, variant=None):
    H = 0.004
    zs = [i * H for i in range(-3000, 3001)]
    ws = [H * fexp(-z * z / 2) / fsqrt(2 * fpi) for z in zs]
    Phi = lambda t: 0.5 * ferfc(-t / fsqrt(2))
    ph = lambda t: fexp(-t * t / 2) / fsqrt(2 * fpi)
    def P(psi):
        s = fsqrt(psi); return sum(w * ftanh(s * z) ** 2 for z, w in zip(zs, ws))
    def R(q):
        c = fsqrt(q / (1 - q)); f = alpha / (1 - q) if variant != "no-prefactor" else alpha
        return f * sum(w * (ph(c * z) / Phi(c * z)) ** 2 for z, w in zip(zs, ws))
    q = 0.56
    for _ in range(300):
        qn = P(R(q))
        if abs(qn - q) < 1e-16: break
        q = qn
    return q

def selftest():
    ok = True
    one = E("P", {"sp": DI.q(0)}, even=True)            # th(0)^2 = 0: integral 0
    p1 = Phi_point(D(1)); p0 = Phi_point(D(0))
    lp0 = integrate("LP", {"c": DI.q(0)}, False)       # log Phi(0) = -log 2, so the integral is -log 2 (1 - 2 Psi(10))
    with _ctx(ROUND_HALF_EVEN, 100):           # the comparisons themselves at 100 digits
        ref1 = D("0.841344746068542948585232545632")       # 30 digits
        pm = Phi_point(D(-1))
        t1 = (p1.lo <= ref1 + D(10) ** -30 and ref1 - D(10) ** -30 <= p1.hi and p0.lo <= D("0.5") <= p0.hi
              and (p1 + pm).lo <= 1 <= (p1 + pm).hi and (p1 + pm).hi - (p1 + pm).lo < D(10) ** -50)
        ln2 = D(2).ln()
        t2 = abs(lp0.lo + ln2) < D(10) ** -15 and abs(lp0.hi + ln2) < D(10) ** -15
        t3 = one.hi < D(10) ** -18
    print("self-test Phi(0) = 1/2, Phi(1) to 30 digits, Phi(1) + Phi(-1) = 1: %s; int log Phi(0) phi = -log 2 (Taylor-model quadrature): %s; "
          "int th(0)^2 phi = 0: %s" % (t1, t2, t3))
    return ok and t1 and t2 and t3

def at2024(alpha):
    """[AT2024]: GD(alpha, 1/2) < -log 2 ?"""
    print("80a program A, [AT2024]: GD(%s, 1/2) = %s E log H(u) + (1 + log(1/2))/2 against -log 2" % (alpha, alpha))
    if not selftest():
        print("VERDICT FAILED (self-test)"); return
    elh = E("LP", {"c": DI.q(1)})                      # E log Phi(Z) = E log H(-Z) = E log H(u)
    print("E log H(u) by quadrature in [%.18f, %.18f] (exactly -1: H(u) is uniform on (0, 1))" % (float(elh.lo), float(elh.hi)))
    ln2 = ln_di(DI.q(2))
    gd = DI.q(alpha) * elh + (DI.q(1) - ln2) / 2
    gap = gd + ln2                                      # GD(alpha, 1/2) + log 2
    print("GD(%s, 1/2) + log 2 in [%.12e, %.12e]" % (alpha, float(gap.lo), float(gap.hi)))
    print("information: the same with E log H(u) = -1 exactly: alpha - 1/2 - (log 2)/2 > 0 iff alpha > 1/2 + (log 2)/2 = %.12f"
          % float(((DI.q(1) + ln2) / 2).lo))
    ok = gap.hi < 0
    print("row [AT2024] capacity <= %s: GD(%s) <= GD(%s, 1/2) < -log 2: %s" % (alpha, alpha, alpha, ok))
    print("VERDICT PASS" if ok else "VERDICT REFUSED: GD(alpha, 1/2) is not below -log 2")

def main():
    a = sys.argv[1:]
    opt = {"--alpha-lo": "0.833078599", "--alpha-hi": "0.833078600"}
    for i in range(0, len(a), 2):
        opt[a[i]] = a[i + 1]
    variant = opt.get("--variant")
    if "--at-alpha" in opt:
        return at2024(Q(opt["--at-alpha"]))
    if variant:
        print("CONTROL: variant", variant)
    print("80a program A: claimed alpha_star in (%s, %s): G_star(lo) > 0 > G_star(hi)" % (opt["--alpha-lo"], opt["--alpha-hi"]))
    if not selftest():
        print("VERDICT FAILED (self-test)"); return
    qlb, qub = Q("0.56394907949"), Q("0.56394908030")
    d = Q(1, 10 ** 13)
    res = {}
    for name in ("--alpha-lo", "--alpha-hi"):
        alpha = DI.q(Q(opt[name]))
        qt = Q(repr(float_fixed_point(float(Q(opt[name])), variant)))
        q1, q2 = qt - d, qt + d
        f1 = P_of(R_of(DI.q(q1), alpha, variant)) - q1
        f2 = P_of(R_of(DI.q(q2), alpha, variant)) - q2
        bracket = f1.lo > 0 and f2.hi < 0
        inside = qlb < q1 and q2 < qub
        Qi = DI.q(q1).hull(DI.q(q2))
        psi = R_of(Qi, alpha, variant)
        G = G_of(alpha, Qi, psi, variant)
        print("alpha = %s: q~ = %.16f; P(R(q~ - d)) - (q~ - d) in [%.3e, %.3e], P(R(q~ + d)) - (q~ + d) in [%.3e, %.3e] "
              "(d = 1e-13): fixed point bracketed %s, inside (q_lb, q_ub) %s"
              % (opt[name], float(qt), float(f1.lo), float(f1.hi), float(f2.lo), float(f2.hi), bracket, inside))
        print("   q_star in [%s, %s]; psi_star in [%s, %s]" % (str(+Qi.lo)[:20], str(+Qi.hi)[:20], str(psi.lo)[:20], str(psi.hi)[:20]))
        print("   G_star in [%.6e, %.6e]" % (float(G.lo), float(G.hi)))
        res[name] = (bracket and inside, G)
    ok_lo = res["--alpha-lo"][0] and res["--alpha-lo"][1].lo > 0
    ok_hi = res["--alpha-hi"][0] and res["--alpha-hi"][1].hi < 0
    print("G_star(%s) > 0: %s; G_star(%s) < 0: %s" % (opt["--alpha-lo"], ok_lo, opt["--alpha-hi"], ok_hi))
    if ok_lo and ok_hi:
        print("RESULT: alpha_star in (%s, %s) (with Proposition 1.3's monotonicity, not re-checked)" % (opt["--alpha-lo"], opt["--alpha-hi"]))
        lo, hi = Q(opt["--alpha-lo"]), Q(opt["--alpha-hi"])
        rounds = Q("0.8325") <= lo and hi < Q("0.8335")
        print("rows [H2024] and [DS2025]: alpha_star ~ 0.833 (rounds to 0.833): %s" % rounds)
        print("VERDICT PASS" if rounds else "VERDICT REFUSED: the interval does not round to 0.833")
    else:
        print("VERDICT REFUSED: the claimed interval is not established")

if __name__ == "__main__":
    main()
