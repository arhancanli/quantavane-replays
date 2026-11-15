#!/usr/bin/env python3
"""80a, program B (written separately from program A): the Krauth-Mezard threshold alpha_star (kappa = 0).

Same definition and claim as program A (see there): G_star(a_lo) > 0 > G_star(a_hi), with
    P(psi) = E th(psi^(1/2) Z)^2,  R_alpha(q) = alpha/(1-q) E m(cZ)^2,  m = phi/Phi,  c = (q/(1-q))^(1/2),
    G(alpha, q, psi) = -(1-q) psi/2 + E log(2 ch(psi^(1/2) Z)) + alpha E log Phi(cZ).
A value from its definition, not a certificate; Ding and Sun's Proposition 1.3 (monotonicity, uniqueness) not
re-checked.

Route B:
  * The fixed point in psi: psi = R_alpha(P(psi)); a floating-point iteration gives psi~, then rigorously
    R(P(psi~ - e)) > psi~ - e and R(P(psi~ + e)) < psi~ + e, so psi_star is in S = [psi~ - e, psi~ + e] and
    q_star = P(psi_star) is in P(S); G is enclosed with S and P(S) as interval parameters.
  * E th^2 and E log 2ch in z on [0, 9.5] (even integrands); E m(cZ)^2 and E log Phi(cZ) in t = cz on [-11, 11]
    with the weight phi(t/c)/c. Each cell [a, a + h] (h = 1/32): the Taylor polynomial at the LEFT end a (order 13)
    integrated over [0, h], plus a Lagrange remainder from the 14th coefficient enclosed over the cell.
    Series by their differential equations: th' = 1 - th^2, (log 2ch)' = th, m' = -m (t + m), (log Phi)' = m,
    (phi(t/c))' = -(t/c^2) phi(t/c).
  * Phi: the erf series (alternating) for |t| <= 4, Laplace's continued fraction Psi(y) = phi(y)/(y + 1/(y + 2/(y + ...)))
    for y > 4 (successive approximants bracket the value); then m(t) = 1/K(-t) and log Phi(t) = log phi(t) + log K(-t)
    on the far left, with no subtraction.
Arithmetic: integer fixed point at 2^-200, interval ends rounded outward; exp (halving and squaring), log (atanh
series), sqrt (isqrt) and pi (Machin) are this program's own. Standard library only.
Also [AT2024] 0.847: GD(alpha, 1/2) = alpha E log H(u) + 1/2 + log(1/2)/2 < -log 2 (see program A); here
E log H(u) = E log Phi(t) under phi(t), by this program's quadrature in t with c = 1.
Usage: c80a_b.py [--alpha-lo A] [--alpha-hi A] [--variant no-prefactor|sign-reversed] [--at-alpha A]
"""
import sys
from fractions import Fraction as Q
from math import isqrt, tanh as ftanh, exp as fexp, sqrt as fsqrt, erfc as ferfc, pi as fpi

F = 200
ONE = 1 << F

def fl(a):
    return a >> F                         # floor of a / 2^F
def ce(a):
    return -((-a) >> F)

class FI:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        self.lo = lo; self.hi = lo if hi is None else hi
    @staticmethod
    def q(v):
        v = Q(v)
        return FI((v.numerator << F) // v.denominator, -((-(v.numerator << F)) // v.denominator))
    def __add__(self, o):
        if not isinstance(o, FI): o = FI.q(o)
        return FI(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__
    def __neg__(self):
        return FI(-self.hi, -self.lo)
    def __sub__(self, o):
        if not isinstance(o, FI): o = FI.q(o)
        return FI(self.lo - o.hi, self.hi - o.lo)
    def __rsub__(self, o):
        return FI.q(o) - self
    def __mul__(self, o):
        if isinstance(o, int):
            return FI(self.lo * o, self.hi * o) if o >= 0 else FI(self.hi * o, self.lo * o)
        if not isinstance(o, FI): o = FI.q(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return FI(fl(min(p)), ce(max(p)))
    __rmul__ = __mul__
    def recip(self):
        if self.lo > 0:
            return FI((ONE << F) // self.hi, -((-(ONE << F)) // self.lo))
        if self.hi < 0:
            return -((-self).recip())
        raise ArithmeticError("1/x with 0 in x")
    def __truediv__(self, o):
        if isinstance(o, int) and o > 0:
            return FI(self.lo // o, -((-self.hi) // o))
        if not isinstance(o, FI): o = FI.q(o)
        return self * o.recip()
    def hull(self, o):
        return FI(min(self.lo, o.lo), max(self.hi, o.hi))
    def mag(self):
        return max(abs(self.lo), abs(self.hi))
    def f(self):
        return self.lo / ONE

ZERO = FI(0)

# ---------------- elementary functions (this program's own) ----------------
def _exp_point(a):
    """exp(a / 2^F) as an interval"""
    s = max(0, abs(a).bit_length() - F + 12)
    y = FI(a >> s, -((-a) >> s))           # |y| < 2^-11
    tot = FI(ONE); term = FI(ONE)
    for k in range(1, 26):
        term = term * y / k
        tot = tot + term
    r = 2 * (y.mag() ** 26 // (ONE ** 25)) // 1 + 1   # 2|y|^26/26! < 2|y|^26 (scaled, rounded up)
    tot = tot + FI(-r, r)
    for _ in range(s):
        tot = tot * tot
    return tot
def exp_fi(x):
    return FI(_exp_point(x.lo).lo, _exp_point(x.hi).hi)

def _atanh_series(u):
    """2 atanh(u) for an interval 0 <= u <= 1/3"""
    tot = ZERO; u2 = u * u; p = u; k = 0
    while True:
        tot = tot + p / (2 * k + 1)
        p = p * u2; k += 1
        if p.hi * 9 < 1 << 10:             # remaining tail <= p/(1 - u^2) < 2 p
            tot = tot + FI(0, 2 * p.hi + 2)
            return tot * 2
LOG2 = _atanh_series(FI.q(Q(1, 3)))
def _log_point(a):
    assert a > 0
    e = a.bit_length() - 1 - F            # a/2^F = 2^e y, 1 <= y < 2
    yl = a << (-e) if e < 0 else a >> e
    yh = yl if e <= 0 else -((-a) >> e)
    y = FI(yl, yh)
    u = (y - ONE * 0 - FI(ONE)) / (y + FI(ONE))
    u = FI(max(u.lo, 0), u.hi)
    return LOG2 * e + _atanh_series(u)
def log_fi(x):
    return FI(_log_point(x.lo).lo, _log_point(x.hi).hi)

def sqrt_fi(x):
    lo = isqrt(max(x.lo, 0) << F)
    hi = isqrt(x.hi << F)
    if hi * hi < x.hi << F: hi += 1
    return FI(lo, hi)

def _arctan_inv(m):
    tot = ZERO; k = 0; mm = m * m; den = m
    while True:
        t = FI.q(Q(1, (2 * k + 1) * den))
        if t.hi < 4:
            return tot + FI(-4, 4)
        tot = tot + (t if k % 2 == 0 else -t)
        den *= mm; k += 1
PI = _arctan_inv(5) * 16 - _arctan_inv(239) * 4
SQRT2 = sqrt_fi(FI(2 * ONE))
INV_SQRT_2PI = sqrt_fi(PI * 2).recip()
LOG_SQRT_2PI = log_fi(PI * 2) / 2

def phi_fi(x):
    x2 = x * x
    if x.lo < 0 < x.hi:
        x2 = FI(0, x2.hi)
    return exp_fi(-x2 / 2) * INV_SQRT_2PI

# ---------------- the normal tail ----------------
def _K(y, N=None):
    """Laplace's continued fraction K(y) = Psi(y)/phi(y) = 1/(y + 1/(y + 2/(y + ...))), y > 4 a point; brackets"""
    def approx(N):
        v = FI(y)
        for k in range(N, 0, -1):
            v = FI(y) + FI.q(k) * v.recip()
        return v.recip()
    N = N or 120
    while True:
        a, b = approx(N), approx(N + 1)
        lo, hi = min(a.lo, b.lo), max(a.hi, b.hi)
        if hi - lo < 1 << (F - 170) or N > 4000:
            return FI(lo, hi)
        N *= 2

def _Psi_series(y):
    """Psi(y) = (1 - erf(y/sqrt 2))/2 for a point 0 <= y <= 4: the alternating erf series"""
    x = FI(y) / SQRT2
    x2 = x * x
    term = x; tot = ZERO; n = 0
    while True:
        tot = tot + term / (2 * n + 1)
        nxt = -(term * x2) / (n + 1)
        n += 1
        if n > x2.hi // ONE + 2 and nxt.mag() < 1 << 8:
            tot = tot + FI(-nxt.mag(), nxt.mag())          # alternating, decreasing: |tail| <= first omitted term
            break
        term = nxt
    erf = tot * 2 * sqrt_fi(PI).recip()
    return (FI(ONE) - erf) / 2

def Psi_point(y):
    """Psi(y) = P(Z > y) for a point y >= 0 (scaled integer)"""
    if y <= 4 * ONE:
        return _Psi_series(y)
    return phi_fi(FI(y)) * _K(y)

def mills_point(t):
    """m(t) = phi(t)/Phi(t) at a point t"""
    if t < -4 * ONE:
        return _K(-t).recip()
    Ph = Psi_point(-t) if t < 0 else FI(ONE) - Psi_point(t)
    return phi_fi(FI(t)) / Ph

def logPhi_point(t):
    if t < -4 * ONE:
        return -(FI(t) * FI(t)) / 2 - LOG_SQRT_2PI + log_fi(_K(-t))
    Ph = Psi_point(-t) if t < 0 else FI(ONE) - Psi_point(t)
    return log_fi(Ph)

# ---------------- series (lists of FI) ----------------
def conv(a, b, k):
    s = ZERO
    for j in range(k + 1):
        s = s + a[j] * b[k - j]
    return s

def th_ser(sp, z0, n):
    """th(sp (z0 + u))"""
    x = sp * z0
    th0 = FI(ONE) - (exp_fi(x * 2) + FI(ONE)).recip() * 2
    y = [th0]
    for k in range(n - 1):
        v = (FI(ONE) if k == 0 else ZERO) - conv(y, y, k)
        y.append(sp * v / (k + 1))
    return y

def phi_ser_z(z0, n):
    """phi(z0 + u): (k+1) f_(k+1) = -(z0 f_k + f_(k-1))"""
    f = [phi_fi(z0)]
    for k in range(n - 1):
        v = z0 * f[k] + (f[k - 1] if k else ZERO)
        f.append(-v / (k + 1))
    return f

def phi_ser_t(t0, c, n):
    """phi((t0 + u)/c)/c: (k+1) f_(k+1) = -(t0 f_k + f_(k-1))/c^2"""
    ic = c.recip(); ic2 = ic * ic
    f = [phi_fi(t0 * ic) * ic]
    for k in range(n - 1):
        v = t0 * f[k] + (f[k - 1] if k else ZERO)
        f.append(-(v * ic2) / (k + 1))
    return f

def mills_interval(t0):
    if t0.lo == t0.hi:
        return mills_point(t0.lo)
    return FI(mills_point(t0.hi).lo, mills_point(t0.lo).hi)          # m decreasing

def mills_ser(t0, n):
    """m(t0 + u): (k+1) m_(k+1) = -((t m)_k + (m^2)_k), t = t0 + u"""
    m = [mills_interval(t0)]
    for k in range(n - 1):
        tm = t0 * m[k] + (m[k - 1] if k else ZERO)
        m.append(-(tm + conv(m, m, k)) / (k + 1))
    return m

def logPhi_interval(t0):
    if t0.lo == t0.hi:
        return logPhi_point(t0.lo)
    return FI(logPhi_point(t0.lo).lo, logPhi_point(t0.hi).hi)

def integrand(kind, x0, n, par):
    if kind == "P":
        y = th_ser(par["sp"], x0, n); w = phi_ser_z(x0, n)
        y2 = [conv(y, y, k) for k in range(n)]
        return [conv(y2, w, k) for k in range(n)]
    if kind == "LC":
        sp = par["sp"]; y = th_ser(sp, x0, n); w = phi_ser_z(x0, n)
        x = sp * x0
        l = [log_fi(exp_fi(x) + exp_fi(-x))] + [sp * y[k] / (k + 1) for k in range(n - 1)]
        return [conv(l, w, k) for k in range(n)]
    if kind == "R":                       # in t: m(t)^2 phi(t/c)/c
        m = mills_ser(x0, n); w = phi_ser_t(x0, par["c"], n)
        m2 = [conv(m, m, k) for k in range(n)]
        return [conv(m2, w, k) for k in range(n)]
    if kind == "LP":                      # in t: log Phi(t) phi(t/c)/c
        m = mills_ser(x0, n); w = phi_ser_t(x0, par["c"], n)
        l = [logPhi_interval(x0)] + [m[k] / (k + 1) for k in range(n - 1)]
        return [conv(l, w, k) for k in range(n)]
    raise ValueError(kind)

ORDER = 13
H = Q(1, 32)
ZL = Q(19, 2)            # z in [0, 9.5] for the even integrands
TL = 11                  # t in [-11, 11] for the others

def integrate(kind, par):
    n = ORDER + 1
    if kind in ("P", "LC"):
        a0, a1 = Q(0), ZL
    else:
        a0, a1 = Q(-TL), Q(TL)
    ncell = int((a1 - a0) / H)
    tot = ZERO
    hp = [FI.q(H ** (j + 1) / (j + 1)) for j in range(n + 1)]
    for i in range(ncell):
        a = a0 + H * i
        g = integrand(kind, FI.q(a), n, par)
        for j in range(n):
            tot = tot + g[j] * hp[j]
        gi = integrand(kind, FI.q(a).hull(FI.q(a + H)), n + 1, par)
        e = (gi[n].mag() * hp[n].hi >> F) + 1
        tot = tot + FI(-e, e)
    return tot * 2 if kind in ("P", "LC") else tot

def tails(kind, par):
    """|z| > 9.5 (P, LC) or |t| > 11, i.e. |z| > 11/c (R, LP)"""
    if kind in ("P", "LC"):
        Lz = FI.q(ZL)
    else:
        Lz = FI.q(TL) / par["c"]
        Lz = FI(Lz.lo, Lz.lo)                # a point below 11/c: the tail beyond it is larger
    ph = phi_fi(Lz); Ps = ph / Lz            # Psi(L) <= phi(L)/L
    if kind == "P":
        return FI(0, (Ps * 2).hi)
    if kind == "LC":
        return FI(0, (Ps * 2 * FI.q(Q(7, 10)) + ph * 2 * par["sp"]).hi)
    c = par["c"]; cL2 = c * c * (ph * Lz + Ps)
    if kind == "R":                        # (phi/Phi)(t) <= |t| + 1/|t| <= |t| + 1 for t <= -1 (Gordon); <= 1 for t >= 0
        return FI(0, (cL2 + c * ph * 2 + Ps * 2).hi)
    if kind == "LP":                       # |log Phi(t)| <= t^2/2 + |t| + 2 for t <= -1; <= 1 for t >= 0
        b = cL2 / 2 + c * ph + Ps * 3
        return FI(-b.hi, 0)
    raise ValueError(kind)

def E(kind, par):
    return integrate(kind, par) + tails(kind, par)

def P_of(psi):
    return E("P", {"sp": sqrt_fi(psi)})
def c_of(q):
    return sqrt_fi(q / (FI(ONE) - q))
def R_of(q, alpha, variant=None):
    v = E("R", {"c": c_of(q)}) * alpha
    return v if variant == "no-prefactor" else v / (FI(ONE) - q)
def G_of(alpha, q, psi, variant=None):
    lc = E("LC", {"sp": sqrt_fi(psi)})
    lp = E("LP", {"c": c_of(q)})
    t = alpha * lp
    return -((FI(ONE) - q) * psi) / 2 + lc + (-t if variant == "sign-reversed" else t)

def float_psi(alpha, variant):
    h = 0.005
    zs = [i * h for i in range(-2400, 2401)]
    ws = [h * fexp(-z * z / 2) / fsqrt(2 * fpi) for z in zs]
    Phi = lambda t: 0.5 * ferfc(-t / fsqrt(2))
    ph = lambda t: fexp(-t * t / 2) / fsqrt(2 * fpi)
    def P(psi):
        s = fsqrt(psi); return sum(w * ftanh(s * z) ** 2 for z, w in zip(zs, ws))
    def R(q):
        c = fsqrt(q / (1 - q)); fac = alpha if variant == "no-prefactor" else alpha / (1 - q)
        return fac * sum(w * (ph(c * z) / Phi(c * z)) ** 2 for z, w in zip(zs, ws))
    psi = 2.5
    for _ in range(400):
        nps = R(P(psi))
        if abs(nps - psi) < 1e-15: break
        psi = nps
    return psi

def selftests():
    ok = True
    # the two routes for Psi agree at y = 4 (series) and y = 4 + 2^-F (continued fraction)
    a = _Psi_series(4 * ONE); b = phi_fi(FI(4 * ONE)) * _K(4 * ONE)
    t1 = a.lo - (1 << 60) <= b.hi and b.lo <= a.hi + (1 << 60)
    # pi, log 2, exp(1)
    t2 = abs(PI.f() - 3.141592653589793) < 1e-15 and abs(LOG2.f() - 0.6931471805599453) < 1e-15 \
        and abs(_exp_point(ONE).lo / ONE - 2.718281828459045) < 1e-15
    # Psi(1) = 0.158655253931457051414767454367962077522087033273395609012605... (36 digits)
    p1 = Psi_point(ONE); ref = FI.q(Q("0.158655253931457051414767454367962077"))
    t3 = abs(p1.lo - ref.lo) < 1 << (F - 110)
    # int phi = 1 - 2 Psi(9.5) through the quadrature: log Phi kernel with c -> int over t of log Phi(t) phi(t/c)/c ...
    # simpler: E th(0)^2 = 0
    z = integrate("P", {"sp": FI(0)})
    t4 = z.hi < 1 << (F - 60)
    print("self-test Psi(4): series and continued fraction agree: %s; pi, log 2, e: %s; Psi(1) to 36 digits: %s; "
          "E th(0)^2 = 0: %s" % (t1, t2, t3, t4))
    return ok and t1 and t2 and t3 and t4

def at2024(alpha):
    print("80a program B, [AT2024]: is GD(%s, 1/2) below -log 2?" % alpha)
    if not selftests():
        print("VERDICT FAILED (self-test)"); return
    elh = E("LP", {"c": FI(ONE)})
    print("E log H(u) by quadrature in [%.18f, %.18f]" % (elh.f(), elh.hi / ONE))
    half_log_half = -LOG2 / 2
    gd = FI.q(alpha) * elh + FI.q(Q(1, 2)) + half_log_half
    gap = gd + LOG2
    print("GD(%s, 1/2) + log 2 in [%.12e, %.12e]" % (alpha, gap.f(), gap.hi / ONE))
    ok = gap.hi < 0
    print("row [AT2024] capacity <= %s: %s" % (alpha, ok))
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
    print("80a program B: claimed alpha_star in (%s, %s): G_star(lo) > 0 > G_star(hi)" % (opt["--alpha-lo"], opt["--alpha-hi"]))
    if not selftests():
        print("VERDICT FAILED (self-test)"); return
    qlb, qub = FI.q(Q("0.56394907949")), FI.q(Q("0.56394908030"))
    e = Q(1, 10 ** 12)
    res = {}
    for name in ("--alpha-lo", "--alpha-hi"):
        al = FI.q(Q(opt[name]))
        pt = Q(repr(float_psi(float(Q(opt[name])), variant)))
        s1, s2 = FI.q(pt - e), FI.q(pt + e)
        g1 = R_of(P_of(s1), al, variant) - s1
        g2 = R_of(P_of(s2), al, variant) - s2
        br = g1.lo > 0 and g2.hi < 0
        S = s1.hull(s2)
        qs = P_of(S)
        inside = qs.lo > qlb.hi and qs.hi < qub.lo
        G = G_of(al, qs, S, variant)
        print("alpha = %s: psi~ = %.15f; R(P(psi~ - e)) - (psi~ - e) in [%.3e, %.3e], R(P(psi~ + e)) - (psi~ + e) in "
              "[%.3e, %.3e] (e = 1e-12): bracketed %s" % (opt[name], float(pt), g1.f(), g1.hi / ONE, g2.f(), g2.hi / ONE, br))
        print("   q_star = P(psi_star) in [%.16f, %.16f], inside (q_lb, q_ub): %s" % (qs.f(), qs.hi / ONE, inside))
        print("   G_star in [%.6e, %.6e]" % (G.f(), G.hi / ONE))
        res[name] = (br and inside, G)
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
