#!/usr/bin/env python3
"""Independent replay of the finite certificate behind the 57c row [W2026], U >= 306121/500000 = 0.612242 (Wikstrom,
arXiv:2610.34580v1, Theorem 1.2 and Lemma 6.1; the row is recorded in PR #311), written from the paper's LaTeX source
alone: its verifiers, run records and witness files were neither read nor run.

The statement checked (Lemma 6.1). Put R = 306121/500000 and, for 0 <= s <= 1/3, a = u + iv with |a| <= b(s) =
2/3 - s^2/2, q = |a|^2, the centre c = s(A0 u - i B0 v), A0 = 2/5 - q/3 - s^2/5, B0 = 25/24 - 3q/2 - s^2/4. For a
direction theta, h = Re(c e^-i theta), t = Re(c) sin theta - Im(c) cos theta, p = R + h - t^2/(2R) - t^4/(8R^3). For a
witness z = (alpha + i beta, sqrt3 (gamma + i delta)): A = alpha^2 - beta^2, B = alpha gamma - beta delta,
C = gamma^2 - delta^2, E = alpha delta + beta gamma, F = 2 gamma delta, N = alpha^2 + beta^2 + 3(gamma^2 + delta^2),
d + ie = a e^(3 i theta), V = N + 3s(B cos 2theta - E sin 2theta) - (9/2)(C d - F e), P = 8p^3 V - 4p^2 A - 6pB - 3C.
Lemma 6.1: every (s, a, theta) admits a witness in (2^-16 Z)^4 with A, B, C >= 0 and P < 0. With the paper's
symmetries (c(s, conj a) = conj c, c(s, -a) = -c, and the same for the omitted-value sets) it suffices that
u, v >= 0; the paper's coordinates s = x0/3, a = b(s) x1 e^(i pi x2/2), theta = 2 pi x3, x in [0, 1]^4.

How (our own branch and bound). For a box of x we propose a witness at its midpoint: the matrix M~(p) = H + iG of
the paper's unitary congruence (section 5.2), the 4 x 4 real symmetric [[I + H, -G], [-G, I - H]], an eigenvector of
its smallest eigenvalue by our own Jacobi iteration, scaled and rounded to 2^-16, together with its conjugate and
both with |beta| clipped to |alpha| and |delta| to |gamma| (the eigenvector can lie just outside the cone A, B, C >= 0);
of those passing the exact sign conditions, the one with the smallest T at the midpoint is used ((1, 0, 1, 0) if
none passes). These are proposals only: any witness passing the exact checks is a legitimate one. Then, with that witness fixed, A, B, C >= 0 is checked exactly (integers), and P is
enclosed over the whole box by interval arithmetic: floating point with every operation rounded outward by one ulp
(math.nextafter); the cosines and sines of each box's angle intervals from their midpoints, computed with Python's
decimal at 30 digits (Taylor series, pi by Machin), widened by the half-width (Lipschitz) and by 1e-25. Where this
natural enclosure is not negative, the mean-value form P(m) + sum_j dP/dx_j(box) [-h_j, h_j] is used as well, with
the derivatives by interval forward-mode differentiation (the witness fixed). A box with P < 0 throughout is closed;
otherwise it is bisected along its widest coordinate. The run is complete when no box is left.

Also checked: the paper's side facts used with Lemma 5.2 and the polynomial radial bound (section 5.1): the centre
weights lie in [0, 25/24] on the parameter set, |c| <= 25/108 = L < R, and p >= R - L - L^2/(2R) - L^4/(8R^3) > 0
(exact rationals); the formula for T(rho) = |z|^2 + Re(z^T M~ z) against the matrix itself at 2,000 random points
(floating point, a check of the algebra); and the obstruction
of section 6.3: the root r_* of F(r) = 176r^4 - 96r^3 + 24r^2 - 24r + 3 in (0.6, 0.7) lies in
(0.612242772171, 0.612242772172) (exact sign changes, F convex).

Red controls (--controls), at R = 0.6122428, just above r_*, where Lemma 6.1 is false: at the obstruction point
x = (1, 0, 0, 0) (s = 1/3, a = 0, theta = 0), exactly, 192 R^4 det(I + M) = F(R) and I + M, I - M are positive
definite, so no witness exists there; our witness and enclosure at that point give P > 0 (and P < 0 at the row's R);
and on a small corner box around that point the search, which closes at the row's R, goes past depth 120 at
R = 0.6122428 within the same box budget.
Usage: python3 -I w2026_check.py [--controls]. Python 3.9+, standard library only."""
import math
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction as Fr
from math import inf, nextafter

getcontext().prec = 32
R0 = Fr(306121, 500000)
SQ3 = math.sqrt(3.0)
CORNER_K, CORNER_BUDGET = 10, 2 * 10 ** 6                    # the red control's corner box and its box budget


# --------------------------------------------------------------------- outward-rounded float intervals
def dn(x):
    return nextafter(x, -inf)


def up(x):
    return nextafter(x, inf)


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    @staticmethod
    def c(v):
        """an exact rational or float constant, enclosed (an interval is returned as it is)"""
        if isinstance(v, I):
            return v
        if isinstance(v, Fr):
            f = float(v)
            return I(dn(f), up(f))
        return I(v, v)

    def __add__(a, b):
        b = b if isinstance(b, I) else I.c(b)
        return I(dn(a.lo + b.lo), up(a.hi + b.hi))

    __radd__ = __add__

    def __sub__(a, b):
        b = b if isinstance(b, I) else I.c(b)
        return I(dn(a.lo - b.hi), up(a.hi - b.lo))

    def __rsub__(a, b):
        return I.c(b) - a

    def __neg__(a):
        return I(-a.hi, -a.lo)

    def __mul__(a, b):
        b = b if isinstance(b, I) else I.c(b)
        ps = (a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi)
        return I(dn(min(ps)), up(max(ps)))

    __rmul__ = __mul__

    def sq(a):
        if a.lo >= 0:
            return I(dn(a.lo * a.lo), up(a.hi * a.hi))
        if a.hi <= 0:
            return I(dn(a.hi * a.hi), up(a.lo * a.lo))
        return I(0.0, up(max(a.lo * a.lo, a.hi * a.hi)))

    def __truediv__(a, b):
        b = b if isinstance(b, I) else I.c(b)
        assert b.lo > 0 or b.hi < 0
        qs = (a.lo / b.lo, a.lo / b.hi, a.hi / b.lo, a.hi / b.hi)
        return I(dn(min(qs)), up(max(qs)))


# --------------------------------------------------------------------- trigonometry at midpoints, in decimal
def _pi():
    def atn(n):
        x = D(1) / n
        x2, t, s, k = x * x, x, x, 1
        while True:
            t *= -x2
            k += 2
            if abs(t) < D("1e-40"):
                return s
            s += t / k
    return 16 * atn(5) - 4 * atn(239)


PI = _pi()
PI_I = I(dn(float(PI)), up(float(PI)))   # pi enclosed: the decimal is exact to 30 digits, float(PI) rounds it
EPS_TRIG = 1e-25


def dcos_sin(m):
    """cos and sin of a float m (any real), by Taylor series in decimal after reduction to [-pi, pi]."""
    x = D(m)
    k = (x / (2 * PI)).to_integral_value()
    x -= 2 * PI * k
    x2 = x * x
    c, s, tc, ts, n = D(1), x, D(1), x, 0
    while True:
        n += 2
        tc = -tc * x2 / (n * (n - 1))
        ts = -ts * x2 / ((n + 1) * n)
        c += tc
        s += ts
        if abs(tc) < D("1e-31") and abs(ts) < D("1e-31"):
            return c, s


def cos_sin_iv(ang):
    """enclosures of cos and sin over the angle interval ang = [lo, hi] (floats)"""
    m = (ang.lo + ang.hi) / 2
    w = up(max(m - ang.lo, ang.hi - m)) + EPS_TRIG
    c, s = dcos_sin(m)
    cf, sf = float(c), float(s)
    cl, ch = max(-1.0, dn(dn(cf) - w)), min(1.0, up(up(cf) + w))
    sl, sh = max(-1.0, dn(dn(sf) - w)), min(1.0, up(up(sf) + w))
    return I(cl, ch), I(sl, sh)


# --------------------------------------------------------------------- interval forward-mode derivatives
class ID:
    """an interval value with its gradient (four intervals) over the box: for the mean-value enclosure"""
    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v, self.g = v, g

    @staticmethod
    def c(x):
        z = I(0.0, 0.0)
        return ID(I.c(x) if not isinstance(x, I) else x, (z, z, z, z))

    def __add__(a, b):
        b = b if isinstance(b, ID) else ID.c(b)
        return ID(a.v + b.v, tuple(x + y for x, y in zip(a.g, b.g)))

    __radd__ = __add__

    def __sub__(a, b):
        b = b if isinstance(b, ID) else ID.c(b)
        return ID(a.v - b.v, tuple(x - y for x, y in zip(a.g, b.g)))

    def __rsub__(a, b):
        return ID.c(b) - a

    def __neg__(a):
        return ID(-a.v, tuple(-x for x in a.g))

    def __mul__(a, b):
        b = b if isinstance(b, ID) else ID.c(b)
        return ID(a.v * b.v, tuple(x * b.v + a.v * y for x, y in zip(a.g, b.g)))

    __rmul__ = __mul__

    def sq(a):
        two = I(2.0, 2.0)
        return ID(a.v.sq(), tuple(two * a.v * x for x in a.g))

    def __truediv__(a, b):
        b = b if isinstance(b, ID) else ID.c(b)
        q = a.v / b.v
        b2 = b.v.sq()
        return ID(q, tuple((x * b.v - a.v * y) / b2 for x, y in zip(a.g, b.g)))


def cos_sin_ad(ang):
    """ang an ID whose value is an angle interval; cos and sin with their gradients"""
    c, s = cos_sin_iv(ang.v)
    return ID(c, tuple(-(s * g) for g in ang.g)), ID(s, tuple(c * g for g in ang.g))


# --------------------------------------------------------------------- the problem
def params(box, R, ad=False):
    """quantities over a box [x0, x1, x2, x3] (pairs of floats in [0, 1]); with ad=True, with their gradients"""
    if ad:
        z = I(0.0, 0.0)
        X = [ID(I(lo, hi), tuple(I(1.0, 1.0) if k == j else z for k in range(4))) for j, (lo, hi) in enumerate(box)]
        C = ID.c
        trig = cos_sin_ad
    else:
        X = [I(lo, hi) for lo, hi in box]
        C = I.c
        trig = cos_sin_iv
    x0, x1, x2, x3 = X
    s = x0 * C(Fr(1, 3))
    bs = C(Fr(2, 3)) - s.sq() * C(Fr(1, 2))
    mod = bs * x1                                          # |a|
    phi = x2 * C(PI_I) * C(Fr(1, 2))                       # arg a in [0, pi/2]
    th = x3 * C(PI_I) * C(Fr(2))                           # theta
    cphi, sphi = trig(phi)
    cth, sth = trig(th)
    c2, s2 = trig(th * C(Fr(2)))
    cps, sps = trig(phi + th * C(Fr(3)))                   # arg of a e^(3 i theta)
    u, v = mod * cphi, mod * sphi
    q = mod.sq()
    A0 = C(Fr(2, 5)) - q * C(Fr(1, 3)) - s.sq() * C(Fr(1, 5))
    B0 = C(Fr(25, 24)) - q * C(Fr(3, 2)) - s.sq() * C(Fr(1, 4))
    rec, imc = s * A0 * u, -(s * B0 * v)
    h = rec * cth + imc * sth
    t = rec * sth - imc * cth
    Rv = C(R)
    t2 = t.sq()
    p = Rv + h - t2 / (Rv * C(Fr(2))) - t2.sq() / (Rv * Rv * Rv * C(Fr(8)))
    d, e = mod * cps, mod * sps
    return s, p, c2, s2, d, e, C


def P_value(box, wit, R, ad=False):
    s, p, c2, s2, d, e, C = params(box, R, ad)
    al, be, ga, de = (C(Fr(x, 1 << 16)) for x in wit)
    A = al.sq() - be.sq()
    B = al * ga - be * de
    Cc = ga.sq() - de.sq()
    E = al * de + be * ga
    F = C(Fr(2)) * ga * de
    N = al.sq() + be.sq() + C(Fr(3)) * (ga.sq() + de.sq())
    V = N + C(Fr(3)) * s * (B * c2 - E * s2) - C(Fr(9, 2)) * (Cc * d - F * e)
    p2 = p.sq()
    return C(Fr(8)) * p2 * p * V - C(Fr(4)) * p2 * A - C(Fr(6)) * p * B - C(Fr(3)) * Cc


def P_enclosure(box, wit, R, centred=True):
    """the natural enclosure of P over the box, intersected with the mean-value form
    P(m) + sum_j dP/dx_j(box) [-h_j, h_j] (derivatives by interval forward mode, the witness fixed)."""
    nat = P_value(box, wit, R)
    if not centred or nat.hi < 0:
        return nat
    m = [((lo + hi) / 2, (lo + hi) / 2) for lo, hi in box]
    pm = P_value(m, wit, R)
    g = P_value(box, wit, R, ad=True).g
    hi = pm.hi
    for (lo_, hi_), dj in zip(box, g):
        hw = up(max((hi_ - lo_) / 2, 0.0))
        hi = up(hi + up(max(abs(dj.lo), abs(dj.hi)) * hw))
    return I(nat.lo, min(nat.hi, hi))


def jacobi_min_vec(S):
    """eigenvector of the smallest eigenvalue of a small real symmetric matrix (floats; a proposal only)"""
    n = len(S)
    a = [row[:] for row in S]
    V = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(60):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-24:
            break
        for p_ in range(n):
            for q_ in range(p_ + 1, n):
                if abs(a[p_][q_]) < 1e-300:
                    continue
                th = (a[q_][q_] - a[p_][p_]) / (2 * a[p_][q_])
                t = (1.0 if th >= 0 else -1.0) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1)
                s_ = t * c
                for k in range(n):
                    akp, akq = a[k][p_], a[k][q_]
                    a[k][p_], a[k][q_] = c * akp - s_ * akq, s_ * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p_][k], a[q_][k]
                    a[p_][k], a[q_][k] = c * apk - s_ * aqk, s_ * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p_], V[k][q_]
                    V[k][p_], V[k][q_] = c * vkp - s_ * vkq, s_ * vkp + c * vkq
    i = min(range(n), key=lambda k: a[k][k])
    return [V[k][i] for k in range(n)]


def signs_ok(w):
    a_, b_, g_, d_ = w
    return a_ * a_ - b_ * b_ >= 0 and a_ * g_ - b_ * d_ >= 0 and g_ * g_ - d_ * d_ >= 0


def propose(box, R):
    x = [(lo + hi) / 2 for lo, hi in box]
    s = x[0] / 3
    mod = (2 / 3 - s * s / 2) * x[1]
    phi, th = math.pi * x[2] / 2, 2 * math.pi * x[3]
    u, v = mod * math.cos(phi), mod * math.sin(phi)
    q = mod * mod
    A0 = 2 / 5 - q / 3 - s * s / 5
    B0 = 25 / 24 - 1.5 * q - s * s / 4
    rec, imc = s * A0 * u, -s * B0 * v
    Rf = float(R)
    h = rec * math.cos(th) + imc * math.sin(th)
    t = rec * math.sin(th) - imc * math.cos(th)
    rho = Rf + h - t * t / (2 * Rf) - t ** 4 / (8 * Rf ** 3)
    m11 = complex(-1 / (2 * rho), 0)
    m12 = SQ3 * (s * complex(math.cos(2 * th), math.sin(2 * th)) / 2 - 1 / (8 * rho * rho))
    m22 = -1.5 * mod * complex(math.cos(phi + 3 * th), math.sin(phi + 3 * th)) - 1 / (8 * rho ** 3)
    H = [[m11.real, m12.real], [m12.real, m22.real]]
    G = [[m11.imag, m12.imag], [m12.imag, m22.imag]]
    S = [[(1.0 if i == j else 0.0) + H[i][j] for j in range(2)] + [-G[i][j] for j in range(2)] for i in range(2)]
    S += [[-G[i][j] for j in range(2)] + [(1.0 if i == j else 0.0) - H[i][j] for j in range(2)] for i in range(2)]
    x1, x2, y1, y2 = jacobi_min_vec(S)
    al, be, ga, de = x1, y1, x2 / SQ3, y2 / SQ3

    def T_mid(w):                                          # T(rho)/|z|^2 at the midpoint, floats (for choosing)
        a_, b_, g_, d_ = (c / 65536 for c in w)
        z1, z2 = complex(a_, b_), SQ3 * complex(g_, d_)
        n = abs(z1) ** 2 + abs(z2) ** 2
        return (n + (m11 * z1 * z1 + 2 * m12 * z1 * z2 + m22 * z2 * z2).real) / n

    def rnd(vals):
        sc = max(abs(c) for c in vals)
        return tuple(int(round(c / sc * (1 << 16))) for c in vals)

    # candidates: the eigenvector, its conjugate, and both with |beta| clipped to |alpha| and |delta| to |gamma|
    # (the eigenvector can sit just outside the cone A, B, C >= 0); any candidate passing the exact sign conditions
    # is a legitimate witness, and the one with the smallest T at the midpoint is proposed
    cands = []
    for b2, d2 in ((be, de), (-be, -de)):
        cands.append((al, b2, ga, d2))
        cands.append((al, math.copysign(min(abs(b2), abs(al)), b2), ga, math.copysign(min(abs(d2), abs(ga)), d2)))
    ws = [rnd(c) for c in cands if max(abs(x) for x in c) > 0]
    ws = [w for w in ws if signs_ok(w)]
    if ws:
        return min(ws, key=T_mid)
    return (1 << 16, 0, 1 << 16, 0)


def certify(R, max_boxes, max_depth=120, log=print, root=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0))):
    stack = [(root, 0)]
    visited = closed = deepest = 0
    while stack:
        box, depth = stack.pop()
        visited += 1
        deepest = max(deepest, depth)
        if visited > max_boxes or depth > max_depth:
            return False, visited, closed, deepest, box
        w = propose(box, R)
        assert signs_ok(w)                                   # checked exactly, in integers
        if P_enclosure(box, w, R).hi < 0:
            closed += 1
            continue
        j = max(range(4), key=lambda k: box[k][1] - box[k][0])
        lo, hi = box[j]
        mid = (lo + hi) / 2
        for half in ((lo, mid), (mid, hi)):
            nb = list(box)
            nb[j] = half
            stack.append((tuple(nb), depth + 1))
        if log and visited % 200000 == 0:
            log(f"    ... {visited:,} boxes visited, {len(stack):,} open")
    return True, visited, closed, deepest, None


def side_facts(R):
    ok = True
    L = Fr(25, 108)
    # weights on the parameter set: q <= b(s)^2 <= 4/9 and s <= 1/3
    A0_min = Fr(2, 5) - Fr(4, 9) / 3 - Fr(1, 9) / 5
    B0_min = Fr(25, 24) - Fr(3, 2) * Fr(4, 9) - Fr(1, 9) / 4
    ok &= A0_min > 0 and B0_min > 0 and Fr(2, 5) <= Fr(25, 24)
    ok &= Fr(1, 3) * Fr(25, 24) * Fr(2, 3) == L and L < R
    pmin = R - L - L * L / (2 * R) - L ** 4 / (8 * R ** 3)
    ok &= pmin > 0
    return ok, A0_min, B0_min, pmin


def obstruction():
    def Fp(r):
        return 176 * r ** 4 - 96 * r ** 3 + 24 * r ** 2 - 24 * r + 3
    lo, hi = Fr(612242772171, 10 ** 12), Fr(612242772172, 10 ** 12)
    d1 = lambda r: 704 * r ** 3 - 288 * r ** 2 + 48 * r - 24
    # F'' = 2112 r^2 - 576 r + 48 has discriminant 576^2 - 4*2112*48 < 0, so F is convex; F'(0.6) > 0
    convex = 576 ** 2 - 4 * 2112 * 48 < 0
    return Fp(lo) < 0 < Fp(hi) and convex and d1(Fr(3, 5)) > 0 and Fp(Fr(3, 5)) < 0 < Fp(Fr(7, 10))


def obstruction_point(R):
    """the obstruction point of section 6.3, x = (1, 0, 0, 0): s = 1/3, a = 0, theta = 0. There c = 0 and p = R, and
    for z = x + iy the form is T = x^T (I + M) x + y^T (I - M) y with the real symmetric M = M_{1/3,0}(R) (entries
    -1/(2R), sqrt3 (1/6 - 1/(8R^2)) and -1/(8R^3)). Exact: whether I + M and I - M are both positive definite, in
    which case T > 0 and P = 8R^3 T > 0 for every witness, so no witness exists there."""
    m11, m22 = -1 / (2 * R), -1 / (8 * R ** 3)
    b2 = 3 * (Fr(1, 6) - 1 / (8 * R ** 2)) ** 2                # the off-diagonal entry, squared
    det_plus = (1 + m11) * (1 + m22) - b2
    plus = 1 + m11 > 0 and det_plus > 0
    minus = 1 - m11 > 0 and (1 - m11) * (1 - m22) - b2 > 0
    return plus and minus, det_plus


def formula_check(n=2000):
    """T(rho) of the paper's formula against |z|^2 + Re(z^T M~(rho) z) computed from the matrix (floats)."""
    import random
    rnd = random.Random(7)
    worst = 0.0
    for _ in range(n):
        s = rnd.uniform(0, 1 / 3)
        a = complex(rnd.uniform(-0.6, 0.6), rnd.uniform(-0.6, 0.6))
        th, rho = rnd.uniform(0, 2 * math.pi), rnd.uniform(0.3, 0.9)
        al, be, ga, de = (rnd.uniform(-1, 1) for _ in range(4))
        z = [complex(al, be), SQ3 * complex(ga, de)]
        e2, e3 = complex(math.cos(2 * th), math.sin(2 * th)), complex(math.cos(3 * th), math.sin(3 * th))
        M = [[-1 / (2 * rho), SQ3 * (s * e2 / 2 - 1 / (8 * rho ** 2))],
             [SQ3 * (s * e2 / 2 - 1 / (8 * rho ** 2)), -1.5 * a * e3 - 1 / (8 * rho ** 3)]]
        direct = abs(z[0]) ** 2 + abs(z[1]) ** 2 + sum(z[i] * M[i][j] * z[j] for i in range(2) for j in range(2)).real
        A, B, C = al * al - be * be, al * ga - be * de, ga * ga - de * de
        E, F, N = al * de + be * ga, 2 * ga * de, al * al + be * be + 3 * (ga * ga + de * de)
        d, e = (a * e3).real, (a * e3).imag
        V = N + 3 * s * (B * e2.real - E * e2.imag) - 4.5 * (C * d - F * e)
        T = V - A / (2 * rho) - 3 * B / (4 * rho ** 2) - 3 * C / (8 * rho ** 3)
        worst = max(worst, abs(T - direct))
    return worst


def red_controls():
    """the red controls at R = 0.6122428 (just above r_*), where Lemma 6.1 is false; True when every one refuses"""
    ok = True
    Rc = Fr(6122428, 10 ** 7)
    print("Red controls at R = 0.6122428, just above r_* = 0.6122427721..., where Lemma 6.1 is false (each must "
          "refuse):")
    nowit_c, det_c = obstruction_point(Rc)
    nowit_0, det_0 = obstruction_point(R0)
    Fp = lambda r: 176 * r ** 4 - 96 * r ** 3 + 24 * r ** 2 - 24 * r + 3
    same = 192 * Rc ** 4 * det_c == Fp(Rc) and 192 * R0 ** 4 * det_0 == Fp(R0)
    good = nowit_c and not nowit_0 and same
    print(f"  at the obstruction point x = (1, 0, 0, 0) (s = 1/3, a = 0, theta = 0; 192 R^4 det(I + M) = F(R) "
          f"exactly: {same}): I + M and I - M are positive definite, so P > 0 for every witness (exact): "
          f"{nowit_c}; at the row's R, det(I + M) < 0: {not nowit_0}")
    ok &= good
    xs = ((1.0, 1.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0))
    e0 = P_enclosure(xs, propose(xs, R0), R0)
    ec = P_enclosure(xs, propose(xs, Rc), Rc)
    good = e0.hi < 0 and ec.lo > 0
    print(f"  our witness and enclosure at that point: P in [{ec.lo:.6e}, {ec.hi:.6e}], not negative (at the "
          f"row's R: P <= {e0.hi:.6e} < 0): {'refused, as expected' if good else 'NOT REFUSED'}")
    ok &= good
    k = CORNER_K
    Q = ((1.0 - 2.0 ** -k, 1.0), (0.0, 2.0 ** -k), (0.0, 1.0), (0.0, 2.0 ** -(k + 1)))
    g0, v0, _, d0, _ = certify(R0, CORNER_BUDGET, root=Q, log=None)
    gc, vc, _, dc, box = certify(Rc, CORNER_BUDGET, root=Q, log=None)
    where = "" if box is None else f" near x = {tuple(float(f'{(a + b) / 2:.6g}') for a, b in box)}"
    good = g0 and not gc and dc > 120
    print(f"  the corner box [1 - 2^-{k}, 1] x [0, 2^-{k}] x [0, 1] x [0, 2^-{k + 1}] around it, the same budget "
          f"of {CORNER_BUDGET:,} boxes: at the row's R it closes ({v0:,} boxes, deepest {d0}); at R = 0.6122428 "
          f"the search passes depth 120 after {vc:,} boxes{where}: "
          f"{'refused, as expected' if good else 'NOT REFUSED'}")
    ok &= good
    return ok


def main():
    t0 = time.time()
    controls = "--controls" in sys.argv
    ok = True
    print("57c row [W2026]: U >= 306121/500000 = 0.612242 (Wikstrom, Theorem 1.2), through Lemma 6.1")
    good, A0m, B0m, pmin = side_facts(R0)
    print(f"  centre weights >= {float(A0m):.4f} and {float(B0m):.4f} (> 0) and <= 25/24; |c| <= 25/108 < R; "
          f"p >= {float(pmin):.6f} > 0 (exact): {good}")
    ok &= good
    good = obstruction()
    print(f"  the obstruction: F(r) = 176r^4 - 96r^3 + 24r^2 - 24r + 3 is convex and changes sign in "
          f"(0.612242772171, 0.612242772172), its only root in (0.6, 0.7): {good}")
    ok &= good
    worst = formula_check()
    good = worst < 1e-9
    print(f"  T(rho) by the paper's formula equals |z|^2 + Re(z^T M~ z) at 2,000 random points (largest difference "
          f"{worst:.1e}): {good}")
    ok &= good
    print("  Lemma 6.1 over x in [0, 1]^4, our own witnesses and enclosures:")
    good, visited, closed, deepest, _ = certify(R0, 10 ** 8)
    print(f"  complete: {good}; {visited:,} boxes visited, {closed:,} closed (terminal), deepest {deepest} "
          f"(the paper: 187,263 visited, 93,632 terminal, depth 49); {time.time() - t0:.0f} s")
    ok &= good
    if controls:
        ok &= red_controls()
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
