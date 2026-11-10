#!/usr/bin/env python3
"""Replay of the conformal-radius half of the 57c row [vS2026] (W. van Staden), rho(D) >= 1.523726208329 (Proposition 4.1),
written from the paper's LaTeX source (paper/bloch_slit_discs.tex, Section 4: Lemma 4.2, Lemma 4.3, Proposition 4.4,
Lemmas 4.5 to 4.8 and the four steps of Section 4.4) and the data file certify/density_G_n64.json alone (its
coefficients equal the paper's appendix table; we check that). conformal_certificate.py and every other program of the
repository were neither read nor run; our parameters (the ellipses, M, N) are our own choices.

Comparison function (Lemma 4.3, re-derived): in w = log z, U(w) = (L - Re w) + sum_j int Phi(w - w_j(y)) psi_j(y) dmu(y),
Phi(w) = -log|2 sinh(w/2)|, w_j(x) = L + a_j x + i(theta_j + b1_j (a_j x)^2 + b2_j (a_j x)^4), psi_j odd with the
class's Chebyshev coefficients; U is harmonic off the arcs, U = 0 on Re w = L by the reflection, and
U = -Re w + c + o(1) with c = L - (pi/4) sum_j a_j c_{j,1}. If |U| <= delta on the slits then log rho(D) >= c - delta
(Lemma 4.2 with v = u - delta). On slit i, f_i(x) = U(w_i(x)) is odd and equals (4.4)
   f_i(x) = -a_i x + sum_n c_{i,n} (pi/n) T_n(x) - int log|K_i(x,y)| psi_i dmu - sum_{j != i} int log|G_ij(x,y)| psi_j dmu,
G_ij = 2 sinh((w_i(x) - w_j(y))/2), K_i = G_ii/(x - y) (the log kernel's Chebyshev identity, Lemma 4.5, re-derived).
By the two symmetries one slit per class suffices (slits 0, 3, 6, 12, 18 in our numbering).

Bounding sup |f_i| (Section 4.4, our implementation):
 1. kernel bounds, x in the rectangle Q_b (contains the Bernstein ellipse E_r, r = b + sqrt(1 + b^2)) and y in [-1, 1]:
    lambda <= |G~_ij| <= Lambda by branch and bound in complex interval arithmetic (|2 sinh(z/2)|^2 =
    4 sinh^2(Re z/2) + 4 sin^2(Im z/2); K_i = Delta shc(tau Delta), tau = (x - y)/2, by its series); A_ij =
    max |log| of the bounds; B_i = a_i (r + 1/r)/2 + sum_n |c_{i,n}| (pi/2n)(r^n + r^-n) + pi sum_j ||psi_j|| A_ij;
    then the Chebyshev coefficients of f_i are at most 4 B_i r^-k (Lemma 4.6);
 2. the same for x in [-1, 0] and y in Q_b' (ellipse r'), giving A'_ij; the N-point Gauss-Chebyshev rule then errs by at
    most 8 pi A' ||psi_j|| r'^-(2N-1-63)/(r'-1) (Lemmas 4.6, 4.7); N_ij is chosen to make this at most 10^-13;
 3. f_i at the negative zeros of T_M in interval arithmetic (vs_nodes.c), M chosen with 8 B r^(1-M)/(r-1) <= 10^-9;
    positive nodes by oddness;
 4. delta_i = sum_q |gamma_q| + 8 B r^(1-M)/(r-1) (Lemma 4.8) with gamma the discrete Chebyshev transform of the node
    values. We transform the interval midpoints (exactly, in rationals with enclosed cosines) and add
    sqrt(2M) max_k (node half-width): the transform is orthogonal up to scale, so sum_q |gamma_q - gamma~_q| <=
    sqrt(M) ||gamma - gamma~||_2 <= sqrt(2M) max |f_k - f~_k|.
Arithmetic: rationals with series bounds for the constants (pi, logs); binary64 intervals widened by one ulp per
operation (IEEE operations are correctly rounded) and by four ulps for libm's sin, cos, sinh, log, exp (taken to be
within one ulp) in the branch and bound and in vs_nodes.c.

Usage: python3 vs_conformal_check.py DENSITY_JSON [--tex=PAPER_TEX] [--controls]   (needs bin/vs_nodes, from vs_nodes.c;
       the ellipse sizes b, b' are chosen from a fixed list by trying the rigorous branch and bound for each; about 10 s)
"""
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
SHA_DENSITY = "c1164ede5da9ca075398c25050290ea5352c3e583d5f6ff240ec97145f45da73"
R = F("5.001117628756668")
CLASSES = [(F(0), F(1), F(0), F(0)),
           (None, F("1.9999999974092333"), F(0), F(0)),
           (F("0.029930188727091613"), F("2.789342409403121"), F("0.9026921410641155"), F("-0.3301459494660587")),
           (F("1.0157094972537686"), F("3.473381248030611"), F("-1.2656124700253997"), F("0.9794269477754911")),
           (F("0.5757776442359243"), F("4.342874184977691"), F("0.03225521024585182"), F("2.9360081049550737"))]
REPS = [0, 3, 6, 12, 18]
B_CANDS = (0.5, 0.35, 0.25, 0.18, 0.13, 0.09, 0.065, 0.045, 0.032, 0.022, 0.016, 0.012, 0.009, 0.0065, 0.0045, 0.003)
CLAIM = F("1.523726208329")


def need(c, msg):
    if not c:
        raise AssertionError(msg)


# ----------------------------------------------------------------------------- exact rational constants (outward)
P = 200
ONE = 1 << P


def rd(x):
    return F((x.numerator * ONE) // x.denominator, ONE)


def ru(x):
    return F(-((-x.numerator * ONE) // x.denominator), ONE)


def atan_inv(n, terms=90):
    s, prev = F(0), F(0)
    for k in range(terms + 1):
        prev = s
        s += F((-1) ** k, (2 * k + 1) * n ** (2 * k + 1))
    return min(s, prev), max(s, prev)


_a, _b = atan_inv(5), atan_inv(239)
PI = (rd(16 * _a[0] - 4 * _b[1]), ru(16 * _a[1] - 4 * _b[0]))


def log_q(y, terms=200):
    """log(y), rational y > 1 (y <= 6): 2 artanh(u), u = (y-1)/(y+1), with the geometric tail."""
    u = (y - 1) / (y + 1)
    s, p = F(0), u
    for k in range(terms + 1):
        s += p / (2 * k + 1)
        p *= u * u
    tail = p / ((2 * terms + 3) * (1 - u * u))
    return rd(2 * s), ru(2 * (s + tail))


def exp_q(x, terms=80):
    """exp(x), |x| <= 2: Taylor with the Lagrange remainder (e^2 < 8)."""
    s, term = F(0), F(1)
    for k in range(terms + 1):
        s += term
        term = term * x / (k + 1)
    rem = abs(term) * 8
    return rd(s - rem), ru(s + rem)


def iv_mul(a, b):
    p = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return min(p), max(p)


# ----------------------------------------------------------------------------- binary64 intervals
dn = lambda x: math.nextafter(x, -math.inf)
up = lambda x: math.nextafter(x, math.inf)


def dn4(x):
    return dn(dn(dn(dn(x))))


def up4(x):
    return up(up(up(up(x))))


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi

    @staticmethod
    def q(fr):
        """an exact rational, or a rational interval (lo, hi), enclosed in binary64"""
        if isinstance(fr, tuple):
            return I(dn(float(fr[0])), up(float(fr[1])))
        return I(dn(float(fr)), up(float(fr)))

    def __add__(s, o):
        o = o if isinstance(o, I) else I(o)
        return I(dn(s.lo + o.lo), up(s.hi + o.hi))

    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, I) else I(o)
        return I(dn(s.lo - o.hi), up(s.hi - o.lo))

    def __rsub__(s, o):
        return I(o) - s

    def __neg__(s):
        return I(-s.hi, -s.lo)

    def __mul__(s, o):
        o = o if isinstance(o, I) else I(o)
        p = (s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi)
        return I(dn(min(p)), up(max(p)))

    __rmul__ = __mul__

    def sqr(s):
        a, b = abs(s.lo), abs(s.hi)
        if s.lo <= 0 <= s.hi:
            return I(0.0, up(max(a, b) ** 2))
        m, M = min(a, b), max(a, b)
        return I(dn(m * m), up(M * M))

    def div(s, o):
        o = o if isinstance(o, I) else I(o)
        need(o.lo > 0 or o.hi < 0, "division by an interval containing 0")
        p = (s.lo / o.lo, s.lo / o.hi, s.hi / o.lo, s.hi / o.hi)
        return I(dn(min(p)), up(max(p)))

    def mag(s):
        return max(abs(s.lo), abs(s.hi))

    def mig(s):
        return 0.0 if s.lo <= 0 <= s.hi else min(abs(s.lo), abs(s.hi))


PI_I = I.q(PI)


def isinh(a):
    return I(dn4(math.sinh(a.lo)), up4(math.sinh(a.hi)))


def isin2(a):
    if not (a.hi - a.lo < 1.0):
        return I(0.0, 1.0)
    mid = 0.5 * (a.lo + a.hi)
    k = math.floor(mid / (0.5 * PI_I.lo) + 0.5)
    e = (k * (0.5 * PI_I.lo), k * (0.5 * PI_I.hi))
    elo, ehi = dn(dn(min(e))), up(up(max(e)))
    v1, v2 = math.sin(a.lo) ** 2, math.sin(a.hi) ** 2
    lo, hi = min(v1, v2), max(v1, v2)
    if ehi >= a.lo and elo <= a.hi:
        if k % 2 == 0:
            lo = 0.0
        else:
            hi = 1.0
    return I(max(0.0, dn(lo * (1 - 4e-15))), min(1.0, up(hi * (1 + 4e-15))))


class C:
    """rectangular complex interval"""
    __slots__ = ("re", "im")

    def __init__(self, re, im):
        self.re, self.im = re, im

    def __add__(s, o):
        return C(s.re + o.re, s.im + o.im)

    def __sub__(s, o):
        return C(s.re - o.re, s.im - o.im)

    def __mul__(s, o):
        if isinstance(o, I):
            return C(s.re * o, s.im * o)
        return C(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)

    def abs2(s):
        return s.re.sqr() + s.im.sqr()


def shc_c(t, terms=26):
    """shc(t) = sinh(t)/t = sum t^(2k)/(2k+1)! for |t| <= 3, complex interval, with the tail bound
    |t|^(2K+2)/(2K+3)! * 2 added to both parts."""
    t2 = t * t
    s = C(I(1.0), I(0.0))
    term = C(I(1.0), I(0.0))
    for k in range(1, terms + 1):
        term = term * t2
        d = float((2 * k) * (2 * k + 1))
        term = C(term.re.div(I(d)), term.im.div(I(d)))
        s = s + term
    m = math.sqrt(t2.re.mag() ** 2 + t2.im.mag() ** 2) * (1 + 1e-15)     # >= |t|^2
    need(m <= 9.0, "shc series range")
    r = 2.0
    for k in range(terms + 1):
        r *= m
    for k in range(1, 2 * terms + 4):
        r /= k
    r = up(r * 1.0001)
    return C(I(dn(s.re.lo - r), up(s.re.hi + r)), I(dn(s.im.lo - r), up(s.im.hi + r)))


class Slit:
    def __init__(self, k):
        if k < 6:
            cl, r, eps = k // 3, k % 3, 1
        else:
            cl, m = 2 + (k - 6) // 6, (k - 6) % 6
            r, eps = m // 2, (1 if m % 2 == 0 else -1)
        th, t, b1, b2 = CLASSES[cl]
        self.cl = cl
        th_iv = (PI[0] / 3, PI[1] / 3) if th is None else (th, th)
        if eps < 0:
            th_iv = (-th_iv[1], -th_iv[0])
        th_iv = (th_iv[0] + 2 * r * PI[0] / 3, th_iv[1] + 2 * r * PI[1] / 3)
        a = log_q(R / t)
        self.a_q, self.th_q = a, th_iv
        self.b1, self.b2 = eps * b1, eps * b2
        self.B1_q = (min(self.b1 * a[0] ** 2, self.b1 * a[1] ** 2), max(self.b1 * a[0] ** 2, self.b1 * a[1] ** 2))
        self.B2_q = (min(self.b2 * a[0] ** 4, self.b2 * a[1] ** 4), max(self.b2 * a[0] ** 4, self.b2 * a[1] ** 4))
        self.A, self.TH, self.B1, self.B2 = I.q(a), I.q(th_iv), I.q(self.B1_q), I.q(self.B2_q)
        self.af = float(a[0])

    def w_c(self, z):
        """w(z) - L for complex z = p + iq: (a p - 2 B1 p q - B2 Im z^4) + i (Theta + a q + B1 (p^2 - q^2) + B2 Re z^4)"""
        p, q = z.re, z.im
        re2 = p.sqr() - q.sqr()
        im2 = 2 * (p * q)
        re4 = re2.sqr() - im2.sqr()
        im4 = 2 * (re2 * im2)
        return C(self.A * p - self.B1 * im2 - self.B2 * im4, self.TH + self.A * q + self.B1 * re2 + self.B2 * re4)


SLITS = [Slit(k) for k in range(24)]


def g2_enclosure(i, j, x, y):
    """|G~_ij(x, y)|^2 on a box: x, y complex intervals (C)."""
    si = SLITS[i]
    if j != i:
        d = si.w_c(x) - SLITS[j].w_c(y)
        sh = isinh(d.re.div(I(2.0)))
        return I(4.0) * (sh.sqr() + isin2(d.im.div(I(2.0))))
    # K_i = Delta shc(tau Delta), tau = (x - y)/2, Delta = a + i (B1 (x + y) + B2 (x^3 + x^2 y + x y^2 + y^3))
    tau = C((x.re - y.re).div(I(2.0)), (x.im - y.im).div(I(2.0)))
    s = x + y
    x2, y2 = x * x, y * y
    cub = (x2 * x) + (x2 * y) + (x * y2) + (y2 * y)
    beta = C(si.B1 * s.re + si.B2 * cub.re, si.B1 * s.im + si.B2 * cub.im)       # Delta = a + i beta
    delta = C(si.A - beta.im, beta.re)
    t = tau * delta
    sh = shc_c(t)
    return delta.abs2() * sh.abs2()


def bnb(i, j, box, which, max_boxes=400000):
    """lambda^2 <= |G~|^2 <= Lambda^2 over box: which = 'x' (x = p + iq complex, y real) or 'y' (x real, y complex).
    box = ((p0, p1), (q0, q1), (s0, s1)) with (s0, s1) the real variable. Returns (lo, hi, boxes) or None."""
    stack = [box]
    lo, hi, n = math.inf, 0.0, 0
    while stack:
        b = stack.pop()
        n += 1
        if n > max_boxes:
            return None
        (p0, p1), (q0, q1), (s0, s1) = b
        z = C(I(p0, p1), I(q0, q1))
        r_ = C(I(s0, s1), I(0.0))
        e = g2_enclosure(i, j, z, r_) if which == "x" else g2_enclosure(i, j, r_, z)
        if e.lo > 0 and e.hi <= 16 * e.lo:
            lo, hi = min(lo, e.lo), max(hi, e.hi)
            continue
        if max(p1 - p0, s1 - s0, q1 - q0) < 1e-9:
            return None
        # bisect the widest of the three (q usually thin)
        w = [(p1 - p0, 0), (q1 - q0, 1), (s1 - s0, 2)]
        _, d = max(w)
        if d == 0:
            m = 0.5 * (p0 + p1)
            stack += [((p0, m), (q0, q1), (s0, s1)), ((m, p1), (q0, q1), (s0, s1))]
        elif d == 1:
            m = 0.5 * (q0 + q1)
            stack += [((p0, p1), (q0, m), (s0, s1)), ((p0, p1), (m, q1), (s0, s1))]
        else:
            m = 0.5 * (s0 + s1)
            stack += [((p0, p1), (q0, q1), (s0, m)), ((p0, p1), (q0, q1), (m, s1))]
    return lo, hi, n


# ----------------------------------------------------------------------------- the certificate
def run(density_path, coef_override=None, log=print, tol_x=F(1, 10 ** 9), tol_q=F(1, 10 ** 13)):
    t0 = time.time()
    raw = open(density_path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    d = json.loads(raw)
    co = {int(k): [F(v) for v in vs] for k, vs in d["coefficients"].items()}
    if coef_override:
        for (cl, n), v in coef_override.items():
            co[cl][n] = v
    need(all(co[c][n] == 0 for c in co for n in range(0, 64, 2)), "an even-indexed coefficient is not zero")
    norm = {c: sum(abs(v) for v in co[c]) for c in co}               # ||psi|| <= sum |c_n|
    # c = L - (pi/4) sum_j a_j c_{j,1}
    L = log_q(R)
    s_lo = s_hi = F(0)
    for sl in SLITS:
        c1 = co[sl.cl][1]
        lo_, hi_ = iv_mul(sl.a_q, (c1, c1))
        s_lo, s_hi = s_lo + lo_, s_hi + hi_
    t_lo, t_hi = iv_mul((PI[0] / 4, PI[1] / 4), (s_lo, s_hi))
    c_iv = (L[0] - t_hi, L[1] - t_lo)
    log("  sha256 %s; c = L - (pi/4) sum a_j c_j1 in [%.16f, %.16f]" % (sha[:16] + "...", float(c_iv[0]), float(c_iv[1])))
    deltas = {}
    for i in REPS:
        ti = time.time()
        si = SLITS[i]
        cl = si.cl
        # Step 1: for each candidate b, rigorous kernel bounds for every j (branch and bound); keep the b with the
        # smallest M (a candidate whose branch and bound does not close is skipped)
        best = None
        nbox = 0
        for b in B_CANDS:
            hw = math.sqrt(1 + b * b)
            r_q = F(b + math.sqrt(1 + b * b)) * (1 - F(1, 10 ** 9))   # slightly below b + sqrt(1 + b^2)
            Qp, Qq = up(hw * (1 + 1e-12)), up(b * (1 + 1e-12))
            need((r_q + 1 / r_q) / 2 <= F(Qp) and (r_q - 1 / r_q) / 2 <= F(Qq), "ellipse not inside the rectangle")
            A = {}
            for j in range(24):
                res = bnb(i, j, ((-Qp, Qp), (-Qq, Qq), (-1.0, 1.0)), "x", max_boxes=40000)
                if res is None:
                    break
                lo2, hi2, nb_ = res
                nbox += nb_
                A[j] = max(abs(0.5 * dn4(math.log(lo2))), abs(0.5 * up4(math.log(hi2)))) * (1 + 1e-12) + 1e-15
            if len(A) < 24:
                continue
            Bq = F(si.a_q[1]) * (r_q + 1 / r_q) / 2
            for n in range(1, 64, 2):
                Bq += abs(co[cl][n]) * PI[1] / (2 * n) * (r_q ** n + r_q ** -n)
            Bq += PI[1] * sum(norm[SLITS[j].cl] * F(A[j]) for j in range(24))
            Mf = 1 + math.log(8 * float(Bq) / ((float(r_q) - 1) * float(tol_x))) / math.log(float(r_q))
            if best is None or Mf < best[0]:
                best = (Mf, b, r_q, Bq)
        need(best is not None, "no ellipse closed for slit %d" % i)
        _, b, r_q, Bq = best
        r = float(r_q)
        M = int(math.ceil(best[0])) + 2
        M += M % 2
        rem_x = 8 * Bq / (r_q - 1) / r_q ** (M - 1)
        need(rem_x <= tol_x, "interpolation remainder")
        rem_x = up(float(ru(rem_x)))
        # Step 2: per j, the y-ellipse with the smallest N among the candidates whose branch and bound closes
        N, qerr = {}, {}
        nbox2 = 0
        for j in range(24):
            bestj = None
            for bp in B_CANDS:
                hwp = math.sqrt(1 + bp * bp)
                rp = F(bp + math.sqrt(1 + bp * bp)) * (1 - F(1, 10 ** 9))
                Qp2, Qq2 = up(hwp * (1 + 1e-12)), up(bp * (1 + 1e-12))
                need((rp + 1 / rp) / 2 <= F(Qp2) and (rp - 1 / rp) / 2 <= F(Qq2), "y ellipse not inside the rectangle")
                res = bnb(i, j, ((-Qp2, Qp2), (-Qq2, Qq2), (-1.0, 0.0)), "y", max_boxes=40000)
                if res is None:
                    continue
                lo2, hi2, nb_ = res
                nbox2 += nb_
                Ap = F(max(abs(0.5 * dn4(math.log(lo2))), abs(0.5 * up4(math.log(hi2)))) * (1 + 1e-12) + 1e-15)
                pre = 8 * PI[1] * Ap * norm[SLITS[j].cl] / (rp - 1)
                nf = (64 + math.log(float(pre) / float(tol_q)) / math.log(float(rp))) / 2
                if bestj is None or nf < bestj[0]:
                    bestj = (nf, rp, pre)
            need(bestj is not None, "no y ellipse closed for the pair %d, %d" % (i, j))
            nf, rp, pre = bestj
            nn = max(40, int(math.ceil(nf)) + 1)
            qe = pre / rp ** (2 * nn - 1 - 63)
            need(qe <= tol_q, "quadrature remainder")
            N[j] = nn
            qerr[j] = up(float(ru(qe)))
        # Step 3: the nodes, in C
        lines = ["%d %d" % (M, 24), "%s %s" % (PI_I.lo.hex(), PI_I.hi.hex()), "%d" % i]
        for c_ in range(5):
            lines.append(" ".join(float(v).hex() for v in co[c_]))
            need(all(F(float(v)) == v for v in co[c_]), "a coefficient is not a binary64 number")
        for j in range(24):
            sj = SLITS[j]
            lines.append("%d %s %s %s %s %s %s %s %s %d" % (sj.cl, sj.A.lo.hex(), sj.A.hi.hex(), sj.B1.lo.hex(),
                         sj.B1.hi.hex(), sj.B2.lo.hex(), sj.B2.hi.hex(), sj.TH.lo.hex(), sj.TH.hi.hex(), N[j]))
        out = subprocess.run([os.path.join(HERE, "bin", "vs_nodes")], input="\n".join(lines) + "\n",
                             capture_output=True, text=True, check=True).stdout.split("\n")
        vals = {}
        for ln in out:
            if ln.strip():
                k, lo_s, hi_s = ln.split()
                vals[int(k)] = (float.fromhex(lo_s), float.fromhex(hi_s))
        need(len(vals) == M // 2, "node count")
        qsum = F(0)
        for v in qerr.values():
            qsum += F(v)
        # node values with the quadrature remainders; oddness for the positive nodes
        f_lo, f_hi = {}, {}
        for k, (lo_, hi_) in vals.items():
            f_lo[k], f_hi[k] = F(lo_) - qsum, F(hi_) + qsum
            f_lo[M - 1 - k], f_hi[M - 1 - k] = -f_hi[k], -f_lo[k]
        maxnode = max(max(abs(f_lo[k]), abs(f_hi[k])) for k in range(M))
        halfw = max((f_hi[k] - f_lo[k]) / 2 for k in range(M))
        mids = [(f_lo[k] + f_hi[k]) / 2 for k in range(M)]
        # Step 4: gamma~_q = (2/M) sum_k mid_k cos(q (2k+1) pi/(2M)) (q >= 1), in binary64 intervals from the exact
        # midpoints (each rounded outward first); oddness: only odd q are non-zero
        cosv = []
        for m in range(4 * M):
            ang = I(m) * PI_I
            ang = ang.div(I(float(2 * M)))
            mid_ = 0.5 * (ang.lo + ang.hi)
            w_ = up(max(mid_ - ang.lo, ang.hi - mid_))
            cv = math.cos(mid_)
            cosv.append(I(dn(dn4(cv) - w_), up(up4(cv) + w_)))
        midI = [I.q(mv) for mv in mids]
        # f odd and the nodes symmetric: gamma_q = 0 for even q (paired terms cancel exactly, the midpoints being exactly
        # odd), and for odd q gamma_q = (4/M) sum over the negative nodes k >= M/2
        sum_abs = 0.0
        for q_ in range(1, M, 2):
            acc_lo = acc_hi = 0.0
            for k in range(M // 2, M):
                p_ = midI[k] * cosv[(q_ * (2 * k + 1)) % (4 * M)]
                acc_lo, acc_hi = dn(acc_lo + p_.lo), up(acc_hi + p_.hi)
            g = (I(acc_lo, acc_hi) * I(4.0)).div(I(float(M)))
            sum_abs = up(sum_abs + g.mag())
        need(all(mids[k] == -mids[M - 1 - k] for k in range(M)), "midpoints not odd")
        corr = math.sqrt(2 * M) * float(halfw) * (1 + 1e-12)
        delta_i = F(up(up(sum_abs + corr))) + F(rem_x) * (1 + F(1, 10 ** 9))
        deltas[i] = delta_i
        log("  slit %2d (class %d): b = %.3f (r = %.5f), B = %.2f, M = %d, N = %d..%d; max |f| at nodes <= %.3e, "
            "node half-width <= %.1e, sum |gamma| <= %.3e, interpolation remainder <= %.1e, quadrature <= %.1e; "
            "delta_%d <= %.4e  [%d + %d boxes, %.0f s]"
            % (i, cl, b, r, float(Bq), M, min(N.values()), max(N.values()), float(maxnode), float(halfw), sum_abs,
               float(rem_x), float(qsum), i, float(delta_i), nbox, nbox2, time.time() - ti))
    delta = max(deltas.values())
    e = exp_q(c_iv[0] - delta)
    log("  delta <= %.4e; log rho >= c - delta >= %.13f; rho(D) >= %.12f" % (float(delta), float(c_iv[0] - delta),
                                                                          float(e[0])))
    return e[0], delta, sha, time.time() - t0


def table_check(json_path, tex_path):
    """the 32 odd rows x 5 classes of the paper's appendix table (read as text) against the data file"""
    import re
    co = json.load(open(json_path))["coefficients"]
    tex = open(tex_path, encoding="utf-8").read()
    tab = tex[tex.index("\\label{tab:coef}") - 12000:tex.index("\\label{tab:coef}")]
    rows = re.findall(r"^(\d+) & (.*?)\\\\$", tab, re.M)
    bad = 0
    for n, rest in rows:
        vals = [float(x.strip().strip("$").replace("\\mathrm{e}{-}", "e-")) for x in rest.split("&")]
        bad += sum(co[str(c)][int(n)] != vals[c] for c in range(5))      # the printed decimals round to the file's doubles
    return len(rows), bad


def main():
    path = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    controls = "--controls" in sys.argv
    print("57c [vS2026], the conformal radius: rho(D) >= 1.523726208329 (Proposition 4.1)")
    tex = [a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--tex=")]
    if tex:
        nrow, bad = table_check(path, tex[0])
        print("  the data file against the paper's appendix table: %d rows x 5 classes; the printed decimals round to the file's binary64 values, which we use exactly: %d differences" % (nrow, bad))
        need(nrow == 32 and bad == 0, "the data file differs from the appendix table")
    rho, delta, sha, tt = run(path)
    ok_sha = sha == SHA_DENSITY
    print("  density file sha256 %s (the paper's appendix prints f2d1d4e6...843b; the coefficients are compared with "
          "its table separately)" % ("as pinned" if ok_sha else "DIFFERENT"))
    ok = rho >= CLAIM and ok_sha
    print("  rho(D) >= 1.523726208329: %s" % (rho >= CLAIM))
    lam = 1 + F(1, 10 ** 9)
    ratio_ok = lam / rho < F("0.65628589")
    print("  with R(D) <= 1 + 10^-9: U <= R(D)/rho(D) <= %.12f < 0.65628589: %s" % (float(lam / rho), ratio_ok))
    ok &= ratio_ok
    print("  run time %.0f s" % tt)
    if controls:
        print("Red controls (each must be refused):")
        # c_{0,1} lowered by 1% (class 0, the largest coefficient): c drops, the comparison function no longer vanishes
        rho2, d2, _, _ = run(path, coef_override={(0, 1): F(0.2658439564499184 * 0.99)}, log=lambda s: None)
        c1 = not (rho2 >= CLAIM)
        print("  1. c_{0,1} (class 0) lowered by 1%%: delta <= %.2e, rho >= %.9f, so 1.523726208329 is not certified: %s"
              % (float(d2), float(rho2), "refused" if c1 else "ACCEPTED (bad)"))
        ok &= c1
        rho3, d3, _, _ = run(path, coef_override={(2, 5): F(0.000755213124641526 * 1.01)}, log=lambda s: None)
        c2 = not (rho3 >= CLAIM)
        print("  2. c_{2,5} (class 2) raised by 1%%: delta <= %.2e, rho >= %.9f: %s"
              % (float(d3), float(rho3), "refused" if c2 else "ACCEPTED (bad)"))
        ok &= c2
        try:
            run(path, coef_override={(1, 2): F(1, 2 ** 30)}, log=lambda s: None)
            c3 = False
        except AssertionError:
            c3 = True
        print("  3. an even-indexed coefficient made non-zero (psi no longer odd, so U is not odd under the reflection): %s"
              % ("refused" if c3 else "ACCEPTED (bad)"))
        ok &= c3
        c4 = not (rho >= F("1.5237262489"))
        print("  4. the bound 1.5237262489, above the floating-point value 1.523726248805 of the conformal radius: "
              "%s" % ("not certified, as it must not be" if c4 else "CERTIFIED (bad)"))
        ok &= c4
    print("RESULT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
