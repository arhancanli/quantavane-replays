#!/usr/bin/env python3
"""Independent replay of the 22b row [Claude2026], C_22b <= 7.2012719* (PR #267): for every odd n >= 59 an explicit
C^{1,1} curve of type T(2, n) with thickness >= 1 and length <= 7.2012719 n + 19.649.

Written from the paper alone (YouhuaLi/tight-alternating-knot-constant @ debf7c92: paper/main.tex and numbers.tex read
as text, design.json read as data). The one constant that is only in code, the seven decimal coefficients of the even
polynomial q in the cap phase, was taken from code/capaxis.py at that commit by the coordinator's ruling of 12 Nov
2026: Python's ast parser, the single dict literal QC, the string argument of each of its seven Fr(...) calls; nothing
executed and no other part of any program read (take_q.py repeats that extraction). No program of the repository
was read or run.

What is checked, for the exact decimal parameters (rationals) and our own reading of Section 4 (the cap phase of
cap 1 uses (E1, F1), cap 2 uses (E2, F2); every phase offset follows from continuity along the lap):
  exact   the junction identities (X, Y, q, beta, S5 and the phase profiles at tau = 0, 1), the closing phase, c,
          the curvatures of the helix (H4), (H3), L_s(59), and the sign of the axis curvature (Sturm sequences);
  (H1-2)  the helix lemma by bisection with the mean-value form in rigorous interval arithmetic (ia22.py);
  (V2)    in the top zone (both entry phases) and the bottom zone: the curvature of the free pieces, and every pair
          of points with a point on a free piece either has an upper bound < pi on its arc separation or a certified
          squared distance > 4 (element pairs, mean-value form, bisection);
  (V3)    the separation constants: cap 1 above L_s - delta_t, cap 2 below y_cap2, the swap radii, L_s(n) for n >= 59;
  (V4)    the axis radius of curvature and the strand radii (the reach condition of the knot-type lemma);
  (V5)    rigorous upper bounds for the strand lengths of the free pieces, hence C <= 19.649.
Red controls (--controls): a1 at the helix optimum, W = 2, E1 = -19.5, m_o = 0, the Kim-Oh-Huh helix, a perturbed q.

Usage: python3 -I c22b_check.py [--controls] [--sources DIR]     (Python 3.9+, standard library only)
"""
import json
import os
import sys
import time
from fractions import Fraction as Fr
from math import inf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ia22 import Iv, Jet, ivq, poly_jet, sincos, PI, PI_LO, PI_HI, dn, up   # noqa: E402

# ---------------------------------------------------------------------------------------------------------------
# data
DESIGN = {"a1": "1.934062", "a2": "0.075486", "b": "0.322723", "R": "6.06", "E1": "-20.5", "F1": "0",
          "E2": "-20.53", "F2": "-0.42", "mo": "-0.015", "mi": "0.54", "Sw": "6.76", "Ew": "-10.41", "Fw": "0",
          "mw": "0.224", "W": "5.5"}                       # design.json at debf7c92 (sha256 a900eed1...a25a569)
Q_STRINGS = {0: "2.516", 2: "-2.4075", 4: "-1.0513", 6: "0.1583", 8: "9.5169", 10: "-6.5804", 12: "1.634"}
#   ^ capaxis.py at debf7c92 (sha256 559e5a1c...5d90), the seven Fr('...') arguments of QC, taken with ast only
X_COEF = {1: Fr(2517, 1000), 3: Fr(-3511, 1000), 5: Fr(1653, 500), 7: Fr(-1647, 1000), 9: Fr(67, 200)}  # main.tex
Y_COEF = {0: Fr(51, 125), 2: Fr(-843, 1000), 4: Fr(493, 500), 6: Fr(-23, 40), 8: Fr(3, 125)}           # main.tex
H_OVER_R = Fr("3.786")
PRINTED = dict(c_up=Fr("7.2012719"), c_lo=Fr("7.2012718"), C=Fr("19.649"), n0=59)


def pmul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i + j] = out.get(i + j, 0) + a * b
    return {k: v for k, v in out.items() if v}


def padd(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}


def pscale(p, c):
    return {k: v * c for k, v in p.items() if v * c}


def pder(p):
    return {k - 1: v * k for k, v in p.items() if k}


def pint(p):
    return {k + 1: v / (k + 1) for k, v in p.items()}


def peval(p, x):
    return sum(v * x ** k for k, v in p.items())


BETA = {2: Fr(16), 3: Fr(-32), 4: Fr(16)}            # 16 t^2 (1 - t)^2
S5 = {3: Fr(10), 4: Fr(-15), 5: Fr(6)}
H3 = {2: Fr(3), 3: Fr(-2)}                            # 3 t^2 - 2 t^3
G4 = {2: Fr(1), 3: Fr(-2), 4: Fr(1)}                  # t^2 (1 - t)^2


class Design:
    """exact parameters and derived constants; overrides give the red controls"""

    def __init__(self, over=None, qover=None):
        p = {k: Fr(v) for k, v in DESIGN.items()}
        p.update({k: Fr(v) for k, v in (over or {}).items()})
        self.p = p
        for k, v in p.items():
            setattr(self, k, v)
        self.q = {k: Fr(v) for k, v in Q_STRINGS.items()}
        self.q.update({k: Fr(v) for k, v in (qover or {}).items()})
        self.H = H_OVER_R * self.R
        lam = pint(self.q)                                # Lambda(v) = int_{-1}^v q
        self.lam = padd(lam, {0: -peval(lam, Fr(-1))})
        self.T0 = self.R * peval(self.lam, Fr(1)) / self.b
        self.Thf = 2 * self.T0 + self.E1 + self.E2 + self.Ew
        self.Xd, self.Yd = pder(X_COEF), pder(Y_COEF)
        # interval constants
        I = ivq
        self.iv = {k: I(v) for k, v in p.items()}
        self.iv.update(H=I(self.H), Rb=I(self.R / self.b), Swb=I(self.Sw / self.b), Wb=I(self.W / self.b))
        sa = (I(self.a1 * self.a1 + self.b * self.b)).sqrt()
        sb = (I(self.a2 * self.a2 + self.b * self.b)).sqrt()
        self.ell = (sa + sb) / I(self.b)

    def theta1(self, case):
        """phase at the entry of cap 1, mod 2 pi: (n pi - Theta_f)/2 with n = 1 or 3 (mod 4)"""
        return PI * (0.5 if case == 1 else 1.5) - ivq(self.Thf / 2)

    def bottom_phases(self):
        th5 = PI                                              # end of cap 2: n pi, n odd
        th4 = th5 - ivq(self.T0 + self.E2)
        th3 = th4 - ivq(self.Sw / self.b + self.Ew)
        return th3, th4, th5


# ---------------------------------------------------------------------------------------------------------------
# pieces: each returns jets (order K) of P_A and P_B over the base interval of tau, and the axis point
class Piece:
    def __init__(self, d, kind, th0, free, **kw):
        self.d, self.kind, self.th0, self.free = d, kind, th0, free
        self.__dict__.update(kw)

    def jets(self, base, K):
        d, iv = self.d, self.d.iv
        t = Jet.var(base, K)
        if self.kind == "col":     # axis (x0, y0 + dir L t), T = (0, dir), N = (-dir, 0)
            th = t * ivq(self.L / d.b) + self.th0
            s, c = th.sincos()
            ux = c * float(-self.dir)
            ax = Jet.const(ivq(self.x0), K)
            ay = t * ivq(self.dir * self.L) + ivq(self.y0)
            rA, rB = ivq(self.rA), ivq(self.rB)
            PA = (ax + ux * rA, ay, s * rA)
            PB = (ax - ux * rB, ay, -(s * rB))
            return PA, PB
        if self.kind == "swap":    # axis (2R, y0 + Sw (1 - t)), T = (0, -1), N = (1, 0)
            th = t * iv["Swb"] + poly_jet(H3, t) * iv["Ew"] + poly_jet(G4, t) * iv["Fw"] + self.th0
            s, c = th.sincos()
            be, s5 = poly_jet(BETA, t), poly_jet(S5, t)
            rA = s5 * ivq(d.a2 - d.a1) + be * iv["mw"] + iv["a1"]
            rB = s5 * ivq(d.a1 - d.a2) + be * iv["mw"] + iv["a2"]
            ax = Jet.const(ivq(2 * d.R), K)
            ay = t * ivq(-d.Sw) + ivq(self.y0 + d.Sw)
            PA = (ax + c * rA, ay, s * rA)
            PB = (ax - c * rB, ay, -(s * rB))
            return PA, PB
        # caps: sg = +1 for cap 1, -1 for cap 2 (cap 2 is cap 1 turned by pi)
        sg = 1 if self.kind == "cap1" else -1
        v = t * 2.0 - 1.0
        X, Y = poly_jet(X_COEF, v), poly_jet(Y_COEF, v)
        Xd, Yd = poly_jet(d.Xd, v), poly_jet(d.Yd, v)
        ax = X * ivq(sg * d.R) + ivq(d.R)
        ay = Y * ivq(sg * d.H) + ivq(self.y0)
        dx, dy = Xd * ivq(2 * sg * d.R), Yd * ivq(2 * sg * d.H)    # d alpha / d tau
        sp = (dx * dx + dy * dy).sqrt()
        isp = sp.inv()
        Tx, Ty = dx * isp, dy * isp
        E, F = (iv["E1"], iv["F1"]) if sg == 1 else (iv["E2"], iv["F2"])
        th = poly_jet(d.lam, v) * iv["Rb"] + poly_jet(H3, t) * E + poly_jet(G4, t) * F + self.th0
        s, c = th.sincos()
        be = poly_jet(BETA, t)
        if sg == 1:
            rA, rB = (be * iv["mo"] + 1.0) * iv["a1"], (be * iv["mi"] + 1.0) * iv["a2"]
        else:
            rA, rB = (be * iv["mi"] + 1.0) * iv["a2"], (be * iv["mo"] + 1.0) * iv["a1"]
        ux, uy = c * (-Ty), c * Tx                     # N = (-T_y, T_x, 0)
        PA = (ax + ux * rA, ay + uy * rA, s * rA)
        PB = (ax - ux * rB, ay - uy * rB, -(s * rB))
        return PA, PB


def top_zone(d, case):
    th1 = d.theta1(case)
    th2 = th1 + ivq(d.T0 + d.E1)
    W = d.W
    return [Piece(d, "col", th1 - ivq(W / d.b), False, x0=Fr(0), y0=-W, dir=1, L=W, rA=d.a1, rB=d.a2),
            Piece(d, "cap1", th1, True, y0=Fr(0)),
            Piece(d, "col", th2, False, x0=2 * d.R, y0=Fr(0), dir=-1, L=W, rA=d.a1, rB=d.a2)]


def bottom_zone(d):
    th3, th4, th5 = d.bottom_phases()
    W = d.W
    return [Piece(d, "col", th3 - ivq(W / d.b), False, x0=2 * d.R, y0=d.Sw + W, dir=-1, L=W, rA=d.a1, rB=d.a2),
            Piece(d, "swap", th3, True, y0=Fr(0)),
            Piece(d, "cap2", th4, True, y0=Fr(0)),
            Piece(d, "col", th5, False, x0=Fr(0), y0=Fr(0), dir=1, L=W, rA=d.a2, rB=d.a1)]


# ---------------------------------------------------------------------------------------------------------------
# helpers on 3-vectors of intervals
def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    s = dot(a, a)
    return Iv(max(s.lo, 0.0), s.hi).sqrt()


class ElementCache:
    def __init__(self):
        self.memo = {}

    def get(self, piece, t0, t1):
        key = (id(piece), t0, t1)
        r = self.memo.get(key)
        if r is None:
            m = 0.5 * (t0 + t1)
            PcA, PcB = piece.jets(Iv(m), 0)
            JA, JB = piece.jets(Iv(t0, t1), 1)
            h = 0.5 * (t1 - t0)
            r = []
            for Pc, J in ((PcA, JA), (PcB, JB)):
                c = tuple(j.c[0] for j in Pc)
                dP = tuple(j.c[1] for j in J)
                box = tuple(c[i] + dP[i] * Iv(-h, h) for i in range(3))
                r.append((c, dP, h, box))
            self.memo[key] = r
        return r


def pair_lower_bound(e, f):
    """lower bound of |P(s) - P(t)|^2 over s in e, t in f (mean-value form, as in the paper's Section 5)"""
    ce, de, he, be = e
    cf, df, hf, bf = f
    gap = 0.0
    for i in range(3):
        g = max(be[i].lo - bf[i].hi, bf[i].lo - be[i].hi, 0.0)
        gap += g * g
    if dn(gap) > 4.0:
        return dn(gap), True
    D0 = dot(vsub(ce, cf), vsub(ce, cf))
    delta = vsub(be, bf)
    ge, gf = dot(delta, de), dot(delta, df)
    D = D0 + ge * Iv(-2 * he, 2 * he) - gf * Iv(-2 * hf, 2 * hf)
    return D.lo, False


def point_D_upper(piece_e, se, te, piece_f, sf, tf):
    """upper bound of |P(s) - P(t)|^2 at the centres (a certified witness when below 4)"""
    PeA, PeB = piece_e.jets(Iv(te), 0)
    PfA, PfB = piece_f.jets(Iv(tf), 0)
    pe = tuple(j.c[0] for j in (PeA if se == 0 else PeB))
    pf = tuple(j.c[0] for j in (PfA if sf == 0 else PfB))
    return dot(vsub(pe, pf), vsub(pe, pf)).hi


# ---------------------------------------------------------------------------------------------------------------
# (V2): zone pairs
def check_zone(pieces, label, log, n_elem=None, max_depth=16, stop_on_fail=True):
    t0 = time.time()
    cache = ElementCache()
    n_elem = n_elem or {}
    # elements: dyadic tau grids; arc order along each strand: pieces in order, tau increasing
    elems = []          # (strand, piece index, t0, t1, top index on its strand)
    for strand in (0, 1):
        k = 0
        for pi_, pc in enumerate(pieces):
            n = n_elem.get((pc.kind, pc.free), 2048 if pc.kind.startswith("cap") else 512)
            for j in range(n):
                elems.append((strand, pi_, j / n, (j + 1) / n, k))
                k += 1
    # arc-length upper bounds and prefix sums per strand
    lens = {0: [], 1: []}
    lens_lo = {0: [], 1: []}
    data = []
    for (strand, pi_, a, b, k) in elems:
        r = cache.get(pieces[pi_], a, b)[strand]
        data.append(r)
        n1 = norm(r[1])
        lens[strand].append(up(up(b - a) * n1.hi))
        lens_lo[strand].append(dn(dn(b - a) * n1.lo))
    cum_up, cum_lo = {}, {}
    for s in (0, 1):
        cu, cl = [0.0], [0.0]
        for x in lens[s]:
            cu.append(up(cu[-1] + x))
            cl.append(dn(cl[-1] + x))
        cum_up[s], cum_lo[s] = cu, cl
    low_up, low_lo = {}, {}      # prefix sums of the lower bounds (rounded up and down)
    for s in (0, 1):
        cu, cl = [0.0], [0.0]
        for x in lens_lo[s]:
            cu.append(up(cu[-1] + x))
            cl.append(dn(cl[-1] + x))
        low_up[s], low_lo[s] = cu, cl

    def sep_lower(si, ki, sj, kj):
        if si != sj:
            return inf                    # different laps of the curve: far apart along it
        lo_k, hi_k = min(ki, kj), max(ki, kj)
        if hi_k - lo_k < 2:
            return 0.0
        return dn(low_lo[si][hi_k] - low_up[si][lo_k + 1])
    # candidate pairs by a grid on the element boxes expanded by 1
    cell = 2.0
    grid = {}
    boxes = []
    for idx, r in enumerate(data):
        bx = r[3]
        lo = tuple(int((bx[i].lo - 1.001) // cell) for i in range(3))
        hi = tuple(int((bx[i].hi + 1.001) // cell) for i in range(3))
        boxes.append((lo, hi))
        for gx in range(lo[0], hi[0] + 1):
            for gy in range(lo[1], hi[1] + 1):
                for gz in range(lo[2], hi[2] + 1):
                    grid.setdefault((gx, gy, gz), []).append(idx)
    stats = dict(pairs=0, discarded_arc=0, helix_helix=0, accepted=0, bisections=0, min_lb=inf, depth=0)
    fails = []
    for key, lst in grid.items():
        for ii in range(len(lst)):
            i = lst[ii]
            for jj in range(ii + 1, len(lst)):
                j = lst[jj]
                # test each pair once: in the cell holding the lower corner of the two expanded boxes' overlap
                (li, hi_), (lj, hj) = boxes[i], boxes[j]
                corner = (max(li[0], lj[0]), max(li[1], lj[1]), max(li[2], lj[2]))
                if corner != key or any(max(li[q], lj[q]) > min(hi_[q], hj[q]) for q in range(3)):
                    continue
                si, pi_i, ai, bi, ki = elems[i]
                sj, pi_j, aj, bj, kj = elems[j]
                if not pieces[pi_i].free and not pieces[pi_j].free:
                    stats["helix_helix"] += 1
                    continue
                if si == sj:
                    lo_k, hi_k = min(ki, kj), max(ki, kj)
                    if up(cum_up[si][hi_k + 1] - cum_lo[si][lo_k]) < 3.14159:    # < pi
                        stats["discarded_arc"] += 1
                        continue
                stats["pairs"] += 1
                # branch and bound on this pair
                stack = [(ai, bi, aj, bj, 0)]
                while stack:
                    a1_, b1_, a2_, b2_, dep = stack.pop()
                    e = cache.get(pieces[pi_i], a1_, b1_)[si] if dep else data[i]
                    f = cache.get(pieces[pi_j], a2_, b2_)[sj] if dep else data[j]
                    lb, by_box = pair_lower_bound(e, f)
                    if lb > 4.0:
                        stats["accepted"] += 1
                        if not by_box:
                            stats["min_lb"] = min(stats["min_lb"], lb)
                        stats["depth"] = max(stats["depth"], dep)
                        continue
                    if dep >= max_depth:
                        ub = point_D_upper(pieces[pi_i], si, 0.5 * (a1_ + b1_), pieces[pi_j], sj, 0.5 * (a2_ + b2_))
                        fails.append((lb, ub, si, pi_i, a1_, b1_, sj, pi_j, a2_, b2_))
                        if stop_on_fail:
                            stack = []
                            break
                        continue
                    # a certified violation ends the search at once
                    ub = point_D_upper(pieces[pi_i], si, 0.5 * (a1_ + b1_), pieces[pi_j], sj, 0.5 * (a2_ + b2_))
                    if ub < 4.0 and sep_lower(si, ki, sj, kj) >= 3.1416:
                        fails.append((lb, ub, si, pi_i, a1_, b1_, sj, pi_j, a2_, b2_))
                        stack = []
                        break
                    stats["bisections"] += 1
                    m1, m2 = 0.5 * (a1_ + b1_), 0.5 * (a2_ + b2_)
                    w1 = (b1_ - a1_) * norm(e[1]).hi
                    w2 = (b2_ - a2_) * norm(f[1]).hi
                    if w1 >= 0.5 * w2 and w2 >= 0.5 * w1:
                        stack += [(a1_, m1, a2_, m2, dep + 1), (a1_, m1, m2, b2_, dep + 1),
                                  (m1, b1_, a2_, m2, dep + 1), (m1, b1_, m2, b2_, dep + 1)]
                    elif w1 > w2:
                        stack += [(a1_, m1, a2_, b2_, dep + 1), (m1, b1_, a2_, b2_, dep + 1)]
                    else:
                        stack += [(a1_, b1_, a2_, m2, dep + 1), (a1_, b1_, m2, b2_, dep + 1)]
                if fails and stop_on_fail:
                    break
            if fails and stop_on_fail:
                break
        if fails and stop_on_fail:
            break
    ok = not fails
    log(f"  (V2) {label}: {len(elems)} elements, {stats['pairs']:,} element pairs tested (plus "
        f"{stats['discarded_arc']:,} discarded by arc separation < pi and {stats['helix_helix']:,} helix-helix pairs "
        f"left to (H1)-(H3)), {stats['bisections']:,} bisections, depth <= {stats['depth']}, "
        f"{len(cache.memo):,} element evaluations; smallest certified lower bound of |P(s) - P(t)|^2 "
        f"(mean-value form): {stats['min_lb']:.10f}; {'all pairs certified' if ok else 'FAILED'} "
        f"({time.time() - t0:.0f} s)")
    for fl in fails[:3]:
        lb, ub, si, pi_i, a1_, b1_, sj, pi_j, a2_, b2_ = fl
        log(f"    failing pair: strand {'AB'[si]} on {pieces[pi_i].kind} tau in [{a1_:.7f}, {b1_:.7f}], strand "
            f"{'AB'[sj]} on {pieces[pi_j].kind} tau in [{a2_:.7f}, {b2_:.7f}]: lower bound {lb:.8f}, at the centres "
            f"|P(s) - P(t)|^2 <= {ub:.8f}{' < 4 (a certified violation)' if ub < 4 else ''}")
    return ok, stats, fails


# ---------------------------------------------------------------------------------------------------------------
# curvature, heights, lengths over free pieces (adaptive bisection on tau)
def curvature_bound(piece, target, log_label, n0=2048, max_depth=12):
    """certified upper bound of the curvature of both strands over the piece (< target required)"""
    worst = [0.0, 0.0]
    stack = [(j / n0, (j + 1) / n0, 0) for j in range(n0)]
    evals = 0
    while stack:
        a, b, dep = stack.pop()
        JA, JB = piece.jets(Iv(a, b), 2)
        evals += 1
        bad = False
        vals = []
        for s, J in enumerate((JA, JB)):
            d1 = tuple(j.c[1] for j in J)
            d2 = tuple(j.c[2] * 2.0 for j in J)
            num = norm(cross(d1, d2)).hi
            den = norm(d1).lo
            k = up(num / dn(den * dn(den * den))) if den > 0 else inf
            vals.append(k)
        if max(vals) >= target and dep < max_depth:
            m = 0.5 * (a + b)
            stack += [(a, m, dep + 1), (m, b, dep + 1)]
            continue
        for s in (0, 1):
            worst[s] = max(worst[s], vals[s])
    return worst, evals


def height_bound(piece, coord, sign, n0=1024, target=None, max_depth=14):
    """certified bound of sign * P_coord over both strands: max of sign * P (upper bound)"""
    best = -inf
    stack = [(j / n0, (j + 1) / n0, 0) for j in range(n0)]
    while stack:
        a, b, dep = stack.pop()
        m = 0.5 * (a + b)
        PcA, PcB = piece.jets(Iv(m), 0)
        JA, JB = piece.jets(Iv(a, b), 1)
        h = 0.5 * (b - a)
        hi = -inf
        for Pc, J in ((PcA, JA), (PcB, JB)):
            val = Pc[coord].c[0] + J[coord].c[1] * Iv(-h, h)
            hi = max(hi, (val * float(sign)).hi)
        if target is not None and hi > target and dep < max_depth:
            stack += [(a, m, dep + 1), (m, b, dep + 1)]
            continue
        best = max(best, hi)
    return best


def length_bound(piece, strand_set=(0, 1), tol=2.5e-4, n0=4096, max_depth=8):
    """rigorous bounds for the strand lengths over the piece: int |P'| <= 2h |P'(m)| + (2h)^3/24 sup (|P'|)''
    with (|P'|)'' <= (|P''|^2 + |P'| |P'''|)/|P'|; refined until each element's error term is below tol * 2h"""
    hi = [Iv(0.0), Iv(0.0)]
    lo = [Iv(0.0), Iv(0.0)]
    err_total = [0.0, 0.0]
    evals = 0
    stack = [(j / n0, (j + 1) / n0, 0) for j in range(n0)]
    while stack:
        a, b, dep = stack.pop()
        m, h = 0.5 * (a + b), 0.5 * (b - a)
        JmA, JmB = piece.jets(Iv(m), 1)
        JA, JB = piece.jets(Iv(a, b), 3)
        evals += 1
        res = []
        for s, (Jm, J) in enumerate(((JmA, JA), (JmB, JB))):
            sp_m = norm(tuple(j.c[1] for j in Jm))
            d1 = tuple(j.c[1] for j in J)
            d2 = tuple(j.c[2] * 2.0 for j in J)
            d3 = tuple(j.c[3] * 6.0 for j in J)
            n1 = norm(d1)
            g2 = (dot(d2, d2) + n1 * norm(d3)) / Iv(n1.lo)
            e = up(up(8.0 * h * h * h / 24.0) * g2.hi)
            res.append((sp_m, e))
        if max(r[1] for r in res) > tol * 2 * h and dep < max_depth:
            stack += [(a, m, dep + 1), (m, b, dep + 1)]
            continue
        for s in strand_set:
            sp_m, e = res[s]
            hi[s] = hi[s] + Iv(up(2 * h * sp_m.hi + e))
            lo[s] = lo[s] + Iv(dn(2 * h * sp_m.lo - e))
            err_total[s] = up(err_total[s] + e)
    return [Iv(lo[s].lo, hi[s].hi) for s in (0, 1)], err_total, evals


# ---------------------------------------------------------------------------------------------------------------
# exact checks
def sturm_count(p, a, b):
    """number of distinct real roots of the polynomial p (dict) in (a, b], by a Sturm sequence over Q"""
    deg = max(k for k, v in p.items() if v)
    P0 = [p.get(k, Fr(0)) for k in range(deg + 1)]

    def trim(c):
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        return c

    def pdiv_rem(n, d):
        n = list(n)
        while len(n) >= len(d) and any(n):
            f = n[-1] / d[-1]
            sh = len(n) - len(d)
            for i in range(len(d)):
                n[sh + i] -= f * d[i]
            n.pop()
        return trim(n or [Fr(0)])
    P1 = trim([P0[k] * k for k in range(1, len(P0))])
    seq = [trim(P0), P1]
    while len(seq[-1]) > 1 or seq[-1][0] != 0:
        r = pdiv_rem(seq[-2], seq[-1])
        if all(x == 0 for x in r):
            break
        seq.append([-x for x in r])

    def var(x):
        vals = [sum(c * x ** i for i, c in enumerate(s)) for s in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))
    return var(a) - var(b)


def exact_checks(d, log):
    ok = True
    m1, p1 = Fr(-1), Fr(1)
    Xd2, Yd2 = pder(d.Xd), pder(d.Yd)
    conds = [("X(+-1) = +-1", peval(X_COEF, p1) == 1 and peval(X_COEF, m1) == -1),
             ("X'(+-1) = X''(+-1) = 0", all(peval(q, x) == 0 for q in (d.Xd, Xd2) for x in (m1, p1))),
             ("Y(+-1) = 0, Y'(+-1) = -+1", peval(Y_COEF, p1) == 0 and peval(Y_COEF, m1) == 0
              and peval(d.Yd, p1) == -1 and peval(d.Yd, m1) == 1),
             ("q even of degree 12, q(+-1) = H/R = 3.786", max(d.q) == 12 and all(k % 2 == 0 for k in d.q)
              and peval(d.q, p1) == H_OVER_R and peval(d.q, m1) == H_OVER_R),
             ("beta, beta' and the derivatives of 3t^2-2t^3 and t^2(1-t)^2 vanish at t = 0, 1; S5(0) = 0, S5(1) = 1, "
              "S5'(0) = S5'(1) = 0; 3t^2-2t^3 = 0, 1 at t = 0, 1",
              all(peval(BETA, x) == 0 and peval(pder(BETA), x) == 0 and peval(pder(H3), x) == 0
                  and peval(pder(G4), x) == 0 and peval(pder(S5), x) == 0 for x in (Fr(0), Fr(1)))
              and peval(S5, Fr(0)) == 0 and peval(S5, Fr(1)) == 1 and peval(H3, Fr(1)) == 1
              and peval(G4, Fr(1)) == 0)]
    for name, good in conds:
        log(f"  exact: {name}: {good}")
        ok &= good
    # phase rate at the cap ends: d theta/ds = (2 (R/b) q(+-1)) / (2 H |Y'(+-1)|) = 1/b
    rate = (2 * d.R / d.b * peval(d.q, p1)) / (2 * d.H * abs(peval(d.Yd, p1)))
    good = rate == 1 / d.b
    log(f"  exact: d theta/ds at both ends of each cap = 1/b (so the strands are C^1 at all four junctions): {good}")
    ok &= good
    log(f"  exact: Lambda(1) = {float(peval(d.lam, p1)):.12f}..., T0 = R Lambda(1)/b = {float(d.T0):.10f}..., "
        f"Theta_f = 2 T0 + E1 + E2 + Ew = {float(d.Thf):.10f}... (rationals)")
    # closing: total phase 2 L_s/b + Theta_f = n pi with L_s = (b/2)(n pi - Theta_f): an identity
    # sign of the axis curvature on a cap: N(v) = X'(v) Y''(v) - Y'(v) X''(v) <= 0 (R H > 0), zeros only at +-1
    Nv = padd(pmul(d.Xd, Yd2), pscale(pmul(d.Yd, Xd2), Fr(-1)))
    one_m_v2 = {0: Fr(1), 2: Fr(-1)}
    q_, k = Nv, 0
    while True:          # divide out (1 - v^2) while it divides exactly
        quo, rem = pdivmod(q_, one_m_v2)
        if rem:
            break
        q_, k = quo, k + 1
    ends = peval(q_, Fr(-1)) != 0 and peval(q_, Fr(1)) != 0
    nroots = sturm_count(q_, Fr(-1), Fr(1)) if ends else -1
    good = ends and nroots == 0 and peval(Nv, Fr(0)) < 0
    log(f"  exact: the axis curvature of a cap has the sign of X'Y'' - Y'X'' = (1 - v^2)^{k} Q(v) with Q of degree "
        f"{max(q_)} having no root in [-1, 1] (Sturm) and X'Y'' - Y'X'' < 0 at v = 0: the caps turn one way only, "
        f"so the axis is a convex curve: {good}")
    ok &= good
    return ok


def pdivmod(n, d):
    n = dict(n)
    dd = max(d)
    quo = {}
    while n and max(n) >= dd:
        k = max(n)
        f = n[k] / d[dd]
        quo[k - dd] = f
        for e, c in d.items():
            n[k - dd + e] = n.get(k - dd + e, 0) - f * c
        n = {e: c for e, c in n.items() if c}
    return quo, n


# ---------------------------------------------------------------------------------------------------------------
# helix lemma
def helix_lemma(d, log, target_ab=4.0, target_a=4.0, quiet=False):
    a1, a2, b = d.iv["a1"], d.iv["a2"], d.iv["b"]
    ok = True
    two_over_b = 2 / d.b
    # (H4) and c, exactly
    k1, k2 = d.a1 / (d.a1 ** 2 + d.b ** 2), d.a2 / (d.a2 ** 2 + d.b ** 2)
    good4 = k1 < 1 and k2 < 1
    if not quiet:
        log(f"  (H4) curvatures a_i/(a_i^2 + b^2) = {float(k1):.6f}, {float(k2):.6f} (exact rationals): both < 1: {good4}")
    else:
        log(f"  (H4): curvatures {float(k1):.6f}, {float(k2):.6f}: {'both < 1' if good4 else 'NOT both < 1: refused'}")
    ok &= good4

    def f_ab(U):
        s, c = sincos(U)
        return (a1 * a1 + a2 * a2) + (a1 * a2 * 2.0) * c + (b * b) * U.sqr(), (a1 * a2 * -2.0) * s + (b * b * 2.0) * U

    def f_a(U):
        s, c = sincos(U)
        return (a1 * a1 * 2.0) * (1.0 - c) + (b * b) * U.sqr(), (a1 * a1 * 2.0) * s + (b * b * 2.0) * U

    def minimize(f, lo, hi, target, n0=64, max_depth=44):
        """certified lower bound of f on [lo, hi] (bisection, mean-value form), and the least point value seen"""
        best_lb, best_pt = inf, (inf, None)
        stack = []
        w = (hi - lo) / n0
        for j in range(n0):
            stack.append((lo + j * w, lo + (j + 1) * w if j < n0 - 1 else hi, 0))
        while stack:
            a, bb, dep = stack.pop()
            m = 0.5 * (a + bb)
            fm, _ = f(Iv(m))
            _, fd = f(Iv(a, bb))
            enc = fm + fd * Iv(a - m, bb - m)
            if fm.hi < best_pt[0]:
                best_pt = (fm.hi, m)
            if enc.lo <= target and dep < max_depth and fm.hi >= 4.0:
                stack += [(a, m, dep + 1), (m, bb, dep + 1)]
                continue
            best_lb = min(best_lb, enc.lo)
        return best_lb, best_pt
    lb1, pt1 = minimize(f_ab, 0.0, float(Fr(two_over_b)) + 1e-9, target_ab)
    # the range |u| <= 2/b (trivial beyond: f_AB >= (a1 - a2)^2 + b^2 u^2 > 4)
    g1 = lb1 > 4.0
    u_min = ((ivq(d.a1 ** 2 + d.b ** 2)).sqrt()).inv() * PI
    lb2, pt2 = minimize(f_a, u_min.lo, float(Fr(two_over_b)) + 1e-9, target_a)
    g2 = lb2 > 4.0
    # (H3): pi / sqrt(a2^2 + b^2) > 2/b, so the range is empty and f_B >= b^2 u^2 > 4 there
    u3 = (PI / (ivq(d.a2 ** 2 + d.b ** 2)).sqrt()).lo
    g3 = u3 > float(two_over_b)
    if not quiet:
        log(f"  (H1) f_AB(u) > 4 for all u: certified min >= {lb1:.10f} on |u| <= 2/b (least point value "
            f"{pt1[0]:.10f} at u = {pt1[1]:.6f}); |u| >= 2/b is trivial: {g1}")
        log(f"  (H2) f_A(u) > 4 for |u| sqrt(a1^2 + b^2) >= pi: certified min >= {lb2:.10f} (least point value "
            f"{pt2[0]:.10f} at u = {pt2[1]:.6f}): {g2}")
        log(f"  (H3) pi/sqrt(a2^2 + b^2) = {u3:.6f}... > 2/b = {float(two_over_b):.6f}..., so (H3) is trivial: {g3}")
    else:
        log(f"  (H1): certified min of f_AB >= {lb1:.10f}; least point value {pt1[0]:.10f} at u = {pt1[1]:.6f}"
            f"{' < 4: a certified violation, refused' if pt1[0] < 4 else ''}")
    ok &= g1 and g2 and g3
    return ok, (lb1, pt1, lb2, pt2)


def c_bounds(d):
    """c = pi (sqrt(a1^2 + b^2) + sqrt(a2^2 + b^2)) in exact rationals (integer square roots, Machin bounds)"""
    def sqrt_bounds(x, digits=40):
        S = 10 ** digits
        n = x.numerator * S * S // x.denominator
        r = isqrt(n)
        return Fr(r, S), Fr(r + 1, S)
    s1 = sqrt_bounds(d.a1 ** 2 + d.b ** 2)
    s2 = sqrt_bounds(d.a2 ** 2 + d.b ** 2)
    return PI_LO * (s1[0] + s2[0]), PI_HI * (s1[1] + s2[1])


def isqrt(n):
    from math import isqrt as _i
    return _i(n)


# ---------------------------------------------------------------------------------------------------------------
def full_check(d, log, controls_mode=False):
    ok = True
    t0 = time.time()
    log("Exact identities (rationals):")
    ok &= exact_checks(d, log)
    log("Helix lemma (Section 3):")
    good, _ = helix_lemma(d, log, target_ab=4.0000055, target_a=4.0000257)
    ok &= good
    c_lo, c_hi = c_bounds(d)
    good = PRINTED["c_lo"] <= c_lo and c_hi <= PRINTED["c_up"]
    log(f"  c = pi (sqrt(a1^2 + b^2) + sqrt(a2^2 + b^2)) in [{float(c_lo):.12f}, {float(c_hi):.12f}] (exact bounds), "
        f"inside [7.2012718, 7.2012719]: {good}")
    ok &= good
    log("Zones (V2): curvature of the free pieces and every pair with a point on a free piece:")
    zones = [("top zone, n = 1 (mod 4)", top_zone(d, 1)), ("top zone, n = 3 (mod 4)", top_zone(d, 3)),
             ("bottom zone", bottom_zone(d))]
    for label, pcs in zones:
        for pc in pcs:
            if pc.free:
                (kA, kB), ev = curvature_bound(pc, 0.8, label)
                good = max(kA, kB) < 0.8
                log(f"  (V2) {label}, {pc.kind}: curvature of strand A <= {kA:.5f}, strand B <= {kB:.5f} "
                    f"(certified, {ev} interval evaluations): < 0.8: {good}")
                ok &= good
        good, stats, fails = check_zone(pcs, label, log)
        ok &= good
        if not good:
            return False
    log("Separation (V3) and the knot type (V4):")
    top1, top3, bot = zones[0][1], zones[1][1], zones[2][1]
    dt = max(height_bound(top1[1], 1, -1, target=1e-4), height_bound(top3[1], 1, -1, target=1e-4))
    yc2 = height_bound(bot[2], 1, 1, target=1e-4)
    good = dt <= 0.00012 and yc2 <= 0.00011
    log(f"  (S1) cap 1 (both strands, both phases) lies above y = L_s - {max(dt, 0):.2e}; (S2) cap 2 lies below y = "
        f"{max(yc2, 0):.2e} (certified); within the paper's delta_t = 0.00012 and y_cap2 = 0.00011: {good}")
    ok &= good
    # swap radii: exact polynomials in t; bisection
    rmax = swap_radius_bound(d)
    good = rmax <= 1.9511
    log(f"  (S2) swap strand radii <= {rmax:.6f} (certified): <= 1.9511: {good}; so 2R - a1 - 1.9511 = "
        f"{float(2 * d.R - d.a1) - 1.9511:.4f} > 2 separates the swap from column 1")
    ok &= good
    good = 2 * d.R - 2 * d.a1 > 2 and d.W - Fr("0.00012") > 2 and d.W - Fr("0.00011") > 2
    log(f"  (S1)-(S3) 2R - 2a1 = {float(2 * d.R - 2 * d.a1):.6f} > 2, W - delta_t > 2, W - y_cap2 > 2: {good}")
    ok &= good
    ls59 = (PI * 59.0 - ivq(d.Thf)) * ivq(d.b / 2)
    need = max(float(d.Sw + d.W), 0.00012 + float(max(d.Sw, Fr("0.00011"))) + 2, float(2 * d.R))
    good = ls59.lo > need
    ls57 = (PI * 57.0 - ivq(d.Thf)) * ivq(d.b / 2)
    log(f"  (S4) L_s(59) = (b/2)(59 pi - Theta_f) in [{ls59.lo:.6f}, {ls59.hi:.6f}] > max(S_w + W, delta_t + "
        f"max(S_w, y_cap2) + 2, 2R) = {need:.6f} (and L_s grows with n): {good}; L_s(57) = {ls57.mid():.4f}")
    ok &= good
    rho = axis_radius_bound(d)
    good = rho >= 6.0042
    log(f"  (V4) the axis radius of curvature is >= {rho:.5f} (certified): >= 6.0042: {good}; half-width R = 6.06; "
        f"strand radii <= max(a1, 1.9511) < 6.0042 (caps: a1 (1 + m_o beta) <= a1, a2 (1 + m_i beta) <= 1.54 a2)")
    ok &= good
    log("Length (V5):")
    Lsw, esw, ev1 = length_bound(bot[1])
    Lc2, ec2, ev2 = length_bound(bot[2])
    Cmax = None
    for case, pcs in ((1, top1), (3, top3)):
        Lc1, ec1, ev3 = length_bound(pcs[1])
        Lfree = Lc1[0] + Lc1[1] + Lsw[0] + Lsw[1] + Lc2[0] + Lc2[1]
        C = Lfree - (ivq(d.b * d.Thf + d.Sw) * d.ell)
        log(f"  n = {case} (mod 4): strand lengths cap 1 A {Lc1[0].lo:.6f}..{Lc1[0].hi:.6f}, B {Lc1[1].lo:.6f}.."
            f"{Lc1[1].hi:.6f}; swap A {Lsw[0].lo:.6f}..{Lsw[0].hi:.6f}, B {Lsw[1].lo:.6f}..{Lsw[1].hi:.6f}; cap 2 "
            f"A {Lc2[0].lo:.6f}..{Lc2[0].hi:.6f}, B {Lc2[1].lo:.6f}..{Lc2[1].hi:.6f}; L_free in [{Lfree.lo:.6f}, "
            f"{Lfree.hi:.6f}]; C = L_free - (b Theta_f + S_w) ell in [{C.lo:.6f}, {C.hi:.6f}]")
        Cmax = C.hi if Cmax is None else max(Cmax, C.hi)
    good = Cmax <= float(PRINTED["C"])
    log(f"  (V5) C <= {Cmax:.6f} for both entry phases: <= 19.649: {good}")
    ok &= good
    log(f"  (total {time.time() - t0:.0f} s)")
    return ok


def swap_radius_bound(d):
    rA = padd({0: d.a1}, pscale(S5, d.a2 - d.a1), pscale(BETA, d.mw))
    rB = padd({0: d.a2}, pscale(S5, d.a1 - d.a2), pscale(BETA, d.mw))
    best = -inf
    for p in (rA, rB):
        stack = [(j / 256, (j + 1) / 256, 0) for j in range(256)]
        dp = pder(p)
        while stack:
            a, b, dep = stack.pop()
            m = 0.5 * (a + b)
            vm = poly_jet(p, Jet.var(Iv(m), 0)).c[0]
            vd = poly_jet(dp, Jet.var(Iv(a, b), 0)).c[0]
            enc = vm + vd * Iv(a - m, b - m)
            if enc.hi > 1.9511 and dep < 30:
                stack += [(a, m, dep + 1), (m, b, dep + 1)]
                continue
            best = max(best, enc.hi)
    return best


def axis_radius_bound(d):
    """certified lower bound of the radius of curvature of the caps' axis: (x'^2 + y'^2)^(3/2)/|x'y'' - y'x''|"""
    Xd2, Yd2 = pder(d.Xd), pder(d.Yd)
    best = inf
    stack = [(-1 + 2 * j / 256, -1 + 2 * (j + 1) / 256, 0) for j in range(256)]
    while stack:
        a, b, dep = stack.pop()
        v = Jet.var(Iv(a, b), 0)
        xd = poly_jet(d.Xd, v).c[0] * ivq(d.R)
        yd = poly_jet(d.Yd, v).c[0] * ivq(d.H)
        xdd = poly_jet(Xd2, v).c[0] * ivq(d.R)
        ydd = poly_jet(Yd2, v).c[0] * ivq(d.H)
        num = (xd * xd + yd * yd)
        num = num * num.sqrt()
        den = (xd * ydd - yd * xdd).mag()
        rho = num.lo / up(den) if den > 0 else inf
        if rho < 6.0042 and dep < 30:
            m = 0.5 * (a + b)
            stack += [(a, m, dep + 1), (m, b, dep + 1)]
            continue
        best = min(best, rho)
    return best


def compare_sources(src, log):
    ok = True
    pdesign = os.path.join(src, "design.json")
    if os.path.exists(pdesign):
        dj = json.load(open(pdesign))
        flat = {**dj["helix"], **dj["cap"], **dj["swap"], "W": dj["W"]}
        same = flat == DESIGN
        log(f"  design.json read as data: every parameter string equals ours: {same}")
        ok &= same
    else:
        log("  (design.json not present; comparison skipped)")
    pq = os.path.join(src, "q_from_capaxis.json")
    if os.path.exists(pq):
        qj = {int(k): v for k, v in json.load(open(pq)).items()}
        same = qj == Q_STRINGS
        log(f"  the seven q strings taken from capaxis.py with take_q.py equal ours: {same}")
        ok &= same
    return ok


def main():
    log = print
    ok = True
    t0 = time.time()
    if "--sources" in sys.argv:
        src = sys.argv[sys.argv.index("--sources") + 1]
        log("Sources (read as data):")
        ok &= compare_sources(src, log)
    log("[Claude2026] C_22b <= 7.2012719* (PR #267): T(2, n), n >= 59 odd, thickness >= 1, length <= 7.2012719 n + 19.649")
    d = Design()
    ok &= full_check(d, log)
    if "--controls" in sys.argv:
        ok &= red_controls(log)
    log(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


def red_controls(log):
    log("Red controls (each must be refused):")
    ok = True
    # 1. a1 lowered: to the optimum's value (information: the margin of (H1) shrinks), then a little further (refused)
    d = Design({"a1": "1.9340602790"})
    good, (lb1, pt1, _, _) = helix_lemma(d, lambda s: None, quiet=True)
    log(f"  (information) a1 = 1.9340602790, the optimum's outer radius, with a2 and b unchanged: (H1) still holds, "
        f"least value of f_AB {pt1[0]:.10f} at u = {pt1[1]:.6f} (with the design's a1: 4.0000087...), because b is "
        f"still above the optimum's: the margin falls to about 2e-6")
    d = Design({"a1": "1.934055"})
    good, (lb1, pt1, _, _) = helix_lemma(d, lambda s: None, quiet=True)
    refused = (not good) and pt1[0] < 4
    log(f"  a1 = 1.934055 (7e-6 lower): f_AB({pt1[1]:.6f}) <= {pt1[0]:.10f} < 4, a certified violation of (H1): "
        f"refused: {refused}")
    ok &= refused
    # 2. Kim-Oh-Huh helix
    d = Design({"a1": "1.9498", "a2": "0.1420", "b": "0.3226"})
    k2 = d.a2 / (d.a2 ** 2 + d.b ** 2)
    refused = k2 > 1
    log(f"  the Kim-Oh-Huh helix (a1, a2, b) = (1.9498, 0.1420, 0.3226): inner curvature {float(k2):.5f} > 1 (radius "
        f"{float(1 / k2):.4f}), (H4) fails: refused: {refused}")
    ok &= refused
    # 3. W = 2: the windows are too short for the separation step
    d = Design({"W": "2"})
    refused = not (d.W - Fr("0.00012") > 2)
    log(f"  W = 2: W - delta_t = {float(d.W - Fr('0.00012')):.5f} is not > 2, so (S1) cannot separate cap 1 from the "
        f"column points outside the window: refused: {refused}")
    ok &= refused
    # 4. a perturbed q: q(+-1) != H/R, so d theta/ds jumps at the cap junctions
    d = Design(qover={8: "9.5170"})
    refused = peval(d.q, Fr(1)) != H_OVER_R
    log(f"  q with 9.5170 for 9.5169: q(1) = {float(peval(d.q, Fr(1)))} != 3.786, the phase rate jumps at the "
        f"junctions (not C^1): refused: {refused}")
    ok &= refused
    # 5., 6. a perturbed twist or outer-radius profile in cap 1: a certified pair below 2
    for over, name in (({"E1": "-19.5"}, "E1 = -19.5"), ({"mo": "0"}, "m_o = 0")):
        d = Design(over)
        lines = []
        good, stats, fails = check_zone(top_zone(d, 1), "top zone, n = 1 (mod 4)", lines.append, max_depth=16)
        wit = [f for f in fails if f[1] < 4.0]
        refused = (not good) and bool(wit)
        if wit:
            lb, ub, si, pi_i, a1_, b1_, sj, pi_j, a2_, b2_ = wit[0]
            log(f"  {name}: a certified violation in the top zone: strand {'AB'[si]} at tau = {0.5 * (a1_ + b1_):.7f} "
                f"({['window 1', 'cap 1', 'window 2'][pi_i]}) and strand {'AB'[sj]} at tau = {0.5 * (a2_ + b2_):.7f} "
                f"({['window 1', 'cap 1', 'window 2'][pi_j]}), |P(s) - P(t)|^2 <= {ub:.6f} < 4, arc separation "
                f">= pi: refused: {refused}")
        else:
            log(f"  {name}: no certified violation found: refused: {refused}")
        ok &= refused
    log(f"Red controls all refused: {ok}")
    return ok


if __name__ == "__main__":
    sys.exit(main())
