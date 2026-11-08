"""68a: the candidate lower bound C_68 >= 447/1000* of PR #299 (R. Pinho, the contact-row model), replayed from its
PROOF.md and CONTACT_ROW_LEMMA.md (rpinho/math-certificates @ 09c8423, c68/contact-row) and the two data files
certificate-0447.json and row-cover.json alone (Quantavane replays).  Standard library only; run with python3 -I.
The claimant's verify_measure.py and verify_row.py were neither read nor run.

What is checked, with our own code:
  1. The quantitative hypotheses at c = 447/1000 (w = -log c, A = pi^2/w, q = e^-A, Psi = (1+q)^2/(8q), kappa = 3/25):
     Psi > 6, sech^2(A/2) < 1/3, 16 sinh^2(A/4) < 2(1 - 8 kappa) Psi, tanh^2(A/8) > (1 + 8 kappa) pi^2/64,
     sech^2(A/8) < (1 - 4 kappa)/2  (110-digit decimal arithmetic, each value widened by 1e-80).
  2. The uniform linear-coefficient lemma (5) with z = 23/200 on the whole parameter region, from row-cover.json: the
     tree tiles the root box [0,17/50]^2 x [7/10,1] in (x, y, b) exactly (every split strictly inside, children rebuilt
     from the parent), and every leaf holds its stated condition, decided by our own bounds:
       PSD       1 - b^2 - 2(1-b)x^2 - 2(1+b)y^2 < 0 on the box (its maximum is at the low corner), exact rationals;
       constant  a lower bound of max(S, J) above z, exact rationals (sqrt bounded by integer square roots);
       small     x^2 + y^2 <= delta^2, delta = 10^-6, covered by the analytic small-region argument, whose numerical
                 side conditions are checked here (the 17 x 512 bins with rigorous lower bounds for their angular
                 mass and for log gamma_v);
       grid n    on n equal circle cells, W = max(S, J, T) bounded above and below over the whole cell and box (T
                 by Taylor bounds for cos, sin, cos 2theta around exact cell centres), then exp(mean of lower log W)
                 > z (1 + D^2) with D the triangle-inequality bound for |mean(log W e^-i theta)|.  Floating point
                 with every operation rounded outward by one ulp (math.nextafter); logarithms by Decimal.ln
                 (correctly rounded), widened one ulp; cell centres in 40-digit decimal.
     The two reductions outside the root box (a >= 17/50 gives S > z; b <= 7/10 gives S > z) are checked exactly.
  3. The moment conditions of section 5 for the 8-atom measure plus Lebesgue on [R, 1], R = 99/100: A_k >= 0 for
     every k (explicit up to the k where the atom at R alone suffices), B_k >= 0 for every k (explicit until one piece
     of [R, 1] alone beats all atoms), with the k = 1 term sum w_j eps(t_j) t_j, eps = z g (1 - 1/H)/2.  Our own lower
     bound for g on [R, 1] (the one of our #274 replay: g >= l_i (1 - t) on each of our 40 pieces), not the note's
     cubic; 320-bit fixed point with outward rounding.
  4. A sanity test of the contact-row lemma: random contact configurations (two or three contacts, random weights),
     checking max_i T(theta_i) >= J and the covariance inequality (6) after the rotation (floating point; a test of
     the prose, not part of the certificate).
--controls: z raised to 0.11503 (above the cover's minimum), kappa lowered, a split moved off the tiling, a grid leaf
coarsened, half of an atom's mass moved, c raised to 0.448, the eps term at k = 1 taken with the old 19/200.

What rests on mathematics (read, not computed): the Green/flux input and the angular majorant (1) (INPUT_GREEN_FLUX.md,
INPUT_ANGULAR.md), the minimal enclosing disk and its balanced contacts, the contact-row inequality (3) and the
moment floor J, the outer-function and Parseval steps behind (8), the small-region outer factors and the Jensen
identity, the automorphism identity (9) and the Bergman norm implication of section 6.

usage: python3 -I row_check.py CERTS_ROOT [--controls]
"""
import cmath
import hashlib
import json
import math
import os
import random
import sys
import time
from decimal import Decimal as D, getcontext, localcontext
from fractions import Fraction as F

getcontext().prec = 110
WIDEN = D("1e-80")
P = 320
ONE = 1 << P
CERT_SHA = "d9dd308f94efb019368682367609d9fa71c9982bb7195b27da22bb80d02f4afe"
COVER_SHA = "fda0aa1430313c3f896021073f3a44c3b1b920d2b07451d3b7d58d7d4d0bfc27"
NX = math.nextafter
INF = math.inf


# ------------------------------------------------------------------------------------------------ decimal helpers
def pi_dec():
    def arctan_inv(n):
        x = D(1) / n; x2 = x * x; term = x; s = x; k = 1
        while True:
            term *= -x2; k += 2
            t = term / k
            if abs(t) < D("1e-108"):
                return s
            s += t
    return 16 * arctan_inv(5) - 4 * arctan_inv(239)


PI = pi_dec()


def dec(fr):
    return D(fr.numerator) / D(fr.denominator)


def sin_cos_dec(v):
    """sin v, cos v for |v| <= 4 by their series (remainders below 1e-100)."""
    assert abs(v) <= 4
    s = v; c = D(1); ts = v; tc = D(1); n = 1
    while abs(ts) > D("1e-104") or abs(tc) > D("1e-104") or n < 8:
        ts = -ts * v * v / ((n + 1) * (n + 2))
        tc = -tc * v * v / (n * (n + 1))
        s += ts; c += tc; n += 2
    return s, c


def exp_dec(x):
    return x.exp()


def fx_floor(q):
    return (q.numerator << P) // q.denominator


def fx_ceil(q):
    return -((-q.numerator << P) // q.denominator)


def fx_dec_floor(x):
    return int((x * ONE).to_integral_value(rounding="ROUND_FLOOR"))


def fx_dec_ceil(x):
    return int((x * ONE).to_integral_value(rounding="ROUND_CEILING"))


def isqrt_lo(fr, bits=200):
    """a rational lower bound of sqrt(fr) (fr >= 0)."""
    n = (fr.numerator << (2 * bits)) // fr.denominator
    return F(math.isqrt(n), 1 << bits)


def isqrt_hi(fr, bits=200):
    n = -((-fr.numerator << (2 * bits)) // fr.denominator)
    r = math.isqrt(n)
    if r * r < n:
        r += 1
    return F(r, 1 << bits)


# ------------------------------------------------------------------------------------------------ 1. hypotheses
def envelope(c):
    w = -dec(c).ln()
    A = PI * PI / w
    q = (-A).exp()
    Psi = (1 + q) ** 2 / (8 * q)
    return w, A, Psi


def hypotheses(c, kappa, run):
    w, A, Psi = envelope(c)
    cosh = lambda x: (x.exp() + (-x).exp()) / 2
    sinh = lambda x: (x.exp() - (-x).exp()) / 2
    k = dec(kappa)
    v1 = 1 / cosh(A / 2) ** 2
    v2 = 16 * sinh(A / 4) ** 2; r2 = 2 * (1 - 8 * k) * Psi
    v3 = (sinh(A / 8) / cosh(A / 8)) ** 2; r3 = (1 + 8 * k) * PI * PI / 64
    v4 = 1 / cosh(A / 8) ** 2; r4 = (1 - 4 * k) / 2
    run("hypotheses at c = %s: Psi = %.10f > 6" % (c, Psi), Psi - WIDEN > 6)
    run("sech^2(A/2) = %.3e < 1/3" % v1, v1 + WIDEN < D(1) / 3)
    run("16 sinh^2(A/4) = %.6f < 2(1 - 8 kappa) Psi = %.6f" % (v2, r2), v2 + WIDEN < r2 - WIDEN)
    run("tanh^2(A/8) = %.10f > (1 + 8 kappa) pi^2/64 = %.10f" % (v3, r3), v3 - WIDEN > r3 + WIDEN)
    run("sech^2(A/8) = %.10f < (1 - 4 kappa)/2 = %.10f" % (v4, r4), v4 + WIDEN < r4)
    return w, A, Psi


# ------------------------------------------------------------------------------------------------ 2. the row cover
def dn(x):
    return NX(x, -INF)


def up(x):
    return NX(x, INF)


def fl_lo(fr):
    x = float(fr)
    return x if F(x) <= fr else dn(x)


def fl_hi(fr):
    x = float(fr)
    return x if F(x) >= fr else up(x)


def dec_to_fl_lo(d):
    x = float(d)
    return x if D(x) <= d else dn(x)


def dec_to_fl_hi(d):
    x = float(d)
    return x if D(x) >= d else up(x)


def log_lo(x):
    with localcontext() as ctx:
        ctx.prec = 25
        return dn(float(D(x).ln()))


def log_hi(x):
    with localcontext() as ctx:
        ctx.prec = 25
        return up(float(D(x).ln()))


def mul_iv(a_lo, a_hi, c_lo, c_hi):
    """[a_lo, a_hi] (a >= 0) times [c_lo, c_hi] (any sign), outward."""
    if c_lo >= 0:
        return dn(a_lo * c_lo), up(a_hi * c_hi)
    if c_hi <= 0:
        return dn(a_hi * c_lo), up(a_lo * c_hi)
    return dn(a_hi * c_lo), up(a_hi * c_hi)


TRIG = {}


def trig_table(n):
    """for each cell j: outward float bounds of cos, sin, cos 2theta over [theta_j - h0, theta_j + h0], theta_j =
    (2j+1) pi/n, h0 = pi/n, and float centres (cos, sin) with |error| <= 2^-52."""
    if n in TRIG:
        return TRIG[n]
    with localcontext() as ctx:
        ctx.prec = 40
        pi = +PI
        h0 = pi / n
        out = []
        for j in range(n):
            th = (2 * j + 1) * pi / n
            red = th - pi if th > pi else th          # reduce to (-pi, pi]: sin, cos flip sign
            s, c = sin_cos_dec(red)
            if th > pi:
                s, c = -s, -c
            s2 = 2 * s * c; c2 = c * c - s * s
            ec = abs(s) * h0 + h0 * h0 / 2 + D("1e-35")
            es = abs(c) * h0 + h0 * h0 / 2 + D("1e-35")
            e2 = 2 * abs(s2) * h0 + 2 * h0 * h0 + D("1e-35")
            cl = max(dec_to_fl_lo(c - ec), -1.0); ch = min(dec_to_fl_hi(c + ec), 1.0)
            sl = max(dec_to_fl_lo(s - es), -1.0); sh = min(dec_to_fl_hi(s + es), 1.0)
            c2l = max(dec_to_fl_lo(c2 - e2), -1.0); c2h = min(dec_to_fl_hi(c2 + e2), 1.0)
            out.append((cl, ch, sl, sh, c2l, c2h, float(c), float(s)))
        h0f = dec_to_fl_hi(h0)
    TRIG[n] = (out, h0f)
    return TRIG[n]


def sj_bounds(box, kappa):
    """exact rational lower/upper bounds of S and J on a box (x, y >= 0)."""
    xl, xh, yl, yh, bl, bh = box
    k2 = 2 * kappa
    Slo = xl * xl + yl * yl + k2 * (1 - bh * bh)
    Shi = xh * xh + yh * yh + k2 * (1 - bl * bl)
    al2, ah2 = xl * xl + yl * yl, xh * xh + yh * yh
    if ah2 == 0:
        return Slo, Shi, F(0), F(0)
    ah = isqrt_hi(ah2)
    Jlo = (al2 + ((1 + bl) * xl * xl + (1 - bh) * yl * yl) / (2 * ah)) / (1 + ah)
    if al2 == 0:
        Jhi = ah
    else:
        al = isqrt_lo(al2)
        Jhi = min(ah, (ah2 + ((1 + bh) * xh * xh + (1 - bl) * yh * yh) / (2 * al)) / (1 + al))
    return Slo, Shi, Jlo, Jhi


def grid_leaf(box, n, kappa):
    """returns (ratio lower bound as float, ok) for the scalar criterion on the box with n cells; ok means
    exp(mean l) > z (1 + D^2) is decided later by the caller with exact z."""
    Slo, Shi, Jlo, Jhi = sj_bounds(box, kappa)
    base_lo = fl_lo(max(Slo, Jlo)); base_hi = fl_hi(max(Shi, Jhi))
    xl, xh, yl, yh, bl, bh = [x for x in box]
    xlf, xhf, ylf, yhf, blf, bhf = fl_lo(xl), fl_hi(xh), fl_lo(yl), fl_hi(yh), fl_lo(bl), fl_hi(bh)
    k2l, k2h = fl_lo(2 * kappa), fl_hi(2 * kappa)
    cells, h0 = trig_table(n)
    L = []; Mre = []; Mim = []; Esum = []; absM = 0.0
    for (cl, ch, sl, sh, c2l, c2h, cc, ss) in cells:
        a1, a2 = mul_iv(xlf, xhf, cl, ch)
        b1, b2 = mul_iv(ylf, yhf, sl, sh)
        d1, d2 = mul_iv(blf, bhf, c2l, c2h)            # b cos 2theta
        om_lo, om_hi = dn(1.0 - d2), up(1.0 - d1)      # 1 - b cos 2theta
        k_lo = dn((k2l if om_lo >= 0 else k2h) * om_lo)
        k_hi = up((k2h if om_hi >= 0 else k2l) * om_hi)
        Tlo = dn(dn(a1 + b1) + k_lo)
        Thi = up(up(a2 + b2) + k_hi)
        wlo = max(base_lo, Tlo); whi = max(base_hi, Thi)
        if not wlo > 0:
            return None
        lj = log_lo(wlo); uj = log_hi(whi)
        Mj = (lj + uj) / 2; Ej = max(up(uj - Mj), up(Mj - lj))
        L.append(lj); Mre.append(Mj * cc); Mim.append(-Mj * ss)
        Esum.append(up(Ej + up(abs(Mj) * h0)))
        absM += abs(Mj)
    nn = float(n)
    mean_l = dn(dn(math.fsum(L)) / nn)
    re = math.fsum(Mre); im = math.fsum(Mim)
    err = up(absM * 2 ** -48)                         # products, table and summation errors, generously
    harm = up(up(math.hypot(abs(re) + err, abs(im) + err)) * up(1 + 2 ** -50)) / nn
    Dbound = up(up(harm) + up(math.fsum(Esum) / nn) * up(1 + 2 ** -50))
    Dbound = up(Dbound * (1 + 2 ** -50))
    return mean_l, Dbound


def criterion(mean_l, Dbound, z):
    """exp(mean_l) > z (1 + D^2), decided in 50-digit decimal with margins; returns (ok, ratio)."""
    with localcontext() as ctx:
        ctx.prec = 50
        e = D(mean_l).exp() * (1 - D("1e-40"))
        rhs = dec(z) * (1 + D(Dbound) ** 2) * (1 + D("1e-40"))
        ratio = e / (1 + D(Dbound) ** 2)
        return e > rhs, ratio


def small_region(z, delta, run, quiet=False):
    """the numerical side conditions of the small-region argument (section 3)."""
    ok = True
    ell = F(12, 25) * F(999, 1000)
    ok &= run("small region, b <= 999/1000: p = (6/25)(1 - 999/1000) - delta > 0 and 3/25 - delta/2 > z",
              F(6, 25) * F(1, 1000) - delta > 0 and F(3, 25) - delta / 2 > z)
    X0 = delta * (1 + delta) / ell; d0 = X0 / 2
    ok &= run("small region, b >= 999/1000: ell = 2997/6250, ell/4 > z (a = 0), X0 = a(1+a)/ell < 1 and "
              "4 ell > (1+a)^2 for a <= delta", ell / 4 > z and X0 < 1 and 4 * ell > (1 + delta) ** 2)
    ends = [F(1, 2 ** e) for e in range(17, 0, -1)] + [F(16383, 16384)]
    with localcontext() as ctx:
        ctx.prec = 40
        pi = +PI
        total = D(0); bins = 0; side = True
        for u0, v0 in zip(ends, ends[1:]):
            for i in range(512):
                u = u0 + (v0 - u0) * i / 512; v = u0 + (v0 - u0) * (i + 1) / 512
                gam = 1 / v - delta / (ell * v * v)
                side &= u > d0 and u > 2 * delta / ell and gam > 1
                # angular mass (2/pi)(asin sqrt v - asin sqrt u) >= (2/pi)(v - u)/(2 sqrt(max X(1-X) on [u, v]))
                mx = v * (1 - v) if v <= F(1, 2) else (u * (1 - u) if u >= F(1, 2) else F(1, 4))
                mass_lo = (dec(v - u) / (2 * dec(isqrt_hi(mx)))) * 2 / pi * (1 - D("1e-35"))
                lg = dec(gam).ln() * (1 - D("1e-35"))
                total += mass_lo * lg
                bins += 1
        lhs = (dec(ell) / 16).ln() + total - D("1e-30")
        rhs = dec(z).ln() + D("1e-30")
    ok &= run("small region: %d bins, u > d0, u > 2 delta/ell and gamma_v > 1 on every bin" % bins, side)
    ok &= run("small region: log(ell/16) + sum mass * log gamma_v = %.9f > log z = %.9f" % (lhs, rhs), lhs > rhs)
    return ok


def walk_cover(cover, z, kappa, run, quiet=False, report=None):
    """tile check and leaf conditions; returns True if all hold."""
    delta = F(cover["delta"])
    root = [F(x) for x in cover["root"]]
    stats = {"grid": 0, "psd": 0, "constant": 0, "small": 0}
    bad = []
    worst = [None]
    leaves = []

    def rec(node, box):
        if "kind" in node:
            leaves.append((node, box)); return
        ax = node["axis"]; mid = F(node["mid"])
        lo, hi = box[2 * ax], box[2 * ax + 1]
        if not (lo < mid < hi) or set(node) != {"axis", "mid", "left", "right"}:
            bad.append(("split", box, node.get("mid"))); return
        left = list(box); left[2 * ax + 1] = mid
        right = list(box); right[2 * ax] = mid
        rec(node["left"], left); rec(node["right"], right)

    rec(cover["tree"], root)
    tile_ok = not bad
    # the leaves tile the root: disjoint interiors and total volume equal (they come from a binary split tree)
    vol = sum((b[1] - b[0]) * (b[3] - b[2]) * (b[5] - b[4]) for _, b in leaves)
    rvol = (root[1] - root[0]) * (root[3] - root[2]) * (root[5] - root[4])
    tile_ok &= vol == rvol
    run("row cover: the tree tiles the root box %s exactly (%d leaves, every split strictly inside)" %
        ("[0,17/50]^2 x [7/10,1]" if root == [0, F(17, 50), 0, F(17, 50), F(7, 10), 1] else root, len(leaves)),
        tile_ok and root == [0, F(17, 50), 0, F(17, 50), F(7, 10), 1])
    for node, box in leaves:
        kind = node["kind"]; xl, xh, yl, yh, bl, bh = box
        if kind == "psd" and set(node) == {"kind"}:
            ok = 1 - bl * bl - 2 * (1 - bl) * xl * xl - 2 * (1 + bl) * yl * yl < 0 and bl > xh * xh
        elif kind == "constant" and set(node) == {"kind"}:
            Slo, _, Jlo, _ = sj_bounds(box, kappa)
            ok = max(Slo, Jlo) > z
        elif kind == "small" and set(node) == {"kind"}:
            ok = xh * xh + yh * yh <= delta * delta
        elif kind == "grid" and set(node) == {"kind", "n"} and isinstance(node["n"], int) and node["n"] >= 8:
            res = grid_leaf(box, node["n"], kappa)
            if res is None:
                ok = False
            else:
                ok, ratio = criterion(res[0], res[1], z)
                if worst[0] is None or ratio < worst[0][0]:
                    worst[0] = (ratio, box, node["n"])
        else:
            ok = False
        stats[kind] = stats.get(kind, 0) + 1
        if not ok:
            bad.append((kind, box))
    if report is not None:
        report.update(stats=stats, worst=worst[0], bad=bad)
    return not bad


# ------------------------------------------------------------------------------------------------ 3. moments
def g_H(t, w, Psi):
    """g = 2/(1+H), H = sqrt(1 + Psi/(2 sigma)), sigma = sin^2(pi log t/(2w)); returns (g, H) or (0, None) at t = c^2."""
    v = -PI * t.ln() / (2 * w)
    s, _ = sin_cos_dec(v) if v <= 4 else (None, None)
    sig = s * s
    if sig <= D("1e-90"):
        return D(0), None
    H = (1 + Psi / (2 * sig)).sqrt()
    return 2 / (1 + H), H


def moments(cert, z_k1, run, w, Psi, quiet=False):
    c = F(cert["c"]); R = F(cert["R"])
    t = [F(x) for x in cert["nodes"]]; wts = [F(x) for x in cert["weights"]]
    ok = True
    ok &= run("measure: %d atoms, increasing, t_1 = c^2, t_8 = R = 99/100, weights > 0, sum w_j = R exactly" % len(t),
              len(t) == len(wts) == 8 and t == sorted(t) and len(set(t)) == 8 and t[0] == c * c and t[-1] == R
              and all(x > 0 for x in wts) and sum(wts) == R)
    beta_hi = []; eps_lo = []
    for x in t:
        g, H = g_H(dec(x), w, Psi)
        if H is None:
            beta_hi.append(fx_dec_ceil(D(1))); eps_lo.append(0); continue
        gl = g - WIDEN
        beta_hi.append(fx_dec_ceil(1 - gl))
        e = dec(z_k1) * (g / 2) * (1 - 1 / H) - WIDEN
        eps_lo.append(fx_dec_floor(max(e, D(0))))
    # our partition of [R, 1] (geometric in 1 - t) and our slope bounds g >= l_i (1 - t) on each piece
    ds = [F(1, 100)]
    while ds[-1] > F(1, 10 ** 6):
        ds.append(ds[-1] * F(4, 5))
    ds = [d for d in ds if d > F(1, 10 ** 6)] + [F(0)]
    tau = [1 - d for d in ds]
    slopes = []
    for a_, b_ in zip(tau, tau[1:]):
        dlo, dhi = 1 - dec(b_), 1 - dec(a_)
        U = -PI * dec(a_).ln() / (2 * w)
        if not U * U < 6:
            run("tail piece slope (U^2 < 6)", False); return False
        bi = PI / (2 * w) * (1 + dlo / 2) * (1 - U * U / 6)
        li = 2 * D(2).sqrt() * bi / (D(2).sqrt() * bi * dhi + (2 * bi * bi * dhi * dhi + Psi).sqrt()) - WIDEN
        slopes.append(fx_dec_floor(li))
    lo = lambda x, y: x * y >> P
    hi = lambda x, y: -((-x * y) >> P)
    tl = [fx_floor(x) for x in t]; th = [fx_ceil(x) for x in t]
    wl = [fx_floor(x) for x in wts]; wh = [fx_ceil(x) for x in wts]
    taul = [fx_floor(x) for x in tau]; tauh = [fx_ceil(x) for x in tau]
    Rl, Rh = fx_floor(R), fx_ceil(R)
    pl = [ONE] * 8; ph = [ONE] * 8
    Tl = list(taul); Th = list(tauh)
    Rkl, Rkh = Rl, Rh
    K1 = 0
    while wts[-1] * (K1 + 1) < R:
        K1 += 1
    minI = minII = None; K2 = None; k = 0; tail = None
    while k < K1 or K2 is None:
        k += 1
        pl = [lo(pl[j], tl[j]) for j in range(8)]
        ph = [hi(ph[j], th[j]) for j in range(8)]
        Rkl, Rkh = lo(Rkl, Rl), hi(Rkh, Rh)
        Tl = [lo(Tl[i], taul[i]) for i in range(len(tau))]
        Th = [hi(Th[i], tauh[i]) for i in range(len(tau))]
        target_hi = -((-Rkh) // (k + 1))
        if k < K1:
            s_lo = sum(lo(wl[j], pl[j]) for j in range(8))
            if s_lo < target_hi:
                ok = False; run("A_%d >= 0" % k, False); break
            m = (s_lo - target_hi) / target_hi
            if minI is None or m < minI[0]:
                minI = (m, k)
        if K2 is None:
            s_atoms = sum(hi(wh[j], ph[j]) for j in range(8))
            for i in range(len(tau) - 1):
                a_k_lo = (Tl[i] << P) // tauh[i]
                piece = lo(lo(slopes[i], fx_floor((1 - tau[i + 1]) * (tau[i + 1] - tau[i]))), a_k_lo)
                if piece > s_atoms:
                    K2 = k; tail = (i, piece, s_atoms); break
            if K2 is None:
                s_hi = sum(hi(hi(wh[j], ph[j]), beta_hi[j]) for j in range(8))
                integ = 0
                for i in range(len(tau) - 1):
                    b1l, b1h = Tl[i + 1], Th[i + 1]; a1l, a1h = Tl[i], Th[i]
                    up_lo = b1l // (k + 1) - (-((-a1h) // (k + 1)))
                    dn_hi = -((-(hi(b1h, tauh[i + 1]) - lo(a1l, taul[i]))) // (k + 2))
                    integ += lo(slopes[i], max(up_lo - dn_hi, 0))
                lhs_lo = integ + Rkl // (k + 1)
                if k == 1:
                    lhs_lo += sum(lo(lo(wl[j], tl[j]), eps_lo[j]) for j in range(8))
                if lhs_lo < s_hi:
                    ok = False; run("B_%d >= 0" % k, False); break
                m = (lhs_lo - s_hi) / target_hi
                if minII is None or m < minII[0]:
                    minII = (m, k)
                if k == 1:
                    b1margin = (lhs_lo - s_hi) / ONE
        if k > 10 ** 5:
            ok = False; run("tails found below k = 10^5", False); break
    if ok:
        i, piece, sat = tail
        run("A_k >= 0 for k = 1..%d explicitly (smallest relative margin %.4e at k = %d); from k = %d the atom at R "
            "alone (w_8 (k+1) >= R)" % (K1 - 1, minI[0], minI[1], K1), True)
        run("B_k >= 0 for k = 1..%d explicitly (smallest relative margin %.4e at k = %d; B_1 >= %.6e with the "
            "atomic eps term), and from k = %d the piece [%s, %s] alone gives %.4e > sum w_j t_j^k (%.4e), "
            "for every larger k too" % (K2 - 1, minII[0], minII[1], b1margin, K2, tau[i], tau[i + 1], piece / ONE,
                                       sat / ONE), True)
    return ok


# ------------------------------------------------------------------------------------------------ 4. lemma sanity
def lemma_sanity(kappa, trials=20000, seed=299):
    rng = random.Random(seed)
    worst_row = INF; worst_psd = INF
    for _ in range(trials):
        m = rng.choice((2, 3))
        th = [rng.uniform(0, 2 * math.pi) for _ in range(m)]
        lam = [rng.random() for _ in range(m)]; s = sum(lam); lam = [x / s for x in lam]
        m1 = sum(l * cmath.exp(1j * t) for l, t in zip(lam, th))
        m2 = sum(l * cmath.exp(2j * t) for l, t in zip(lam, th))
        a = abs(m1)
        if a < 1e-9:
            continue
        T = [(cmath.exp(1j * t) * m1.conjugate()).real + 2 * kappa * (1 - (cmath.exp(2j * t) * m2.conjugate()).real)
             for t in th]
        q0 = (1 + (m2 * m1.conjugate() ** 2 / a ** 2).real) / 2
        J = a * (q0 + a) / (1 + a)
        worst_row = min(worst_row, max(T) - J)
        # rotate: m2 -> b >= 0, then reflect so x, y >= 0
        phi = cmath.phase(m2) / 2 if abs(m2) > 1e-12 else 0.0
        mm1 = m1 * cmath.exp(-1j * phi); b = abs(m2)
        x, y = abs(mm1.real), abs(mm1.imag)
        worst_psd = min(worst_psd, 1 - b * b - 2 * (1 - b) * x * x - 2 * (1 + b) * y * y)
        J7 = (a * a + ((1 + b) * x * x + (1 - b) * y * y) / (2 * a)) / (1 + a)
        if abs(J7 - J) > 1e-9:
            return False, worst_row, worst_psd
    return worst_row > -1e-9 and worst_psd > -1e-9, worst_row, worst_psd


# ------------------------------------------------------------------------------------------------ main
class Run:
    def __init__(self, quiet=False):
        self.n = 0; self.fails = []; self.quiet = quiet

    def __call__(self, label, ok):
        self.n += 1
        if not ok:
            self.fails.append(label)
        if not self.quiet:
            print(f"  check {label}: {'PASS' if ok else 'FAIL'}")
        return ok


def full(cert, cover, run, quiet=False, report=None, skip_moments=False, skip_cover=False, z_k1=None):
    z = F(cover["zeta"]); kappa = F(cover["kappa"])
    c = F(cert["c"])
    w, A, Psi = hypotheses(c, kappa, run)
    run("the certificate and the cover use the same z = %s and kappa = %s" % (z, kappa),
        F(cert["zeta"]) == z and F(cert["kappa"]) == kappa and cert["cover_sha256"] == COVER_SHA)
    run("outside the root box S > z: (17/50)^2 = 289/2500 > z and (6/25)(1 - (7/10)^2) = 306/2500 > z",
        F(17, 50) ** 2 > z and F(6, 25) * (1 - F(7, 10) ** 2) > z)
    if not skip_cover:
        rep = {} if report is None else report
        small_region(z, F(cover["delta"]), run)
        okc = walk_cover(cover, z, kappa, run, report=rep)
        st = rep["stats"]; worst = rep["worst"]
        run("row cover: every leaf holds (%d grid, %d PSD, %d constant, %d small)%s" %
            (st.get("grid", 0), st.get("psd", 0), st.get("constant", 0), st.get("small", 0),
             "; smallest grid ratio exp(mean log W)/(1 + D^2) = %.8f at n = %d" % (worst[0], worst[2]) if worst else ""),
            okc)
    if not skip_moments:
        moments(cert, z if z_k1 is None else z_k1, run, w, Psi)


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    if len(args) != 1:
        print(__doc__); sys.exit(2)
    t0 = time.time()
    root = os.path.join(args[0], "pin0447")
    rawc = open(os.path.join(root, "certificate-0447.json"), "rb").read()
    rawv = open(os.path.join(root, "row-cover.json"), "rb").read()
    print("68a C_68 >= 447/1000 (PR #299, contact-row model), from PROOF.md and the two data files")
    run = Run()
    run("certificate-0447.json sha256 = d9dd308f...2afe (the packet's SHA256SUMS)",
        hashlib.sha256(rawc).hexdigest() == CERT_SHA)
    run("row-cover.json sha256 = fda0aa14...fc27 (the packet's SHA256SUMS and the certificate's cover_sha256)",
        hashlib.sha256(rawv).hexdigest() == COVER_SHA)
    cert = json.loads(rawc); cover = json.loads(rawv)
    full(cert, cover, run)
    okl, wr, wp = lemma_sanity(float(F(cover["kappa"])))
    run("contact-row lemma, sanity test on 20,000 random 2- and 3-contact configurations: max_i T(theta_i) - J >= "
        "%.3e >= 0, inequality (6) >= %.3e >= 0, and (7) equals a(q0 + a)/(1 + a)" % (wr, wp), okl)
    print(f"  (main run: {time.time() - t0:.0f} s)")
    if "--controls" in sys.argv:
        print("  -- red controls (each altered input must be refused)")

        def refused(label, cert2, cover2, **kw):
            sub = Run(quiet=True)
            full(cert2, cover2, sub, **kw)
            run("control: " + label + " is refused" + (" (first failing check: %s)" % sub.fails[0][:90]
                                                         if sub.fails else ""), bool(sub.fails))

        cv = json.loads(rawv); cv["zeta"] = "11503/100000"
        ce = json.loads(rawc); ce["zeta"] = "11503/100000"
        refused("z raised to 0.11503 (just above the cover's smallest ratio)", ce, cv, skip_moments=True)
        cv = json.loads(rawv); cv["kappa"] = "29/250"
        ce = json.loads(rawc); ce["kappa"] = "29/250"
        refused("kappa lowered from 3/25 to 29/250", ce, cv, skip_moments=True)
        cv = json.loads(rawv); cv["tree"]["left"]["mid"] = "2/5"
        refused("a split point moved outside its box (17/100 to 2/5, the boxes no longer tile)",
                json.loads(rawc), cv, skip_moments=True)
        cv = json.loads(rawv)

        def coarsen(node):
            if "kind" in node:
                if node["kind"] == "grid":
                    node["n"] = 64; return True
                return False
            return coarsen(node["left"]) or coarsen(node["right"])
        coarsen(cv["tree"])
        refused("the first grid leaf coarsened from 4096 to 64 cells", json.loads(rawc), cv, skip_moments=True)
        ce = json.loads(rawc); wv = [F(x) for x in ce["weights"]]; s = wv[1] / 2; wv[1] -= s; wv[2] += s
        ce["weights"] = [str(x) for x in wv]
        refused("half of atom 2's mass moved to atom 3", ce, json.loads(rawv), skip_cover=True)
        ce = json.loads(rawc); ce["c"] = "56/125"; ce["nodes"][0] = str(F(56, 125) ** 2)
        refused("c raised to 0.448 (first node at the new c^2)", ce, json.loads(rawv), skip_cover=True)
        refused("the k = 1 eps term with the old constant 19/200 instead of 23/200", json.loads(rawc),
                json.loads(rawv), skip_cover=True, z_k1=F(19, 200))
    print(f"  checks: {run.n}, failed: {len(run.fails)}")
    print(f"  run time: {time.time() - t0:.0f} s")
    print("RESULT: ALL CHECKS PASS" if not run.fails else "RESULT: FAIL " + "; ".join(run.fails))
    sys.exit(0 if not run.fails else 1)


if __name__ == "__main__":
    main()
