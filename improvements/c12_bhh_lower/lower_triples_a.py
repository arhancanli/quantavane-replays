#!/usr/bin/env python3
"""MST constant in the plane: certified lower bound with isolated points, pairs and TRIPLES (program A).

Extends lower_pairs_a.py (same program family: its interval arithmetic, pi, erf, exp and the head, point and pair
terms are imported from it). Standard library only.

Bound: beta >= int_0^a chi + int_a^inf (n1 + n2 + n3), n3 = isolated triples per point,
  n3(r) = (1/6) int 1[G_r(0,y2,y3) connected] e^{-A_r(0,y2,y3)} dy2 dy3.
Every connected triple has a vertex adjacent to the other two (the middle of a path: one labelling in three; any
vertex of a triangle). Labelling that vertex 1 (at 0), with a = |y2| < 1, b = |y3| < 1, theta the angle between them
and c = |y2 - y3|, after the r-integral (r = scale):
  int_a^inf n3 = 4 pi int_0^1 int_0^1 int_0^pi  w(c) G3(A(a,b,c)) a b dtheta db da,
  w = 1/2 for a path (c >= 1), 1/6 for a triangle (c < 1),
  G3(A) = int_acut^inf r^4 e^{-A r^2} dr = acut^3 e^{-A acut^2}/(2A) + 3 acut e^{-A acut^2}/(4A^2)
          + 3 sqrt(pi) erfc(acut sqrt A)/(8 A^{5/2}),
  A(a,b,c) = area of the union of three unit discs with mutual distances a, b, c.
By the planar Kneser-Poulsen theorem (Bezdek and Connelly 2002) A does not decrease when no distance decreases, so on a
cell [a0,a1] x [b0,b1] x [t0,t1], A <= A(a1, b1, c+) for any c+ >= max c on the cell (realised by an expanded,
possibly collinear, triangle). G3 decreases in A. On each cell: w >= 1/2 if min c >= 1, else w >= 1/6;
weight int a b dtheta = (a1^2-a0^2)/2 (b1^2-b0^2)/2 (t1-t0).
The union area is computed by Green's theorem over the uncovered arcs of each circle, in interval arithmetic, with
the arc end points' sines and cosines from the geometry (no trigonometric series), and the arc angles from acos.

Usage: python3 -I lower_triples_a.py [--a 2107/2000] [--na 128] [--nt 256] [--nc 128] [--pcells 200000]
                                     [--claim 0.6339] [--controls]
"""
import sys, time, math, os
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # lower_pairs_a.py sits next to this file
import lower_pairs_a as LA
from lower_pairs_a import I, ONE, PI, sqrt, exp_neg, erf_pos, asin_half, fdiv, cdiv


def hull(x, y):
    return I(min(x.lo, y.lo), max(x.hi, y.hi))


def asin_any(t):
    """asin on a rational or an interval inside [0, 1/2] (via asin_half)."""
    return asin_half(t)


def acos_q(t):
    """acos of an exact rational t in [-1, 1]."""
    t = Fraction(t)
    assert -1 <= t <= 1
    if abs(t) <= Fraction(1, 2):
        s = asin_half(I.frac(abs(t)))
        return PI * Fraction(1, 2) - s if t >= 0 else PI * Fraction(1, 2) + s
    w = sqrt(I.frac((1 - abs(t)) / 2))
    if w.hi * 2 > ONE:  # cannot happen for |t| > 1/2 except rounding at the boundary
        w = I(w.lo, ONE // 2)
    s = 2 * asin_half(w)
    return s if t > 0 else PI - s


def green(cx, cy, s, t, ss, cs, st, ct):
    """(1/2)[(t - s) + cx (sin t - sin s) - cy (cos t - cos s)] for an arc of a unit circle centred at (cx, cy)."""
    return ((t - s) + cx * (st - ss) - cy * (ct - cs)) * Fraction(1, 2)


class Arc:
    """Arc [p - al, p + al] of circle i covered by disc j: angle p (interval), half-width al, and the sines and
    cosines of both ends from sin p, cos p, cos al = d/2, sin al = sqrt(1 - d^2/4)."""

    def __init__(self, p, sp, cp, d):
        self.ok = d < 2
        if not self.ok:
            return
        al = acos_q(d / 2)
        ca = I.frac(d / 2)
        sa = sqrt(I.frac(1 - d * d / 4))
        self.lo, self.hi = p - al, p + al
        self.s_lo, self.c_lo = sp * ca - cp * sa, cp * ca + sp * sa
        self.s_hi, self.c_hi = sp * ca + cp * sa, cp * ca - sp * sa


def shifted(arc, k):
    """The same arc with its angles moved by 2 pi k (sines and cosines unchanged)."""
    o = Arc.__new__(Arc)
    o.ok = arc.ok
    o.lo, o.hi = arc.lo + PI * (2 * k), arc.hi + PI * (2 * k)
    o.s_lo, o.c_lo, o.s_hi, o.c_hi = arc.s_lo, arc.c_lo, arc.s_hi, arc.c_hi
    return o


def inter_green(cx, cy, J, K):
    """Green integral over J intersect K (arcs of length <= pi each): K moved by a multiple of 2 pi so that its
    centre is within pi of J's; ambiguous comparisons are resolved by taking the hull of every possibility."""
    results = []
    cj = (J.lo + J.hi) * Fraction(1, 2)
    for k in (-1, 0, 1):
        K2 = shifted(K, k)
        ck = (K2.lo + K2.hi) * Fraction(1, 2)
        dl = ck - cj
        if dl.hi < -PI.hi or dl.lo > PI.hi:  # surely not the representative
            continue
        # start = max(J.lo, K2.lo), end = min(J.hi, K2.hi)
        starts = []
        if not (J.lo.hi < K2.lo.lo):
            starts.append((J.lo, J.s_lo, J.c_lo))
        if not (K2.lo.hi < J.lo.lo):
            starts.append((K2.lo, K2.s_lo, K2.c_lo))
        ends = []
        if not (J.hi.lo > K2.hi.hi):
            ends.append((J.hi, J.s_hi, J.c_hi))
        if not (K2.hi.lo > J.hi.hi):
            ends.append((K2.hi, K2.s_hi, K2.c_hi))
        for (s, ss, cs) in starts:
            for (t, st, ct) in ends:
                if t.hi < s.lo:
                    results.append(I(0))
                elif t.lo > s.hi:
                    results.append(green(cx, cy, s, t, ss, cs, st, ct))
                else:
                    results.append(hull(I(0), green(cx, cy, s, t, ss, cs, st, ct)))
    assert results
    out = results[0]
    for r in results[1:]:
        out = hull(out, r)
    return out


def union3(a, b, c):
    """Enclosure of the area of the union of unit discs at (0,0), (a,0), (x,y) with |y3| = b, |y3 - y2| = c
    (exact rationals, a valid possibly degenerate triangle, a, b, c > 0)."""
    x = (a * a + b * b - c * c) / (2 * a)
    y2 = b * b - x * x
    assert y2 >= 0
    y = sqrt(I.frac(y2))
    # directions: 1->2 angle 0; 1->3 angle p13 = acos(x/b); 2->3 angle p23 = acos((x - a)/c)
    p13 = acos_q(x / b)
    p23 = acos_q((x - a) / c)
    s13, c13 = y / I.frac(b), I.frac(x / b)
    s23, c23 = y / I.frac(c), I.frac((x - a) / c)
    zero, one = I(0), I(ONE)
    # circle 1 (0,0): covered by disc 2 (angle 0, d = a) and disc 3 (p13, d = b)
    circles = [
        ((I(0), I(0)), Arc(zero, zero, one, a), Arc(p13, s13, c13, b)),
        ((I.frac(a), I(0)), Arc(PI, zero, -one, a), Arc(p23, s23, c23, c)),
        ((I.frac(x), y), Arc(p13 + PI, -s13, -c13, b), Arc(p23 + PI, -s23, -c23, c)),
    ]
    total = I(0)
    for (cx, cy), J, K in circles:
        cov = I(0)
        for X in (J, K):
            if X.ok:
                cov = cov + green(cx, cy, X.lo, X.hi, X.s_lo, X.c_lo, X.s_hi, X.c_hi)
        if J.ok and K.ok:
            cov = cov - inter_green(cx, cy, J, K)
        total = total + (PI - cov)
    return total


def G3(acut, A):
    """int_acut^inf r^4 e^{-A r^2} dr for a point acut and an interval A > 0."""
    aI = I.frac(acut)
    e = exp_neg(A * (aI * aI))
    sA = sqrt(A)
    erfc = I(ONE) - erf_pos(aI * sA)
    return (aI * aI * aI) * e / (2 * A) + 3 * aI * e / (4 * A * A) + 3 * sqrt(PI) * erfc / (8 * A * A * sA)


def realise(a1, b1, cq):
    c = max(cq, abs(a1 - b1))
    b = b1
    if c > a1 + b:
        b = c - a1
    return a1, b, c


def cos_int(x):
    """cos on an interval inside [0, pi] (decreasing), Taylor series at each end with the alternating bound."""
    def at(v, low):
        V2 = v * v
        tl = tu = ONE
        s_lo = s_hi = 0
        hist = []
        n = 0
        prev = None
        while True:
            if n % 2 == 0:
                s_lo += tl; s_hi += tu
            else:
                s_lo -= tu; s_hi -= tl
            hist.append((s_lo, s_hi))
            if prev is not None and tu < prev and tu < 4:
                break
            prev = tu
            n += 1
            den = (2 * n - 1) * (2 * n) * ONE * ONE
            tl, tu = fdiv(tl * V2, den), cdiv(tu * V2, den)
        (a_lo, a_hi), (b_lo, b_hi) = hist[-1], hist[-2]
        return min(a_lo, b_lo) if low else max(a_hi, b_hi)
    return I(at(x.hi, True), at(x.lo, False))  # cos decreasing on [0, pi]


def qlo(x):
    return Fraction(x.lo, ONE)


def qhi(x):
    return Fraction(x.hi, ONE)


TABLE = {}


def triples(acut, na, nt, nc, log=True, wpath=Fraction(1, 2)):
    acut = Fraction(acut)
    coses = [cos_int(PI * Fraction(m, nt)) for m in range(nt + 1)]
    table = TABLE.setdefault((acut, na, nc), {})
    total = Fraction(0)
    dth = PI * Fraction(1, nt)
    four_pi = PI * 4
    t0 = time.time()
    for i in range(na):
        a0, a1 = Fraction(i, na), Fraction(i + 1, na)
        for j in range(na):
            b0, b1 = Fraction(j, na), Fraction(j + 1, na)
            wab = (a1 * a1 - a0 * a0) / 2 * (b1 * b1 - b0 * b0) / 2
            corners = [(a0, b0), (a0, b1), (a1, b0), (a1, b1)]
            acc = I(0)
            for m in range(nt):
                k1lo = qlo(coses[m + 1])  # cos(theta1) >= k1lo, so c^2 <= a^2 + b^2 - 2 k1lo a b
                c2max = max(p * p + q * q - 2 * k1lo * p * q for p, q in corners)
                k0hi = qhi(coses[m])  # cos(theta0) <= k0hi, so c^2 >= a^2 + b^2 - 2 k0hi a b
                g = lambda p, q: p * p + q * q - 2 * k0hi * p * q
                cand = [g(p, q) for p, q in corners]
                for (fix_a, val) in ((True, a0), (True, a1), (False, b0), (False, b1)):
                    if fix_a:
                        q = k0hi * val
                        if b0 <= q <= b1: cand.append(g(val, q))
                    else:
                        p = k0hi * val
                        if a0 <= p <= a1: cand.append(g(p, val))
                c2min = min(cand)
                w = wpath if c2min >= 1 else Fraction(1, 6)
                X = c2max * nc * nc
                k = int(X.numerator // X.denominator)
                k = math.isqrt(k)
                while k * k < X:
                    k += 1
                k = max(k, 1)
                key = (min(i + 1, j + 1), max(i + 1, j + 1), k)
                if key not in table:
                    A_ = union3(*realise(Fraction(key[0], na), Fraction(key[1], na), Fraction(k, nc)))
                    table[key] = qlo(G3(acut, I(A_.hi)))  # G3 at the upper end of A: a lower bound
                acc = acc + I.frac(table[key] * w)
            total += qlo(four_pi * dth * acc * wab)
        if log and (i % 8 == 7 or i == na - 1):
            print(f"  rows a < {float(a1):.4f}: triple term >= {float(total):.12f} ({len(table)} table entries, "
                  f"{time.time() - t0:.0f} s)", flush=True)
    return total, len(table)


def main(argv):
    acut = Fraction(2107, 2000); na, nt, nc, pcells = 128, 256, 128, 200000
    claim = Fraction("0.6339"); controls = False
    i = 0
    while i < len(argv):
        k = argv[i]
        if k == "--controls": controls = True; i += 1; continue
        v = argv[i + 1]; i += 2
        if k == "--a": acut = Fraction(v)
        elif k == "--na": na = int(v)
        elif k == "--nt": nt = int(v)
        elif k == "--nc": nc = int(v)
        elif k == "--pcells": pcells = int(v)
        elif k == "--claim": claim = Fraction(v)
        else: raise SystemExit("unknown argument " + k)
    # self-tests of the union routine: two coincident-free configurations with known areas
    eq = union3(Fraction(1), Fraction(1), Fraction(1))  # equilateral side 1: 3pi - 3L(1) + T, T = (pi - sqrt 3)/2
    L1 = 2 * PI / 3 - sqrt(I.frac(3)) * Fraction(1, 2)
    T = (PI - sqrt(I.frac(3))) * Fraction(1, 2)
    ref = 3 * PI - 3 * L1 + T
    assert eq.lo <= ref.hi and ref.lo <= eq.hi and eq.width() < Fraction(1, 10 ** 20), (eq.lof(), ref.lof())
    col = union3(Fraction(1), Fraction(1), Fraction(2))  # collinear 0, 1, 2 (third point at -1): 3pi - 2L(1)
    ref2 = 3 * PI - 2 * L1
    assert col.lo <= ref2.hi and ref2.lo <= col.hi and col.width() < Fraction(1, 10 ** 20), (col.lof(), ref2.lof())
    # G3 against monotone Riemann sums of r^4 e^{-A r^2} on [1, 2] with A = 4 (the tail beyond 2 is below 1e-6)
    gA, ga, steps = I.frac(4), Fraction(1), 3000
    h = Fraction(1, steps)
    lo_s = hi_s = Fraction(0)
    for q in range(steps):
        r0, r1 = ga + h * q, ga + h * (q + 1)  # r^4 e^{-A r^2} decreases for r > sqrt(2/A)
        lo_s += h * Fraction((I.frac(r1 ** 4) * exp_neg(gA * I.frac(r1 * r1))).lo, ONE)
        hi_s += h * Fraction((I.frac(r0 ** 4) * exp_neg(gA * I.frac(r0 * r0))).hi, ONE)
    g = G3(ga, gA)
    assert lo_s <= g.hif() and g.lof() <= hi_s + Fraction(1, 10 ** 6), (float(lo_s), float(g.lof()), float(hi_s))
    print("self-tests: union of three unit discs, equilateral side 1 and collinear spacing 1; G3 against Riemann sums: ok")
    lo2, H, T1, T2 = LA.bound(acut, pcells)
    print(f"cut a = {acut}; head + points + pairs (pairs on {pcells} cells) >= {float(lo2):.12f}")
    T3, nt_ = triples(acut, na, nt, nc)
    lo = lo2 + T3
    print(f"triples (na {na}, nt {nt}, nc {nc}; {nt_} union areas) >= {float(T3):.12f}")
    print(f"certified: beta_MST >= {float(lo):.12f}")
    ok = lo >= claim
    print(f"claim {float(claim)}: {'HOLDS' if ok else 'NOT IMPLIED'}")
    if controls:
        bad = 0
        r = not (lo >= Fraction(int(lo * 10 ** 6) + 1, 10 ** 6)); bad += not r
        print(f"red control 1 (claim one unit above in the 6th decimal): {'refused' if r else 'ACCEPTED (bad)'}")
        T3c, _ = triples(acut, 16, 32, 16, log=False)
        r = not (lo2 + T3c >= claim); bad += not r
        print(f"red control 2 (coarse triple grid 16/32/16 gives {float(lo2 + T3c):.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        r = not (lo2 >= claim); bad += not r
        print(f"red control 3 (no triple term, {float(lo2):.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        T3w, _ = triples(acut, na, nt, nc, log=False, wpath=Fraction(1, 6))
        r = not (lo2 + T3w >= claim); bad += not r
        print(f"red control 4 (paths weighted 1/6 like triangles, {float(lo2 + T3w):.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        if bad:
            print("RED CONTROLS FAILED"); return 2
        print("red controls: 4 of 4 refused")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
