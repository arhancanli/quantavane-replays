#!/usr/bin/env python3
"""MST constant in the plane: certified lower bound with isolated points, pairs and TRIPLES (program B).

Written separately from lower_triples_a.py and sharing no code with it; it extends lower_pairs_b.py (its pi, erf and
the head, point and pair terms). Decimal arithmetic. Same mathematics:
  beta >= int_0^a chi + int_a^inf (n1 + n2) + 4 pi int_0^1 int_0^1 int_0^pi w(c) G3(A(a,b,c)) a b dtheta db da,
w = 1/2 if c >= 1 (a path whose middle vertex is the one at 0), 1/6 if c < 1 (a triangle), G3(A) = int_acut^inf
r^4 e^{-A r^2} dr, A(a,b,c) the area of the union of three unit discs with mutual distances a, b, c.
Here the union area is computed by inclusion-exclusion, A = 3 pi - L(a) - L(b) - L(c) + T, with the lens area
L(d) = 2 acos(d/2) - (d/2) sqrt(4 - d^2) and the triple intersection T by Green's theorem over the arcs of each
circle that lie in both other discs; arc angles from atan2, their sines and cosines from Taylor series.
Monotonicity (planar Kneser-Poulsen, Bezdek and Connelly 2002): A(a,b,c) <= A(a1, b', c') for every expanded
triangle with a1 >= a, b' >= b, c' >= c; G3 decreases in A.
Rounding: the triple part runs at 40 significant digits; the union area is raised by 1e-30, each G3 value and each
cell contribution lowered by 1e-30 (the rounding error of a cell is below 1e-35; the area is a Lipschitz function of
the arc end points, so near-ties in angle comparisons move it by no more than the rounding error).

Usage: python3 -I lower_triples_b.py [--a 2107/2000] [--na 120] [--nt 300] [--nc 160] [--pcells 300000]
                                     [--claim 0.6339] [--controls]
"""
import sys, os, time, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # lower_pairs_b.py sits next to this file
from decimal import Decimal as D, getcontext
import lower_pairs_b as LB

PREC_T = 40
EPS = D("1e-30")


def dec(q):
    q = Fraction(q)
    return D(q.numerator) / D(q.denominator)


def atan(x):
    """atan(x), x >= 0: x > 1 via pi/2 - atan(1/x); two halvings x -> x/(1 + sqrt(1 + x^2)); then the series."""
    if x > 1:
        return LB.PI / 2 - atan(1 / x)
    for _ in range(2):
        x = x / (1 + (1 + x * x).sqrt())
    s = D(0); t = x; n = 0; x2 = x * x
    while True:
        term = t / (2 * n + 1)
        s += term if n % 2 == 0 else -term
        if term < D("1e-45"):
            break
        n += 1; t *= x2
    return 4 * s


def atan2(y, x):
    """Angle of (x, y) in [0, 2 pi)."""
    pi = LB.PI
    if x == 0 and y == 0:
        return D(0)
    if x > 0 and y >= 0:
        return atan(y / x)
    if x <= 0 and y > 0:
        return pi / 2 + atan(-x / y)
    if x < 0 and y <= 0:
        return pi + atan(y / x)
    return 3 * pi / 2 + atan(x / -y)


def acos(t):
    if t <= -1:
        return LB.PI
    if t >= 1:
        return D(0)
    return 2 * atan(((1 - t) / (1 + t)).sqrt())


def sincos(x):
    """sin and cos of any angle: reduce to [-pi, pi], then Taylor."""
    pi = LB.PI
    two = 2 * pi
    x = x - two * (x / two).to_integral_value(rounding="ROUND_FLOOR")
    if x > pi:
        x -= two
    s = D(0); c = D(0); t = D(1); n = 0
    while True:
        r = n % 4
        if r == 0: c += t
        elif r == 1: s += t
        elif r == 2: c -= t
        else: s -= t
        n += 1
        t = t * x / n
        if n > 8 and abs(t) < D("1e-45"):
            break
    return s, c


def lens(d):
    if d >= 2:
        return D(0)
    return 2 * acos(d / 2) - d / 2 * (4 - d * d).sqrt()


def union_ie(a, b, c):
    """Union of unit discs at P1 = (0,0), P2 = (a,0), P3 with |P3| = b, |P3 - P2| = c (a, b, c exact rationals,
    a valid possibly degenerate triangle): 3 pi - L(a) - L(b) - L(c) + T."""
    x = (a * a + b * b - c * c) / (2 * a)
    y2 = b * b - x * x
    X = dec(x); Y = dec(y2).sqrt() if y2 > 0 else D(0)
    P = [(D(0), D(0)), (dec(a), D(0)), (X, Y)]
    dist = {(0, 1): dec(a), (0, 2): dec(b), (1, 2): dec(c)}
    dd = lambda i, j: dist[(min(i, j), max(i, j))]
    T = D(0)
    for i in range(3):
        j, k = [m for m in range(3) if m != i]
        if dd(i, j) >= 2 or dd(i, k) >= 2:
            continue
        cx, cy = P[i]
        pj = atan2(P[j][1] - cy, P[j][0] - cx); aj = acos(dd(i, j) / 2)
        pk = atan2(P[k][1] - cy, P[k][0] - cx); ak = acos(dd(i, k) / 2)
        dl = pk - pj
        while dl > LB.PI: dl -= 2 * LB.PI
        while dl <= -LB.PI: dl += 2 * LB.PI
        s = max(pj - aj, pj + dl - ak); t = min(pj + aj, pj + dl + ak)
        if t > s:
            ss, cs = sincos(s); st, ct = sincos(t)
            T += ((t - s) + cx * (st - ss) - cy * (ct - cs)) / 2
    return 3 * LB.PI - lens(dist[(0, 1)]) - lens(dist[(0, 2)]) - lens(dist[(1, 2)]) + T


def G3(acut, A):
    a = dec(acut)
    e = (-(A * a * a)).exp()
    return a ** 3 * e / (2 * A) + 3 * a * e / (4 * A * A) + 3 * LB.SQPI * (1 - LB.erf(a * A.sqrt())) / (8 * A * A * A.sqrt())


def cos_bounds(m, nt):
    """Rational bounds lo <= cos(m pi / nt) <= hi."""
    s, c = LB.sin_cos(LB.PI * m / nt) if m * 4 <= nt * 4 else (None, None)
    v = Fraction(str(c))
    return v - Fraction(1, 10 ** 32), v + Fraction(1, 10 ** 32)


TABLES = {}


def triple_term(acut, na, nt, nc, wpath=Fraction(1, 2), log=True):
    acut = Fraction(acut)
    cb = [cos_bounds(m, nt) for m in range(nt + 1)]
    getcontext().prec = PREC_T
    tab = TABLES.setdefault((acut, na, nc), {})
    total = D(0)
    pi = LB.PI
    t0 = time.time()
    for i in range(na):
        a0, a1 = Fraction(i, na), Fraction(i + 1, na)
        for j in range(na):
            b0, b1 = Fraction(j, na), Fraction(j + 1, na)
            box = [(a0, b0), (a0, b1), (a1, b0), (a1, b1)]
            cell_sum = D(0)
            for m in range(nt):
                lo1 = cb[m + 1][0]   # cos(theta) >= lo1 on the cell
                hi0 = cb[m][1]       # cos(theta) <= hi0 on the cell
                c2hi = max(p * p + q * q - 2 * lo1 * p * q for p, q in box)
                f = lambda p, q: p * p + q * q - 2 * hi0 * p * q
                vals = [f(p, q) for p, q in box]
                for av in (a0, a1):
                    q = hi0 * av
                    if b0 <= q <= b1: vals.append(f(av, q))
                for bv in (b0, b1):
                    p = hi0 * bv
                    if a0 <= p <= a1: vals.append(f(p, bv))
                c2lo = min(vals)
                w = wpath if c2lo >= 1 else Fraction(1, 6)
                # smallest k with (k/nc)^2 >= c2hi
                Xn = c2hi * nc * nc
                k = math.isqrt(Xn.numerator // Xn.denominator)
                while Fraction(k * k) < Xn:
                    k += 1
                k = max(k, 1)
                aa, bb = (i + 1, j + 1) if i <= j else (j + 1, i + 1)
                key = (aa, bb, k)
                g = tab.get(key)
                if g is None:
                    A_, B_, C_ = Fraction(aa, na), Fraction(bb, na), Fraction(k, nc)
                    C_ = max(C_, abs(A_ - B_))
                    if C_ > A_ + B_:
                        B_ = C_ - A_
                    U = union_ie(A_, B_, C_) + EPS
                    g = G3(acut, U) - EPS
                    tab[key] = g
                cell_sum += g * dec(w)
            wab = dec((a1 * a1 - a0 * a0) / 2 * (b1 * b1 - b0 * b0) / 2)
            total += 4 * pi * (pi / nt) * wab * cell_sum - EPS
        if log and (i % 10 == 9 or i == na - 1):
            print(f"  a < {float(a1):.4f}: triple term >= {total:.12f} ({len(tab)} union areas, {time.time() - t0:.0f} s)",
                  flush=True)
    getcontext().prec = 60
    return total, len(tab)


def main(argv):
    acut = "2107/2000"; na, nt, nc, pcells = 120, 300, 160, 300000
    claim = D("0.6339"); controls = False
    i = 0
    while i < len(argv):
        k = argv[i]
        if k == "--controls": controls = True; i += 1; continue
        v = argv[i + 1]; i += 2
        if k == "--a": acut = v
        elif k == "--na": na = int(v)
        elif k == "--nt": nt = int(v)
        elif k == "--nc": nc = int(v)
        elif k == "--pcells": pcells = int(v)
        elif k == "--claim": claim = D(v)
        else: raise SystemExit("unknown argument " + k)
    getcontext().prec = PREC_T
    # self-tests: equilateral side 1 (T = Reuleaux triangle (pi - sqrt 3)/2), collinear spacing 1, atan and acos
    pi = LB.PI
    L1 = 2 * pi / 3 - D(3).sqrt() / 2
    assert abs(union_ie(Fraction(1), Fraction(1), Fraction(1)) - (3 * pi - 3 * L1 + (pi - D(3).sqrt()) / 2)) < D("1e-35")
    assert abs(union_ie(Fraction(1), Fraction(1), Fraction(2)) - (3 * pi - 2 * L1)) < D("1e-35")
    assert abs(4 * atan(D(1)) - pi) < D("1e-35") and abs(acos(D("0.5")) - pi / 3) < D("1e-35")
    assert abs(atan2(D(-1), D(-1)) - 5 * pi / 4) < D("1e-35")
    # G3 against the closed form's own derivative check: d/dacut G3 = -acut^4 e^{-A acut^2}
    A5 = D(5); h = Fraction(1, 10 ** 8); a = Fraction(2107, 2000)
    num = (G3(a + h, A5) - G3(a - h, A5)) / dec(2 * h)
    exact = -(dec(a) ** 4) * (-(A5 * dec(a) ** 2)).exp()
    assert abs(num - exact) < D("1e-12") * abs(exact), (num, exact)
    getcontext().prec = 60
    print("self-tests: union of three discs (equilateral side 1, collinear), atan, acos, atan2, G3' = -a^4 e^{-Aa^2}: ok")
    lo2, head, t1, t2 = LB.bound(acut, pcells)
    print(f"cut a = {acut}; head + points + pairs (pairs on {pcells} cells) >= {lo2:.12f}")
    T3, nu = triple_term(acut, na, nt, nc)
    lo = lo2 + T3
    print(f"triples (na {na}, nt {nt}, nc {nc}; {nu} union areas) >= {T3:.12f}")
    print(f"certified: beta_MST >= {lo:.12f}")
    ok = lo >= claim
    print(f"claim {claim}: {'HOLDS' if ok else 'NOT IMPLIED'}")
    if controls:
        bad = 0
        c1 = (lo * 10 ** 6).to_integral_value(rounding="ROUND_FLOOR") / 10 ** 6 + D("0.000001")
        r = not (lo >= c1); bad += not r
        print(f"red control 1 (claim {c1}, one unit above in the 6th decimal): {'refused' if r else 'ACCEPTED (bad)'}")
        T3c, _ = triple_term(acut, 12, 30, 12, log=False)
        r = not (lo2 + T3c >= claim); bad += not r
        print(f"red control 2 (coarse grid 12/30/12 gives {lo2 + T3c:.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        r = not (lo2 >= claim); bad += not r
        print(f"red control 3 (no triple term, {lo2:.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        T3w, _ = triple_term(acut, na, nt, nc, wpath=Fraction(1, 6), log=False)
        r = not (lo2 + T3w >= claim); bad += not r
        print(f"red control 4 (paths weighted 1/6, {lo2 + T3w:.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        if bad:
            print("RED CONTROLS FAILED"); return 2
        print("red controls: 4 of 4 refused")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
