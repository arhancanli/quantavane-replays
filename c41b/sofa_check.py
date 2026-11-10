#!/usr/bin/env python3
"""Check of the 41b row C41b <= 2 sqrt(2) - 1* (teorth/optimizationproblems #358, vidigoat, 10 Nov 2026; the ambidextrous
moving sofa), from the proof in tao/certificates/c41b/README.md @ cf620ea alone (check_sofa.py was neither read nor run).

The argument (re-derived by hand). An ambidextrous sofa S lies in the strip H = R x [0, 1] and can turn both ways by
monotone rotations of total angles b1, b2 <= pi/2. If some b_i < pi/4, S also lies in a unit strip at angle 90 - b_i to H,
so area(S) <= 1/cos b_i < sqrt 2. Otherwise, at the moments when the two rotations equal pi/4, S lies in H, in the
left-turn hallway turned by -45 degrees and in the right-turn hallway turned by +45 degrees, each translated. In
s = x - y, t = x + y (areas double, r = sqrt 2) these are D = {0 <= t - s <= 2},
    L1 = {s >= a, t <= b} minus {s > a + r, t < b - r}  and  L2 = {s <= d, t >= c} minus {s < d - r, t > c + r},
and the claim is area(D n L1 n L2) <= 4 sqrt 2 - 2 for all a, b, c, d, hence area(S) <= 2 sqrt 2 - 1.

What this program checks, in exact arithmetic in Q(sqrt 2) unless marked:
  1. the reduction: the two hallways {x <= 1, 0 <= y <= 1} u {0 <= x <= 1, y >= 0} (left) and
     {x <= 1, 0 <= y <= 1} u {0 <= x <= 1, y <= 1} (right), rotated by -45 and +45 degrees and mapped to (s, t), are the
     L-shapes L2 and L1 with arms of width r (their corner points mapped exactly), and H maps to D;
  2. Lemma A (an r-square meets D in area <= 2 sqrt 2 - 1) and Lemma B (g(w) for a w x r rectangle) from the exact
     area of D n rectangle, maximised over the position along the diagonal: for each w the area is a concave piecewise
     quadratic in the offset and its maximum is attained where the two cut corners are equal, which we compare with the
     exact area at that offset and at a grid of 401 offsets; g(w) <= sqrt 2 - 1 + w on [0, r] (h = sqrt2 - 1 + w - g is a
     quadratic on [2 - r, r] with h(r) = h'(r) = 0 and h'' = 1/2, and w (r - 1) <= r - 1 below 2 - r);
  3. the case arithmetic: 2 (r - w) + 2 g(w) <= 2 r + 2 sqrt 2 - 2 = 4 sqrt 2 - 2 for w in [0, r], and 2 delta <= 2 sqrt 2
     for delta <= r, and 2 (2 sqrt 2 - 1) = 4 sqrt 2 - 2;
  4. the decomposition F = P u Q u U u V of the README, tested pointwise (exactly) at 20,000 random points of random
     configurations: every point of F is of one of the four types, and the type conditions hold as stated;
  5. sharpness in the reduced problem: area(F) = 4 sqrt 2 - 2 exactly at (a, b, c, d) = (0, 2r + 1 + u, 1, 2r + u) for
     u = 0, 1/4, 1/2, 1 (area of F by inclusion and exclusion over the arms, each term the exact area of D n rectangle);
  6. (floats, evidence only) a random search and coordinate refinement over (b, c, d) with a = 0: the largest area found
     is 4 sqrt 2 - 2 to within 2e-15 and never above it by more than rounding.
Red controls (each must fail): Lemma B's right side lowered by 1/100; the bound 4 sqrt 2 - 2 - 1/1000 at the sharp
configuration; the left hallway turned the wrong way (both L-shapes of the same type), which no longer gives L1 and L2;
and, as an anchor, a single 45-degree hallway (Hammersley's bound): max area(D n L1) = 4 sqrt 2 in (s, t), i.e. 2 sqrt 2.
Not checked: the setting itself (the Gerver and Kallus-Romik formulation of the moving sofa problem for turns of at most
90 degrees, and that both turns pass through 45 degrees by continuity), which we read.

Usage: python3 sofa_check.py [--controls]     (standard library only; under a minute)
"""
import random
import sys
import time
from fractions import Fraction as F


class Q2:
    """a + b sqrt 2 with rational a, b, exactly"""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        if isinstance(a, Q2):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(s, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(s.a + o.a, s.b + o.b)

    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(s.a - o.a, s.b - o.b)

    def __rsub__(s, o):
        return Q2(o) - s

    def __neg__(s):
        return Q2(-s.a, -s.b)

    def __mul__(s, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(s.a * o.a + 2 * s.b * o.b, s.a * o.b + s.b * o.a)

    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, Q2) else Q2(o)
        n = o.a * o.a - 2 * o.b * o.b
        return s * Q2(o.a / n, -o.b / n)

    def sign(s):
        """exact sign of a + b sqrt 2"""
        if s.b == 0:
            return (s.a > 0) - (s.a < 0)
        if s.a == 0:
            return (s.b > 0) - (s.b < 0)
        if (s.a > 0) == (s.b > 0):
            return 1 if s.a > 0 else -1
        # opposite signs: compare a^2 with 2 b^2
        d = s.a * s.a - 2 * s.b * s.b
        return (1 if s.a > 0 else -1) if d > 0 else (-1 if s.a > 0 else 1) if d < 0 else 0

    def __lt__(s, o):
        return (s - o).sign() < 0

    def __le__(s, o):
        return (s - o).sign() <= 0

    def __gt__(s, o):
        return (s - o).sign() > 0

    def __ge__(s, o):
        return (s - o).sign() >= 0

    def __eq__(s, o):
        return (s - o).sign() == 0

    def __hash__(s):
        return hash((s.a, s.b))

    def __float__(s):
        return float(s.a) + float(s.b) * 2 ** 0.5

    def __repr__(s):
        return "%s + %s sqrt2" % (s.a, s.b)


R = Q2(0, 1)          # sqrt 2
BIG = Q2(10 ** 4)


def mx(x, y):
    return x if x >= y else y


def mn(x, y):
    return x if x <= y else y


def area_D_rect(s0, s1, t0, t1):
    """exact area of D n [s0, s1] x [t0, t1], D = {0 <= t - s <= 2}: the integral over s of
    max(0, min(t1, s + 2) - max(t0, s)), linear between the breakpoints t0 - 2, t0, t1 - 2, t1"""
    s0, s1 = mx(s0, t0 - 2), mn(s1, t1)
    if not (s1 > s0) or not (t1 > t0):
        return Q2(0)
    pts = sorted({s0, s1} | {x for x in (t0 - 2, t0, t1 - 2, t1) if s0 < x < s1}, key=float)
    # floats only order the points; the order is then confirmed exactly
    for u, v in zip(pts, pts[1:]):
        assert u < v
    f = lambda s: mx(Q2(0), mn(t1, s + 2) - mx(t0, s))
    tot = Q2(0)
    for u, v in zip(pts, pts[1:]):
        tot = tot + (v - u) * (f(u) + f(v)) / 2
    return tot


def arms(a, b, c, d):
    V1 = (a, a + R, -BIG, b)
    H1 = (a, BIG, b - R, b)
    C1 = (a, a + R, b - R, b)
    V2 = (d - R, d, c, BIG)
    H2 = (-BIG, d, c, c + R)
    C2 = (d - R, d, c, c + R)
    return [(V1, 1), (H1, 1), (C1, -1)], [(V2, 1), (H2, 1), (C2, -1)]


def area_F(a, b, c, d):
    """area(D n L1 n L2) by inclusion and exclusion: 1_L = 1_V + 1_H - 1_C for each L-shape"""
    A1, A2 = arms(a, b, c, d)
    tot = Q2(0)
    for (p, sp) in A1:
        for (q, sq) in A2:
            s0, s1 = mx(p[0], q[0]), mn(p[1], q[1])
            t0, t1 = mx(p[2], q[2]), mn(p[3], q[3])
            tot = tot + (sp * sq) * area_D_rect(s0, s1, t0, t1)
    return tot


def area_L1(a, b):
    A1, _ = arms(a, b, Q2(0), Q2(0))
    tot = Q2(0)
    for (p, sp) in A1:
        tot = tot + sp * area_D_rect(*p)
    return tot


# ------------------------------------------------------------------------------- 1. the reduction
def reduction_check(wrong_way=False):
    """corner points of the hallways, rotated and mapped to (s, t) = (x - y, x + y), exactly.
    Rotation by -45 degrees: (x, y) -> ((x + y)/r, (y - x)/r); by +45: (x, y) -> ((x - y)/r, (x + y)/r)."""
    def rot(p, sign):
        x, y = p
        if sign < 0:
            X, Y = (x + y) / R, (y - x) / R
        else:
            X, Y = (x - y) / R, (x + y) / R
        return (X - Y, X + Y)
    left = [(Q2(1), Q2(0)), (Q2(1), Q2(1)), (Q2(0), Q2(1)), (Q2(0), Q2(0))]     # corners of the left hallway's elbow
    right = [(Q2(1), Q2(1)), (Q2(1), Q2(0)), (Q2(0), Q2(0)), (Q2(0), Q2(1))]
    L = [rot(p, -1 if not wrong_way else +1) for p in left]
    Rr = [rot(p, +1) for p in right]
    # left: (s, t) = r (x, y), so the elbow is the square [0, r]^2, the horizontal arm x <= 1 becomes s <= r with
    # 0 <= t <= r, and the vertical arm y >= 0 becomes t >= 0 with 0 <= s <= r: the shape L2 with d = r, c = 0
    okL = sorted(L, key=lambda q: (float(q[0]), float(q[1]))) == sorted(
        [(R, Q2(0)), (R, R), (Q2(0), R), (Q2(0), Q2(0))], key=lambda q: (float(q[0]), float(q[1])))
    okL &= all(rot((Q2(-5), Q2(0)), -1)[k] == (Q2(-5) * R, Q2(0))[k] for k in (0, 1)) if not wrong_way else True
    # right: (s, t) = (-r y, r x): the elbow is [-r, 0] x [0, r]; arm x <= 1 -> t <= r with -r <= s <= 0 (down),
    # arm y <= 1 -> s >= -r with 0 <= t <= r (right): the shape L1 with a = -r, b = r
    okR = sorted(Rr, key=lambda q: (float(q[0]), float(q[1]))) == sorted(
        [(-R, R), (Q2(0), R), (-R, Q2(0)), (Q2(0), Q2(0))], key=lambda q: (float(q[0]), float(q[1])))
    # sample points: hallway membership equals L-shape membership of the image (exact), 4,000 points each
    rng = random.Random(1945)
    inLeft = lambda x, y: (x <= 1 and Q2(0) <= y <= Q2(1)) or (Q2(0) <= x <= Q2(1) and y >= 0)
    inRight = lambda x, y: (x <= 1 and Q2(0) <= y <= Q2(1)) or (Q2(0) <= x <= Q2(1) and y <= 1)
    inL2 = lambda s_, t_: (s_ <= R and t_ >= 0) and not (s_ < 0 and t_ > R)          # d = r, c = 0
    inL1 = lambda s_, t_: (s_ >= -R and t_ <= R) and not (s_ > 0 and t_ < 0)          # a = -r, b = r
    for _ in range(4000):
        x, y = Q2(F(rng.randint(-400, 400), 100)), Q2(F(rng.randint(-400, 400), 100))
        sL, tL = rot((x, y), -1 if not wrong_way else +1)
        okL &= inLeft(x, y) == inL2(sL, tL)
        sR, tR = rot((x, y), +1)
        okR &= inRight(x, y) == inL1(sR, tR)
    # H = {0 <= y <= 1}: y = (t - s)/2, so 0 <= t - s <= 2 (D)
    return okL, okR


# ------------------------------------------------------------------------------- 2, 3. lemmas and cases
def rect_max(w, n=400):
    """max over the offset of area(D n [0, w] x [o, o + r]) (a w x r rectangle slid along the diagonal; by the
    translation (s, t) -> (s + u, t + u) every position is one of these): exact at n + 1 offsets and at the balanced
    offset where the two cut corners are equal."""
    best = Q2(0)
    lo, hi = -R, Q2(2) + w            # outside [-r, 2 + w] the rectangle misses the band
    for k in range(n + 1):
        o = lo + (hi - lo) * F(k, n)
        A = area_D_rect(Q2(0), Q2(w), o, o + R)
        best = mx(best, A)
    # balanced offset: the rectangle's t - s range [o - w, o + r] centred in [0, 2]: o - w + o + r = 2
    ob = (Q2(2) + w - R) / 2
    Ab = area_D_rect(Q2(0), Q2(w), ob, ob + R)
    return mx(best, Ab), Ab


def g_formula(w):
    w = Q2(w) if not isinstance(w, Q2) else w
    if w <= 2 - R:
        return w * R
    return w * R - (w + R - 2) * (w + R - 2) / 4


def lemmas(shift=Q2(0)):
    ok = True
    worst = None
    for k in range(0, 41):
        w = R * F(k, 40)
        M, Ab = rect_max(w)
        g = g_formula(w)
        ok &= (M == Ab) and (Ab == g)                      # the balanced offset is the maximum and equals g(w)
        ok &= g <= R - 1 + w - shift                        # Lemma B's bound
        gap = R - 1 + w - g
        worst = gap if worst is None or gap < worst else worst
    # h(w) = sqrt2 - 1 + w - g(w) on [2 - r, r]: h(r) = 0, h'(r) = 1 - r + (2r - 2)/2 = 0, h'' = 1/2
    h_r = R - 1 + R - g_formula(R)
    dh_r = 1 - R + (R + R - 2) / 2
    lemA = g_formula(R) == 2 * R - 1
    return ok and h_r == 0 and dh_r == 0 and lemA, worst


def cases():
    ok = True
    for k in range(0, 101):
        w = R * F(k, 100)
        ok &= 2 * (R - w) + 2 * g_formula(w) <= 2 * R + 2 * R - 2          # case 3, delta = r + w in [r, 2r]
    ok &= 2 * R <= 4 * R - 2                                                   # case 3, delta <= r: 2 delta <= 2 sqrt 2
    ok &= 2 * (2 * R - 1) == 4 * R - 2                                         # cases 1 and 2
    return ok


# ------------------------------------------------------------------------------- 4. decomposition
def decomposition(n=20000, seed=41, max_tries=2000000):
    """n random exact points of F (over random configurations with a = 0), each tested against the four types"""
    rng = random.Random(seed)
    bad = hits = tries = 0
    while hits < n and tries < max_tries:
        tries += 1
        a = Q2(0)
        b, c, d = Q2(F(rng.randint(150, 600), 100)), Q2(F(rng.randint(-150, 300), 100)), Q2(F(rng.randint(0, 600), 100))
        s = Q2(F(rng.randint(-100, 700), 100))
        t = s + Q2(F(rng.randint(0, 200), 100))          # inside D
        inL1 = (s >= a and t <= b) and not (s > a + R and t < b - R)
        inL2 = (s <= d and t >= c) and not (s < d - R and t > c + R)
        if not (inL1 and inL2):
            continue
        hits += 1
        Pt = s <= a + R and t <= c + R
        Qt = t >= b - R and s >= d - R
        Ut = (d - R <= s <= a + R) and (c + R < t < b - R)
        Vt = (a + R < s < d - R) and (b - R <= t <= c + R)
        delta, eta = d - a, b - c
        if not (Pt or Qt or Ut or Vt):
            bad += 1
        if Ut and not (delta <= 2 * R and eta > 2 * R):
            bad += 1
        if Vt and not (delta > 2 * R and eta <= 2 * R):
            bad += 1
    return bad, hits


# ------------------------------------------------------------------------------- 6. floats
def area_float(b, c, d, a=0.0):
    """the same area in binary64 (for the search only)"""
    r, BIGf = 2 ** 0.5, 1e4

    def rect(s0, s1, t0, t1):
        s0, s1 = max(s0, t0 - 2), min(s1, t1)
        if s1 <= s0 or t1 <= t0:
            return 0.0
        pts = sorted({s0, s1} | {x for x in (t0 - 2, t0, t1 - 2, t1) if s0 < x < s1})
        f = lambda x: max(0.0, min(t1, x + 2) - max(t0, x))
        return sum((v - u) * (f(u) + f(v)) / 2 for u, v in zip(pts, pts[1:]))
    A1 = [((a, a + r, -BIGf, b), 1), ((a, BIGf, b - r, b), 1), ((a, a + r, b - r, b), -1)]
    A2 = [((d - r, d, c, BIGf), 1), ((-BIGf, d, c, c + r), 1), ((d - r, d, c, c + r), -1)]
    tot = 0.0
    for p, sp in A1:
        for q, sq in A2:
            tot += sp * sq * rect(max(p[0], q[0]), min(p[1], q[1]), max(p[2], q[2]), min(p[3], q[3]))
    return tot


def float_search(n=200000, seed=358):
    r = 2 ** 0.5
    rng = random.Random(seed)
    best = (0.0, None)
    for _ in range(n):
        b, c, d = rng.uniform(-1, 7), rng.uniform(-3, 5), rng.uniform(-2, 7)
        A = area_float(b, c, d)
        if A > best[0]:
            best = (A, (b, c, d))
    A0, (b, c, d) = best
    step = 0.05
    while step > 1e-10:
        imp = False
        for db, dc, dd in ((step, 0, 0), (-step, 0, 0), (0, step, 0), (0, -step, 0), (0, 0, step), (0, 0, -step)):
            A = area_float(b + db, c + dc, d + dd)
            if A > A0:
                A0, b, c, d, imp = A, b + db, c + dc, d + dd, True
        if not imp:
            step /= 2
    return A0, (b, c, d), 4 * r - 2


def main():
    t0 = time.time()
    controls = "--controls" in sys.argv
    print("41b: C41b <= 2 sqrt 2 - 1 (#358), the reduced problem area(D n L1 n L2) <= 4 sqrt 2 - 2 in (s, t)")
    okL, okR = reduction_check()
    print("  1. reduction: the left hallway turned by -45 degrees is L2 (corner square [0, r]^2): %s; the right hallway "
          "turned by +45 degrees is L1 (corner square [-r, 0] x [0, r]): %s" % (okL, okR))
    okLem, worst = lemmas()
    print("  2. Lemma B: at 41 widths w in [0, r] the largest area of D n (w x r rectangle) over 401 exact offsets and "
          "the balanced offset is g(w), and g(w) <= sqrt 2 - 1 + w (smallest gap %.6f, 0 at w = r); h(r) = h'(r) = 0; "
          "Lemma A g(r) = 2 sqrt 2 - 1: %s" % (float(worst), okLem))
    okC = cases()
    print("  3. cases: 2(r - w) + 2 g(w) <= 4 sqrt 2 - 2 on [0, r] (101 points, exact; it is the quadratic h >= 0 "
          "doubled), 2 delta <= 2 sqrt 2 for delta <= r, 2 (2 sqrt 2 - 1) = 4 sqrt 2 - 2: %s" % okC)
    bad, hits = decomposition()
    print("  4. F = P u Q u U u V with the README's conditions on delta and eta, at %d random exact points of F: %d "
          "exceptions" % (hits, bad))
    sharp = []
    for u in (F(0), F(1, 4), F(1, 2), F(1)):
        A = area_F(Q2(0), 2 * R + 1 + u, Q2(1), 2 * R + u)
        sharp.append(A == 4 * R - 2)
    print("  5. area(F) = 4 sqrt 2 - 2 exactly at (a, b, c, d) = (0, 2r + 1 + u, 1, 2r + u), u = 0, 1/4, 1/2, 1: %s"
          % sharp)
    A0, cfg, target = float_search()
    print("  6. (floats) random search and refinement over (b, c, d), a = 0: largest area %.16f at (%.6f, %.6f, %.6f); "
          "4 sqrt 2 - 2 = %.16f" % (A0, cfg[0], cfg[1], cfg[2], target))
    ok = okL and okR and okLem and okC and bad == 0 and all(sharp) and A0 <= target + 4e-15
    print("  so area(S) <= (4 sqrt 2 - 2)/2 = 2 sqrt 2 - 1 = %.10f on the reduction, and the bound is sharp for it: %s"
          % (float(2 * R - 1), ok))
    if controls:
        print("Red controls (each must fail):")
        c1, _ = lemmas(shift=Q2(F(1, 100)))
        print("  1. Lemma B with its right side lowered by 1/100: %s" % ("fails" if not c1 else "HOLDS (bad)"))
        A = area_F(Q2(0), 2 * R + 1, Q2(1), 2 * R)
        c2 = not (A <= 4 * R - 2 - Q2(F(1, 1000)))
        print("  2. the bound 4 sqrt 2 - 2 - 1/1000 at the sharp configuration: %s" % ("fails" if c2 else "HOLDS (bad)"))
        okLw, _ = reduction_check(wrong_way=True)
        print("  3. the left hallway turned by +45 degrees (the wrong way): it is not L2: %s" % (not okLw))
        best = Q2(0)
        for k in range(0, 81):
            best = mx(best, area_L1(Q2(0), Q2(F(k, 10))))
        c4 = best == 4 * R
        print("  4. anchor: one 45-degree hallway alone, max over b of area(D n L1) = %s = 4 sqrt 2 (Hammersley's "
              "2 sqrt 2 in the plane): %s" % (best, c4))
        ok &= (not c1) and c2 and (not okLw) and c4
    print("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("run time %.0f s" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
