#!/usr/bin/env python3
"""Exact replay of the 73a lower-bound row [CS2019] (teorth/optimizationproblems, main 87b4ebb): Flt(3) >= 2 + sqrt 2,
"achieved by an explicit hollow non-lattice tetrahedron in dimension 3" (Codenotti and Santos, arXiv 1812.00916,
section 3): with the affine lattice L1 = {(a,b,c): a,b,c in 1+2Z, a+b+c in 1+4Z} (a translate of the fcc lattice
L = {a,b,c in 2Z, a+b+c in 4Z}, whose dual is L* = {u in (1/4)Z^3 : u1+u2, u1+u3, u2+u3 in (1/2)Z}), the tetrahedron
    A = (2+r, r, 2+r),  B = (-r, 2+r, -2-r),  C = (-2-r, -r, 2+r),  D = (r, -2-r, -2-r),   r = sqrt 2.

Exact arithmetic in Q(sqrt 2) (pairs of rationals; signs decided by comparing squares).  Checks:
  1. hollow: no point of L1 lies in the open tetrahedron (every point of L1 in a box containing it is tested
     against the four facet inequalities, strictly);
  2. width: for every nonzero u in L* with |u|^2 <= 9 the width max_v u.v - min_v u.v over the vertices is computed;
     the minimum is 2 + sqrt 2.  Every other u has larger width: the tetrahedron contains the regular tetrahedron
     T = conv{(-1,1,1), (-1,-1,-1), (1,-1,1), (1,1,-1)} (checked: its vertices satisfy the facet inequalities), whose
     inscribed ball has centre 0 and radius 1/sqrt 3, so the width in direction u is at least 2|u|/sqrt 3, which
     exceeds 2 + sqrt 2 once |u|^2 > (3/4)(2 + sqrt 2)^2 = 4.5 + 3 sqrt 2 (< 9).
So w(Delta) = 2 + sqrt 2 for a hollow convex body and Flt(3) >= 2 + sqrt 2.  Written from the paper's definitions
(read, not run: the paper has no program); no download.
usage: python3 -I cs2019_check.py [--controls]
"""
import itertools
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


class Q2:
    """a + b sqrt 2 with rational a, b."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a, self.b = F(a), F(b)

    def __add__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Q2) else Q2(o)))

    def __rsub__(self, o):
        return Q2(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.a * o.a + 2 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def sign(self):
        a, b = self.a, self.b
        if a >= 0 and b >= 0:
            return 0 if a == 0 and b == 0 else 1
        if a <= 0 and b <= 0:
            return -1
        # opposite signs: compare a^2 with 2 b^2
        d = a * a - 2 * b * b
        return (1 if a > 0 else -1) if d > 0 else (-1 if a > 0 else 1)

    def __eq__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return self.a == o.a and self.b == o.b

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __le__(self, o):
        return (self - o).sign() <= 0

    def __hash__(self):
        return hash((self.a, self.b))

    def __float__(self):
        return float(self.a) + float(self.b) * 2 ** 0.5

    def __repr__(self):
        return f"{self.a} + {self.b}*sqrt2"


R = Q2(0, 1)
TWO_R = 2 + R


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def sub(u, v):
    return tuple(a - b for a, b in zip(u, v))


def facets(verts):
    """For each facet: (normal n, offset c) with n.x <= c on the tetrahedron (n.x = c on the facet)."""
    out = []
    for i in range(4):
        p, q, s = [verts[j] for j in range(4) if j != i]
        n = cross(sub(q, p), sub(s, p))
        c = dot(n, p)
        if (dot(n, verts[i]) - c).sign() > 0:
            n, c = tuple(-x for x in n), -c
        need((dot(n, verts[i]) - c).sign() < 0, "degenerate tetrahedron")
        out.append((n, c))
    return out


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def interior(fs, x):
    return all((dot(n, x) - c).sign() < 0 for n, c in fs)


def inside(fs, x):
    return all((dot(n, x) - c).sign() <= 0 for n, c in fs)


def width(verts, u):
    vals = [dot(u, v) for v in verts]
    lo = min(vals, key=lambda t: t)
    hi = max(vals, key=lambda t: t)
    return hi - lo


def dual_vectors(maxsq):
    """Nonzero u in L* with |u|^2 <= maxsq, one of each pair +-u."""
    out = []
    m = int(4 * maxsq ** 0.5) + 1
    for a, b, c in itertools.product(range(-m, m + 1), repeat=3):
        if (a, b, c) == (0, 0, 0) or (a, b, c) < (0, 0, 0):
            continue
        if (a + b) % 2 or (a + c) % 2 or (b + c) % 2:   # u = (a,b,c)/4: u1+u2 in Z/2 iff a+b even
            continue
        if F(a * a + b * b + c * c, 16) <= maxsq:
            out.append((F(a, 4), F(b, 4), F(c, 4)))
    return out


def analyse(verts, label, verbose=True):
    fs = facets(verts)
    # 1. hollow: lattice points of L1 in a box around the tetrahedron
    box = 1 + max(int(abs(float(c))) for v in verts for c in v)
    odd = [k for k in range(-box - 1, box + 2) if k % 2]
    pts = [(a, b, c) for a in odd for b in odd for c in odd if (a + b + c - 1) % 4 == 0]
    inner = [p for p in pts if interior(fs, tuple(Q2(t) for t in p))]
    onbd = [p for p in pts if inside(fs, tuple(Q2(t) for t in p)) and p not in inner]
    # T inside
    T = [(-1, 1, 1), (-1, -1, -1), (1, -1, 1), (1, 1, -1)]
    t_in = all(inside(fs, tuple(Q2(t) for t in v)) for v in T)
    # 2. widths
    us = dual_vectors(9)
    ws = [(width(verts, tuple(Q2(t) for t in u)), u) for u in us]
    wmin = min(w for w, _ in ws)
    argmin = [u for w, u in ws if w == wmin]
    if verbose:
        print(f"{label}: {len(pts)} points of L1 in the box |coordinates| <= {box + 1}; in the open tetrahedron: "
              f"{len(inner)}; on its boundary: {len(onbd)}; T inside: {t_in}")
        print(f"  {len(us)} pairs +-u in L* with |u|^2 <= 9; least width {wmin} = {float(wmin):.12f}, attained by "
              f"{len(argmin)} of them: {[tuple(str(x) for x in u) for u in argmin]}")
    return inner, t_in, wmin, argmin


def main():
    verts = [(TWO_R, R, TWO_R), (-R, TWO_R, -TWO_R), (-TWO_R, -R, TWO_R), (R, -TWO_R, -TWO_R)]
    inner, t_in, wmin, argmin = analyse(verts, "Delta (the paper's vertices)")
    bound_ok = (Q2(F(9, 2), 3) <= Q2(9)) and (Q2(F(3, 4)) * TWO_R * TWO_R == Q2(F(9, 2), 3))
    good = not inner and t_in and wmin == TWO_R and len(argmin) == 7 and bound_ok
    print(f"  the tail bound: (3/4)(2 + sqrt 2)^2 = 9/2 + 3 sqrt 2 <= 9: {bound_ok}, so no u outside the search has "
          f"width <= 2 + sqrt 2")
    print(f"RESULT [CS2019]: Delta is hollow with respect to L1 and has lattice width exactly 2 + sqrt 2 = "
          f"{float(TWO_R):.12f} (seven functionals up to sign, as the paper says), so Flt(3) >= 2 + sqrt 2: "
          f"{'HOLDS' if good else 'FAILS'}")
    ok = good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        # vertex A pushed outward by 1/10 in its first coordinate: then a lattice point enters, or the width changes
        bad = [(TWO_R + Q2(F(1, 2)), R, TWO_R)] + verts[1:]
        inner2, _, w2, _ = analyse(bad, "  A moved to (5/2 + sqrt 2, sqrt 2, 2 + sqrt 2)")
        res.append(("A moved outward by 1/2", bool(inner2) or not (w2 == TWO_R)))
        # the family's neighbour (x, y) = (2 + sqrt 2, sqrt 2) with z = 3 instead of 2 + sqrt 2 (off the family)
        bad2 = [(TWO_R, R, Q2(3)), (-R, TWO_R, Q2(-3)), (-TWO_R, -R, Q2(3)), (R, -TWO_R, Q2(-3))]
        inner3, t3, w3, _ = analyse(bad2, "  z = 3 instead of 2 + sqrt 2")
        res.append(("z = 3", bool(inner3) or not t3 or w3 != TWO_R))
        res.append(("a claimed width 2 + sqrt 2 + 1/1000", wmin != TWO_R + Q2(F(1, 1000))))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
