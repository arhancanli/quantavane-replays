#!/usr/bin/env python3
"""Exact replay of the 43a upper-bound row [GP1968] (teorth/optimizationproblems, main 87b4ebb; added in PR #28):
C_43 <= sqrt(3)/2 ~ 0.866 from the equilateral triangle (C_43 = inf over finite planar sets of L_S/L_M, Steiner minimal
tree over minimum spanning tree).

For the unit equilateral triangle A = (0,0), B = (1,0), C = (1/2, sqrt3/2), computed exactly in Q(sqrt 3):
  1. L_M = 2: all three sides have squared length 1, a spanning tree uses two of them;
  2. L_S = sqrt 3: a Steiner minimal tree on three terminals has at most one Steiner point (Gilbert-Pollak, read), so
     L_S = min(L_M, min_S g(S)) with g(S) = |SA| + |SB| + |SC|; g is convex, and at G = (1/2, sqrt3/6) the three unit
     vectors (G - P)/|G - P| sum to (0, 0) exactly (each |G - P|^2 = 1/3), so G is the global minimum and g(G) = sqrt 3
     < 2; a floating-point scan of g over a grid finds nothing below sqrt 3 (anchor);
  3. the ratio L_S/L_M = sqrt(3)/2 = 0.8660254...; the printed 0.866 stands next to the exact expression.
Red controls (--controls): the right isosceles triangle (0,0), (1,0), (0,1) gives a ratio above sqrt(3)/2 (its Fermat
point found numerically); at that triangle's centroid the unit vectors do not sum to 0; a claimed ratio 0.8660 (below
sqrt(3)/2) is refused.
usage: python3 -I steiner_check.py [--controls]
"""
import math
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


# numbers a + b sqrt3 with rational a, b
def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def mul(x, y):
    return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def num(a=0, b=0):
    return (F(a), F(b))


def sq(p, q):
    dx, dy = sub(p[0], q[0]), sub(p[1], q[1])
    return add(mul(dx, dx), mul(dy, dy))


def main():
    ok = True
    A = (num(0), num(0))
    B = (num(1), num(0))
    C = (num(F(1, 2)), num(0, F(1, 2)))
    sides = [sq(A, B), sq(B, C), sq(C, A)]
    LM = 2
    G = (num(F(1, 2)), num(0, F(1, 6)))
    d2 = [sq(G, P) for P in (A, B, C)]
    same = all(d == num(F(1, 3)) for d in d2)
    # unit vectors (G - P)/|G - P| with |G - P| = 1/sqrt3 = sqrt3/3: multiply by sqrt3 (i.e. by (0, 1))
    units = [(mul(sub(G[0], P[0]), num(0, 1)), mul(sub(G[1], P[1]), num(0, 1))) for P in (A, B, C)]
    usum = (add(add(units[0][0], units[1][0]), units[2][0]), add(add(units[0][1], units[1][1]), units[2][1]))
    zero = usum == (num(0), num(0))
    unit_len = all(add(mul(u[0], u[0]), mul(u[1], u[1])) == num(1) for u in units)
    LS = 3 * math.sqrt(F(1, 3))                          # 3 / sqrt 3 = sqrt 3
    pts = [(0.0, 0.0), (1.0, 0.0), (0.5, math.sqrt(3) / 2)]
    g = lambda x, y: sum(math.hypot(x - px, y - py) for px, py in pts)
    scan = min(g(i / 400, j / 400) for i in range(-100, 501) for j in range(-100, 501))
    good = (all(s == num(1) for s in sides) and same and zero and unit_len and scan >= math.sqrt(3) - 1e-12
            and math.sqrt(3) < LM)
    print(f"equilateral triangle, sides squared {[str(s[0]) for s in sides]} -> L_M = {LM}")
    print(f"  G = (1/2, sqrt3/6): |G - P|^2 = 1/3 for all three: {same}; unit vectors of length 1: {unit_len}, sum = "
          f"(0, 0): {zero}; g(G) = 3/sqrt 3 = sqrt 3 = {LS:.12f}; grid scan min g = {scan:.12f}")
    print(f"  L_S = sqrt 3 < L_M = 2; ratio L_S/L_M = sqrt(3)/2 = {math.sqrt(3) / 2:.12f} (printed: sqrt(3)/2, 0.866)")
    print(f"RESULT [GP1968]: the equilateral triangle has L_S/L_M = sqrt(3)/2, so C_43 <= sqrt(3)/2 ~ 0.866: "
          f"{'HOLDS' if good else 'FAILS'}")
    ok &= good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        T = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        h = lambda x, y: sum(math.hypot(x - px, y - py) for px, py in T)
        best = min(h(i / 500, j / 500) for i in range(0, 501) for j in range(0, 501))
        r2 = best / 2.0                                  # MST of the right isosceles triangle: 1 + 1
        print(f"  right isosceles triangle: min g = {best:.6f}, ratio {r2:.6f} > {math.sqrt(3) / 2:.6f}")
        res.append(("right isosceles triangle", r2 > math.sqrt(3) / 2))
        Gc = (1 / 3, 1 / 3)
        us = [((Gc[0] - px) / math.hypot(Gc[0] - px, Gc[1] - py), (Gc[1] - py) / math.hypot(Gc[0] - px, Gc[1] - py)) for px, py in T]
        s = (sum(u[0] for u in us), sum(u[1] for u in us))
        print(f"  its centroid: unit vectors sum to ({s[0]:.4f}, {s[1]:.4f})")
        res.append(("centroid of a non-equilateral triangle", abs(s[0]) + abs(s[1]) > 1e-3))
        res.append(("a claimed 0.8660", F(8660, 10000) ** 2 * 4 < 3))   # 0.8660 < sqrt(3)/2, not implied as a bound
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
