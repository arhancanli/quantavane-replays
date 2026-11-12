#!/usr/bin/env python3
"""Exact replay of the 27a upper-bound row [Had1945], [CR2017] (teorth/optimizationproblems, main 87b4ebb):
chi(R^2) <= 7, "a periodic 7-coloring obtained from a tiling of the plane by small regular hexagons".

The construction as the row describes it and as Cranston and Rabern recall it in words (arXiv 1501.01647v1, e-print
sha256 fad197e5...8b61, introduction, read as text: hexagons, seven colours, cells of one colour more than 1 apart);
Hadwiger's paper was not read. Regular hexagons of side s (circumradius s), the Voronoi cells of the lattice spanned by
a = (sqrt3 s, 0) and b = (sqrt3 s/2, 3s/2); the cell with centre u a + v b gets colour (u + 3v) mod 7; boundary
points go to one of their cells (any fixed rule). Since |u a + v b|^2 = 3 s^2 N(u, v) with N = u^2 + uv + v^2, what is
checked, exactly (rational arithmetic only):
  1. the colour classes are the cosets of L7 = {(u, v): u + 3v = 0 mod 7}, a sublattice of index 7, and
     min N over L7 minus 0 is 7 (enumeration over |u|, |v| <= 12; outside, N >= (3/4) max(|u|,|v|)^2 > 7);
  2. two points of one closed cell are at distance <= 2s; two points of different cells of one colour are at distance
     >= sqrt(21) s - 2s (centre distance minus two circumradii);
  3. at s = 2/5: 2s = 4/5 < 1, and sqrt(21) s - 2s > 1 because 21 s^2 > (1 + 2s)^2 (84/25 > 81/25); so no two points
     at distance exactly 1 share a colour, and chi(R^2) <= 7. The whole admissible range of this crude argument,
     1/(sqrt 21 - 2) = 0.38720... < s < 1/2, is printed, and the true least distance between two cells of one colour
     (polygon distance, floating point, as an anchor) at s = 2/5.
Red controls (--controls): s = 1/2 (a cell of diameter 1) is refused; s = 0.35 is refused, and there the true least
distance between two cells of one colour is below 1 (floating point), so the colouring really fails; the colouring
(u + 2v) mod 7 is refused (L7 then contains (-2, 1) with N = 3).
usage: python3 -I hexagon_check.py [--controls]      (Python 3.9+, standard library only)
"""
import math
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


def min_norm(colour_mult, R=12):
    """min N(u, v) over nonzero (u, v), |u|,|v| <= R, with u + colour_mult * v = 0 mod 7; and the index of L7."""
    best, arg = None, None
    for u in range(-R, R + 1):
        for v in range(-R, R + 1):
            if (u, v) != (0, 0) and (u + colour_mult * v) % 7 == 0:
                n = u * u + u * v + v * v
                if best is None or n < best:
                    best, arg = n, (u, v)
    # outside the box: N >= (3/4) max(|u|, |v|)^2 >= (3/4) (R+1)^2
    tail = F(3, 4) * (R + 1) ** 2
    classes = {(u + colour_mult * v) % 7 for u in range(7) for v in range(7)}
    return best, arg, tail, len(classes)


def hexagon(cx, cy, s):
    return [(cx + s * math.cos(math.pi / 6 + k * math.pi / 3), cy + s * math.sin(math.pi / 6 + k * math.pi / 3))
            for k in range(6)]


def seg_dist(p, a, b):
    ax, ay = a; bx, by = b; px, py = p
    dx, dy = bx - ax, by - ay
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)))
    return math.hypot(px - ax - t * dx, py - ay - t * dy)


def poly_dist(P, Q):
    d = float("inf")
    for X, Y in ((P, Q), (Q, P)):
        for p in X:
            for i in range(6):
                d = min(d, seg_dist(p, Y[i], Y[(i + 1) % 6]))
    return d


def true_min_same_colour(s, colour_mult=3, R=4):
    """Least distance (float) between the cell at 0 and the cells of the same colour with |u|,|v| <= R."""
    r3 = math.sqrt(3)
    H0 = hexagon(0.0, 0.0, s)
    best = float("inf")
    for u in range(-R, R + 1):
        for v in range(-R, R + 1):
            if (u, v) != (0, 0) and (u + colour_mult * v) % 7 == 0:
                cx, cy = u * r3 * s + v * r3 * s / 2, v * 1.5 * s
                best = min(best, poly_dist(H0, hexagon(cx, cy, s)))
    return best


def certify(s, colour_mult=3):
    """The exact argument at side s (a Fraction): returns (ok, lines)."""
    lines = []
    m, arg, tail, idx = min_norm(colour_mult)
    lines.append("colour (u + %dv) mod 7: %d classes; least N on L7 minus 0 is %d at %s (outside the box N >= %s)"
                 % (colour_mult, idx, m, arg, tail))
    ok = idx == 7 and tail > m
    diam_ok = 2 * s < 1
    lines.append("cell diameter 2s = %s < 1: %s" % (2 * s, diam_ok))
    # same-colour separation: sqrt(3 m) s - 2 s > 1  <=>  3 m s^2 > (1 + 2s)^2
    sep_ok = 3 * m * s * s > (1 + 2 * s) ** 2
    lines.append("same-colour cells: centre distance >= sqrt(%d) s, gap >= (sqrt(%d) - 2) s > 1 iff %d s^2 > (1 + 2s)^2: "
                 "%s > %s: %s" % (3 * m, 3 * m, 3 * m, 3 * m * s * s, (1 + 2 * s) ** 2, sep_ok))
    return ok and diam_ok and sep_ok, lines


def main():
    ok = True
    print("27a row [Had1945], [CR2017]: chi(R^2) <= 7 from the hexagon tiling, 7 colours")
    s = F(2, 5)
    good, lines = certify(s)
    for ln in lines:
        print("  s = 2/5: " + ln)
    ok &= good
    lo = 1 / (math.sqrt(21) - 2)
    print("  admissible range of this argument: %.10f < s < 0.5 (2/5 inside: %s)" % (lo, lo < 0.4 < 0.5))
    print("  anchor (floating point): true least distance between two cells of one colour at s = 2/5: %.6f"
          % true_min_same_colour(0.4))
    print("RESULT [Had1945/CR2017]: at s = 2/5 no two points at distance 1 get the same colour, so chi(R^2) <= 7: %s"
          % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g, _ = certify(F(1, 2))
        print("  s = 1/2 (cells of diameter 1): %s" % ("passes (UNEXPECTED)" if g else "refused, as expected"))
        c_ok &= not g
        g, _ = certify(F(35, 100))
        tm = true_min_same_colour(0.35)
        print("  s = 0.35: %s; the true least same-colour distance there is %.6f %s 1"
              % ("passes (UNEXPECTED)" if g else "refused, as expected", tm, "<" if tm < 1 else ">="))
        c_ok &= (not g) and tm < 1
        g, lines = certify(F(2, 5), colour_mult=2)
        print("  colouring (u + 2v) mod 7 at s = 2/5: %s (%s)" % ("passes (UNEXPECTED)" if g else "refused, as expected",
                                                                lines[0]))
        c_ok &= not g
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
