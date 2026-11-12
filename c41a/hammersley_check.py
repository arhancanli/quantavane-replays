#!/usr/bin/env python3
"""Replay of the 41a lower-bound row [Hammersley1968] (teorth/optimizationproblems, main 87b4ebb): A >= pi/2 + 2/pi,
Hammersley's sofa (the shape as it is classically described; Hammersley's note itself was not read).

The sofa S_r (r = 2/pi): the half-stadium D_r = {0 <= y <= 1, |x| <= r} union the quarter disks of radius 1 centred
at (+-r, 0) (outer sides, y >= 0), minus the open half-disk {x^2 + y^2 < r^2, y > 0}. Hallway (world frame):
{0 <= y <= 1, x <= 1} union {0 <= x <= 1, y <= 1}, inner corner at the origin. The motion: rotate S_r clockwise by
t in [0, pi/2] and translate so that the world inner corner sits at the sofa point c(t) = r (cos 2t, sin 2t) on the
hole's circle. In the sofa frame the hallway is {u <= 1, v <= 1, max(u, v) >= 0} with u = <p - c, e1>, v = <p - c, e2>,
e1 = (cos t, sin t), e2 = (-sin t, cos t). What is checked:
  1. (exact, for all t: identities of rational functions of s = tan(t/2) in [0, 1]) <c, e1> = r cos t,
     <c, e2> = r sin t, and the rays from c to (r, 0) and to (-r, 0) point along -e2 and -e1 (Thales), so the open
     quadrant {u < 0, v < 0} meets {y >= 0} only inside the triangle c, (r, 0), (-r, 0), which lies in the closed
     half-disk: no point of S_r is in it;
  2. the support function of D_r in a direction (a, b) with b >= 0 is 1 + r |a| (the maximum over each quarter arc is
     +-r a + cos(phi - theta) <= +-r a + 1, the rectangle corners give less), so max_S u = 1 + r cos t - <c, e1> = 1
     and max_S v = 1: the outer walls touch and are never crossed; this is confirmed numerically on 2001 values of t
     and about 9000 boundary points of S_r (world coordinates, tolerance 1e-12);
  3. (exact, pi as a symbol) area(S_r) = pi/2 + 2r - pi r^2/2, maximal at r = 2/pi, where it equals pi/2 + 2/pi; and
     0 < 2/pi <= 1 (the hole fits under the top). The end positions slide freely along each arm.
Red controls (--controls): r = 6/5 > 1 is refused (the hole would cut the top); the path c(t) = r (cos t, sin t) fails
the identity <c, e1> = r cos t and leaves sofa points outside the hallway; a sofa whose half-stadium is widened to
|x| <= r + 1/10 (hole unchanged) crosses a wall at some t; a claimed area 2.2075 (above pi/2 + 2/pi = 2.20741...) is
refused.
usage: python3 -I hammersley_check.py [--controls]      (Python 3.9+, standard library only)
"""
import math
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


# ---------------------------------------------------------------- polynomials in s (exact), rational functions as (num, den)
def padd(a, b):
    n = max(len(a), len(b))
    r = [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
    while r and r[-1] == 0:
        r.pop()
    return r


def pmul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    while r and r[-1] == 0:
        r.pop()
    return r


def rf(num, den=(1,)):
    return (list(num), list(den))


def radd(x, y):
    return (padd(pmul(x[0], y[1]), pmul(y[0], x[1])), pmul(x[1], y[1]))


def rneg(x):
    return ([-c for c in x[0]], x[1])


def rmul(x, y):
    return (pmul(x[0], y[0]), pmul(x[1], y[1]))


def is_zero(x):
    return not x[0]


ONE_PLUS_S2 = [1, 0, 1]
COS = rf([1, 0, -1], ONE_PLUS_S2)          # cos t = (1 - s^2)/(1 + s^2)
SIN = rf([0, 2], ONE_PLUS_S2)              # sin t = 2s/(1 + s^2)
COS2 = radd(rmul(COS, COS), rneg(rmul(SIN, SIN)))
SIN2 = rmul(rf([2]), rmul(SIN, COS))


def identities(cx, cy):
    """With c = r (cx, cy) (rational functions of s), check <c,e1> = r cos t, <c,e2> = r sin t and the Thales rays."""
    i1 = radd(radd(rmul(cx, COS), rmul(cy, SIN)), rneg(COS))
    i2 = radd(radd(rmul(rneg(cx), SIN), rmul(cy, COS)), rneg(SIN))
    # ray from c to (r, 0) is r (1 - cx, -cy); parallel to -e2 = (sin t, -cos t) iff cross product 0
    cr1 = radd(rmul(radd(rf([1]), rneg(cx)), rneg(COS)), rneg(rmul(rneg(cy), SIN)))
    # ray from c to (-r, 0) is r (-1 - cx, -cy); parallel to -e1 = (-cos t, -sin t)
    cr2 = radd(rmul(radd(rf([-1]), rneg(cx)), rneg(SIN)), rneg(rmul(rneg(cy), rneg(COS))))
    return is_zero(i1), is_zero(i2), is_zero(cr1), is_zero(cr2)


# ---------------------------------------------------------------- numeric containment (anchor)
def sofa_boundary(r, widen=0.0, n=1500):
    pts = []
    w = r + widen
    for k in range(n + 1):
        a = (math.pi / 2) * k / n
        pts.append((w + math.cos(a), math.sin(a)))          # right quarter arc
        pts.append((-w - math.cos(a), math.sin(a)))         # left quarter arc
        pts.append((r * math.cos(math.pi * k / n), r * math.sin(math.pi * k / n)))   # hole arc
        pts.append((-w + 2 * w * k / n, 1.0))               # top
    for k in range(n + 1):
        x = -w - 1 + (w + 1 - r) * k / n
        pts.append((x, 0.0))
        pts.append((-x, 0.0))                               # bottom pieces outside the hole
    return pts


def in_hallway(X, Y, tol=1e-12):
    return (-tol <= Y <= 1 + tol and X <= 1 + tol) or (-tol <= X <= 1 + tol and Y <= 1 + tol)


def worst(r, path, widen=0.0, nt=2000):
    P = sofa_boundary(r, widen)
    bad = 0
    for j in range(nt + 1):
        t = (math.pi / 2) * j / nt
        cx, cy = path(t)
        c, s = math.cos(t), math.sin(t)
        for (px, py) in P:
            # world = R(-t)(p - c)
            dx, dy = px - r * cx, py - r * cy
            X, Y = c * dx + s * dy, -s * dx + c * dy
            if not in_hallway(X, Y):
                bad += 1
    return bad


def main():
    ok = True
    print("41a row [Hammersley1968]: A >= pi/2 + 2/pi from Hammersley's sofa (r = 2/pi)")
    ids = identities(COS2, SIN2)
    print("  1. exact for all t (s = tan(t/2)): <c,e1> = r cos t: %s; <c,e2> = r sin t: %s; ray to (r,0) along -e2: %s; "
          "ray to (-r,0) along -e1: %s" % ids)
    ok &= all(ids)
    r = 2 / math.pi
    bad = worst(r, lambda t: (math.cos(2 * t), math.sin(2 * t)))
    print("  2. numeric anchor: 2001 angles x %d boundary points of S_r, all inside the hallway: %s (violations %d)"
          % (len(sofa_boundary(r)), bad == 0, bad))
    ok &= bad == 0
    # area(r) = pi/2 + 2r - (pi/2) r^2 with pi a symbol; numbers a + b pi + c/pi + d/pi^2 as tuples (a, b, c, d)
    def tmul(x, y):          # product, keeping only the powers pi^1, pi^0, pi^-1, pi^-2 that occur here
        out = {}
        for i, a in zip((0, 1, -1, -2), x):
            for j, b in zip((0, 1, -1, -2), y):
                if a and b:
                    out[i + j] = out.get(i + j, F(0)) + a * b
        assert set(out) <= {0, 1, -1, -2}, out
        return tuple(out.get(k, F(0)) for k in (0, 1, -1, -2))
    tadd = lambda x, y: tuple(a + b for a, b in zip(x, y))
    PI_HALF = (F(0), F(1, 2), F(0), F(0))
    R = (F(0), F(0), F(2), F(0))                       # r = 2/pi
    area4 = tadd(tadd(PI_HALF, tmul((F(2), 0, 0, 0), R)), tmul((F(0), F(-1, 2), F(0), F(0)), tmul(R, R)))
    deriv = tadd((F(2), 0, 0, 0), tmul((F(0), F(-1), F(0), F(0)), R))        # 2 - pi r
    assert deriv == (0, 0, 0, 0)
    area_pi = area4[:3]
    print("  3. area(S_r) = pi/2 + 2r - pi r^2/2; at r = 2/pi: %s + %s pi + %s/pi = pi/2 + 2/pi: %s; d/dr = 2 - pi r = 0 "
          "there (the maximum over r); 0 < 2/pi <= 1: %s" % (area_pi[0], area_pi[1], area_pi[2],
                                                          area_pi == (0, F(1, 2), F(2)), 0 < r <= 1))
    ok &= area_pi == (0, F(1, 2), F(2)) and 0 < r <= 1
    print("     pi/2 + 2/pi = %.12f" % (math.pi / 2 + 2 / math.pi))
    print("RESULT [Hammersley1968]: the sofa moves around the corner (exact path identities, numeric anchor) and has "
          "area pi/2 + 2/pi, so A >= pi/2 + 2/pi = 2.20741609916...: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        print("  r = 6/5: hole fits under the top (r <= 1): %s" % ("yes (UNEXPECTED)" if F(6, 5) <= 1 else
                                                                  "no, refused as expected"))
        c_ok &= not F(6, 5) <= 1
        ids2 = identities(COS, SIN)
        bad2 = worst(r, lambda t: (math.cos(t), math.sin(t)), nt=400)
        print("  path c(t) = r (cos t, sin t): identity <c,e1> = r cos t %s; numeric violations %d: %s"
              % ("holds" if ids2[0] else "fails", bad2, "refused, as expected" if (not ids2[0] and bad2 > 0) else
                 "UNEXPECTED"))
        c_ok &= (not ids2[0]) and bad2 > 0
        bad3 = worst(r, lambda t: (math.cos(2 * t), math.sin(2 * t)), widen=0.1, nt=400)
        print("  half-stadium widened to |x| <= r + 1/10: numeric violations %d: %s"
              % (bad3, "refused, as expected" if bad3 > 0 else "UNEXPECTED"))
        c_ok &= bad3 > 0
        claim = 2.2075 <= math.pi / 2 + 2 / math.pi
        print("  a claimed area 2.2075: %s" % ("holds (UNEXPECTED)" if claim else "refused, as expected"))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
