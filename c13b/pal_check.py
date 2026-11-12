#!/usr/bin/env python3
"""Exact replay of the areas behind the two 13b upper-bound rows [Pal1920] (teorth/optimizationproblems, main
87b4ebb): sqrt(3)/2 = 0.8660254038... (the regular hexagon circumscribing the disk of diameter 1) and
2 - 2/sqrt(3) = 0.8452994616... (that hexagon truncated "via an inscribed regular dodecagon").

Exactly in Q(sqrt 3) (pairs a + b sqrt 3 with rational a, b):
  1. the regular hexagon H with vertices at distance R = 1/sqrt 3 from 0 in the directions k * 60 degrees: every side
     line is at distance exactly 1/2 from 0 (so H circumscribes the disk of diameter 1 and has width 1 in the six
     side directions), and its area by the shoelace formula is sqrt(3)/2;
  2. the regular dodecagon circumscribed about the same disk has its sides on the six side lines of H and on the six
     lines x . u_k = 1/2, u_k the vertex directions of H; cutting two corners of H along two of the latter lines
     (Pal's truncation) leaves a polygon whose exact vertices we compute (line intersections) and whose area by the
     shoelace formula is 2 - 2/sqrt 3; each removed corner is the triangle of area (7 sqrt 3 - 12)/12, and the area does
     not depend on which two corners are cut.
That these sets are universal covers (contain a congruent copy of every set of diameter 1) is Pal's theorem, read, not
checked. The printed decimals 0.8660254038 and 0.8452994616 are compared with the exact values (enclosed through
sqrt 3 in [1.7320508075688772935274463415, ...936]).
Red controls (--controls): the hexagon with R = 1/2 (a smaller one) has side lines at distance sqrt(3)/4 < 1/2 (refused:
it does not contain the disk); cutting at distance 0.45 instead of 1/2 leaves the disk point (1/2, 0) outside (refused);
a claimed area 0.8452 for the truncated hexagon is refused (below 2 - 2/sqrt 3).
usage: python3 -I pal_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv
getcontext().prec = 40


def Q(a, b=0):
    return (F(a), F(b))


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def mul(x, y):
    return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def inv(x):
    n = x[0] * x[0] - 3 * x[1] * x[1]
    return (x[0] / n, -x[1] / n)


def div(x, y):
    return mul(x, inv(y))


def sign(x):
    a, b = x
    if a >= 0 and b >= 0:
        return (a > 0 or b > 0) and 1 or 0
    if a <= 0 and b <= 0:
        return -1
    if a > 0:
        return 1 if a * a > 3 * b * b else -1
    return 1 if 3 * b * b > a * a else -1


def to_dec(x):
    s3 = Decimal(3).sqrt()
    return Decimal(x[0].numerator) / x[0].denominator + Decimal(x[1].numerator) / x[1].denominator * s3


# unit vectors at k * 60 degrees and at 30 + k * 60 degrees
COS60 = [Q(1), Q(F(1, 2)), Q(F(-1, 2)), Q(-1), Q(F(-1, 2)), Q(F(1, 2))]
SIN60 = [Q(0), Q(0, F(1, 2)), Q(0, F(1, 2)), Q(0), Q(0, F(-1, 2)), Q(0, F(-1, 2))]


def hexagon(R):
    return [(mul(R, COS60[k]), mul(R, SIN60[k])) for k in range(6)]


def area(P):
    s = Q(0)
    for i in range(len(P)):
        x1, y1 = P[i]
        x2, y2 = P[(i + 1) % len(P)]
        s = add(s, sub(mul(x1, y2), mul(x2, y1)))
    return mul(s, Q(F(1, 2)))


def line_dist2(p, q):
    """squared distance from 0 to the line through p and q."""
    cr = sub(mul(p[0], q[1]), mul(p[1], q[0]))
    d2 = add(mul(sub(q[0], p[0]), sub(q[0], p[0])), mul(sub(q[1], p[1]), sub(q[1], p[1])))
    return div(mul(cr, cr), d2)


def cut_corner(P, k, h):
    """Replace vertex k of polygon P by the two points where the line x . u_k = h meets its two sides."""
    u = (COS60[k], SIN60[k])
    v = P[k]
    out = []
    for j in range(len(P)):
        if j != k:
            out.append(P[j])
            continue
        prev, nxt = P[j - 1], P[(j + 1) % len(P)]
        pts = []
        for a in (prev, nxt):
            # point v + t (a - v) with (v + t (a - v)) . u = h
            dv = add(mul(v[0], u[0]), mul(v[1], u[1]))
            da = add(mul(a[0], u[0]), mul(a[1], u[1]))
            t = div(sub(Q(h), dv), sub(da, dv))
            pts.append((add(v[0], mul(t, sub(a[0], v[0]))), add(v[1], mul(t, sub(a[1], v[1])))))
        out.extend([pts[0], pts[1]])
    return out


def main():
    ok = True
    R = Q(0, F(1, 3))                     # 1/sqrt 3 = sqrt 3 / 3
    H = hexagon(R)
    d2 = [line_dist2(H[i], H[(i + 1) % 6]) for i in range(6)]
    A = area(H)
    good1 = all(d == Q(F(1, 4)) for d in d2) and A == Q(0, F(1, 2))
    print("13b row [Pal1920] hexagon: all six side lines at squared distance 1/4 from 0: %s; area = %s + %s sqrt3"
          % (all(d == Q(F(1, 4)) for d in d2), A[0], A[1]))
    print("  sqrt(3)/2 = %s...; the printed 0.8660254038 is its rounding to 10 decimals (up, next to the exact value): %s"
          % (str(to_dec(A))[:14], str(to_dec(A).quantize(Decimal("1e-10"))) == "0.8660254038"))
    print("RESULT [Pal1920] sqrt(3)/2: the regular hexagon of width 1 has area exactly sqrt(3)/2: %s"
          % ("HOLDS" if good1 else "FAILS"))
    ok &= good1
    T = cut_by_vertex(cut_corner(H, 0, F(1, 2)), H[2], 2, F(1, 2))       # corners 0 and 2 cut
    AT = area(T)
    tri = sub(A, AT)
    good2 = AT == Q(2, F(-2, 3)) and tri == mul(Q(2), Q(-1, F(7, 12)))
    print("13b row [Pal1920] truncated: %d vertices; area = %s + %s sqrt3 = 2 - 2/sqrt3: %s; removed 2 x (7 sqrt3 - 12)/12: %s"
          % (len(T), AT[0], AT[1], AT == Q(2, F(-2, 3)), tri == mul(Q(2), Q(-1, F(7, 12)))))
    # the other choices of two corners give the same area
    same = all(area(cut_by_vertex(cut_corner(H, 0, F(1, 2)), H[j], j, F(1, 2))) == AT for j in (1, 2, 3, 4, 5))
    print("  cutting corner 0 and any other corner j = 1..5 gives the same area: %s" % same)
    print("  2 - 2/sqrt3 = %s...; the printed 0.8452994616 is its rounding (and truncation) to 10 decimals: %s"
          % (str(to_dec(AT))[:14], str(to_dec(AT).quantize(Decimal("1e-10"))) == "0.8452994616"))
    good2 &= same
    print("RESULT [Pal1920] 2 - 2/sqrt3: the hexagon minus two dodecagon corners has area exactly 2 - 2/sqrt3: %s"
          % ("HOLDS" if good2 else "FAILS"))
    ok &= good2
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        Hs = hexagon(Q(F(1, 2)))
        d = line_dist2(Hs[0], Hs[1])
        print("  hexagon with R = 1/2: side lines at squared distance %s < 1/4: %s" % (d[0], "refused, as expected"
                                                                                      if sign(sub(d, Q(F(1, 4)))) < 0 else "UNEXPECTED"))
        c_ok &= sign(sub(d, Q(F(1, 4)))) < 0
        Tbad = cut_corner(H, 0, F(45, 100))
        inside = contains(Tbad, (Q(F(1, 2)), Q(0)))
        print("  cut at distance 0.45: the disk point (1/2, 0) %s the polygon: %s"
              % ("lies in" if inside else "is outside", "UNEXPECTED" if inside else "refused, as expected"))
        c_ok &= not inside
        claim = sign(sub(AT, Q(F(8452, 10000)))) <= 0
        print("  a claimed area 0.8452 for the truncated hexagon: %s" % ("holds (UNEXPECTED)" if claim else
                                                                        "refused, as expected (the area is larger)"))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


def contains(P, pt):
    """pt in the closed convex polygon P (counterclockwise), exactly."""
    for i in range(len(P)):
        a, b = P[i], P[(i + 1) % len(P)]
        cr = sub(mul(sub(b[0], a[0]), sub(pt[1], a[1])), mul(sub(b[1], a[1]), sub(pt[0], a[0])))
        if sign(cr) < 0:
            return False
    return True


def cut_by_vertex(P, v, k, h):
    """Cut the corner of P at the original hexagon vertex v (direction index k) by the line x . u_k = h."""
    j = P.index(v)
    u = (COS60[k], SIN60[k])
    prev, nxt = P[j - 1], P[(j + 1) % len(P)]
    pts = []
    for a in (prev, nxt):
        dv = add(mul(v[0], u[0]), mul(v[1], u[1]))
        da = add(mul(a[0], u[0]), mul(a[1], u[1]))
        t = div(sub(Q(h), dv), sub(da, dv))
        pts.append((add(v[0], mul(t, sub(a[0], v[0]))), add(v[1], mul(t, sub(a[1], v[1])))))
    return P[:j] + pts + P[j + 1:]


if __name__ == "__main__":
    sys.exit(main())
