#!/usr/bin/env python3
"""13b, lower bound [Elek1994] 0.8257, program B (written separately from program A): branch and bound over the
triangle's position, with no convexity argument.

m(v) = area conv(D u (T + v)), D the disk of radius 1/2 at 0, T the equilateral triangle of side 1 with centroid v
(vertex q_0 = (0, 1/sqrt 3), q_1, q_2 at +-120 degrees). Two lower bounds, valid for every v:
  (i)  m(v) >= pi r^2 + g(|p|) for each vertex p = q_k + v (the hull contains conv(D u {p})), where
       g(d) = (1/2) sqrt(d^2 - 1/4) - (1/4) arccos(1/(2d)) for d > 1/2 and 0 otherwise (increasing);
  (ii) if the caps' angular sectors are pairwise disjoint, the three caps are disjoint and m(v) >= pi r^2 + sum_k g(|p_k|).
       Disjointness is checked on each cell directly: the angle between p_j and p_k (from its cosine, enclosed) must be
       at least arccos(1/(2|p_j|)) + arccos(1/(2|p_k|)) for every v in the cell.
For |v| >= 1/2 some vertex has |p| >= sqrt(1/4 + 1/3) (since sum q_k = 0, some q_k . v >= 0), and (i) gives more than
0.85. The square [-1/2, 1/2]^2 is cut into cells; on a cell, |q_k + v| is enclosed exactly by the nearest and
farthest points of the box, g is increasing, and each cell needs max(i, ii) >= 0.8257 at the lower ends; cells that
fail are split into four, to depth 22. The minimum certified lower bound over all cells is reported.
Arithmetic: Python decimal at 40 digits; square roots and arccos (as 2 arctan of a half-angle tangent, our own series)
are moved outward by 10^-30; g at the nearest point is taken at d rounded down to a multiple of 10^-9 (g increasing),
so its values can be cached. (A first version wrote the disk's area as pi/4 even for the radius control; that
control passed, which exposed the slip; the area is now pi r^2 throughout.) Red controls (--control): the claim 0.8258; a disk of radius 0.49; bound (ii) removed
(the single-cap bound (i) alone gives only about 0.7988 near v = 0, so the search must fail). Standard library only.
"""
import sys
from decimal import Decimal as D, getcontext, ROUND_FLOOR

getcontext().prec = 40
PAD = D(10) ** -30
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
R = D("0.49") if CONTROL == "radius" else D("0.5")
TARGET = D("0.8258") if CONTROL == "claim" else D("0.8257")


def atan_small(x):
    s = D(0); p = x; k = 0
    while abs(p) > D(10) ** -38:
        s += p / (2 * k + 1) if k % 2 == 0 else -(p / (2 * k + 1))
        p *= x * x; k += 1
    return s


def atan(x):
    m = 1
    while x > D("0.2"):
        x = x / (1 + (1 + x * x).sqrt()); m *= 2
    return atan_small(x) * m


PI = 16 * atan_small(D(1) / 5) - 4 * atan_small(D(1) / 239)


def arccos(x):   # 0 <= x <= 1
    if x >= 1:
        return D(0)
    return 2 * atan(((1 - x) / (1 + x)).sqrt())


GCACHE = {}


def g_lo(d):
    """a lower bound for g at any distance >= d"""
    if d <= R:
        return D(0)
    key = d.quantize(D("1e-9"), rounding=ROUND_FLOOR)
    if key <= R:
        return D(0)
    if key not in GCACHE:
        GCACHE[key] = R * (key * key - R * R).sqrt() - R * R * arccos(R / key) - PAD
    return GCACHE[key]


S3 = D(3).sqrt()
Q = [(D(0), 1 / S3), (-S3 / 2 / S3, -1 / (2 * S3)), (S3 / 2 / S3, -1 / (2 * S3))]   # (0, 1/sqrt3), (-1/2, -1/(2 sqrt3)), (1/2, ...)


def box_dist(cx, cy, x0, x1, y0, y1):
    """min and max distance from the point (cx, cy) to the box"""
    dx = max(x0 - cx, D(0), cx - x1); dy = max(y0 - cy, D(0), cy - y1)
    fx = max(abs(cx - x0), abs(cx - x1)); fy = max(abs(cy - y0), abs(cy - y1))
    return (dx * dx + dy * dy).sqrt() - PAD, (fx * fx + fy * fy).sqrt() + PAD


def cell_bound(x0, x1, y0, y1):
    dmin = []; dmax = []
    for qx, qy in Q:
        a, b = box_dist(-qx, -qy, x0, x1, y0, y1)   # |q + v| = distance from v to -q
        dmin.append(max(a, D(0))); dmax.append(b)
    b1 = PI * R * R + max(g_lo(d) for d in dmin) - PAD
    if CONTROL == "no-ii":
        return b1
    # sector test for each pair: the angle between p_j, p_k is at least arccos(cos_max) with
    # cos = (|p_j|^2 + |p_k|^2 - 1)/(2 |p_j| |p_k|) <= (dmax_j^2 + dmax_k^2 - 1)/(2 dmin_j dmin_k) if the numerator >= 0,
    # or <= (dmin_j^2 + dmin_k^2 - 1)/(2 dmax_j dmax_k) if the numerator is negative even at the upper ends
    disjoint = True
    for j in range(3):
        for k in range(j + 1, 3):
            num_hi = dmax[j] ** 2 + dmax[k] ** 2 - 1
            if num_hi >= 0:
                if dmin[j] * dmin[k] == 0:
                    disjoint = False; break
                cmax = num_hi / (2 * dmin[j] * dmin[k]) + PAD
            else:
                cmax = num_hi / (2 * dmax[j] * dmax[k]) + PAD
            if cmax >= 1:
                disjoint = False; break
            gam = arccos(max(cmax, D(-1))) if cmax > -1 else PI
            alj = arccos(min(R / dmax[j], D(1))) + PAD
            alk = arccos(min(R / dmax[k], D(1))) + PAD
            if gam - PAD < alj + alk:
                disjoint = False; break
        if not disjoint:
            break
    if not disjoint:
        return b1
    b2 = PI * R * R + sum(g_lo(d) for d in dmin) - PAD
    return max(b1, b2)


def main():
    assert D("3.141592653589793238462643383279502884") < PI + PAD and PI - PAD < D("3.141592653589793238462643383279502885")
    # far field: |v| >= 1/2 gives a vertex with |p|^2 >= 1/4 + 1/3
    far = PI * R * R + g_lo((D(1) / 4 + D(1) / 3).sqrt() - PAD)
    print("radius %s; far field |v| >= 1/2: m >= %s (>= target %s: %s)" % (R, str(far)[:12], TARGET, far >= TARGET))
    ok = far >= TARGET
    stack = [(D("-0.5"), D("0.5"), D("-0.5"), D("0.5"), 0)]
    cells = 0; worst = None; failed = None; maxdepth = 0
    while stack:
        x0, x1, y0, y1, dep = stack.pop()
        b = cell_bound(x0, x1, y0, y1)
        if b >= TARGET:
            cells += 1
            worst = b if worst is None or b < worst else worst
            maxdepth = max(maxdepth, dep)
            continue
        if dep >= 22:
            failed = (x0, x1, y0, y1, b); break
        xm = (x0 + x1) / 2; ym = (y0 + y1) / 2
        stack += [(x0, xm, y0, ym, dep + 1), (xm, x1, y0, ym, dep + 1), (x0, xm, ym, y1, dep + 1), (xm, x1, ym, y1, dep + 1)]
    if failed:
        print("NOT CERTIFIED: a cell at depth 22 near (%s, %s) has lower bound %s < %s"
              % (str(failed[0])[:10], str(failed[2])[:10], str(failed[4])[:12], TARGET))
        ok = False
    else:
        print("every cell of [-1/2, 1/2]^2 certified: %d cells, depth up to %d, smallest certified lower bound %s"
              % (cells, maxdepth, str(worst)[:14]))
    print("so area conv(D u (T + v)) >= %s for every v: %s" % (TARGET, ok))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
