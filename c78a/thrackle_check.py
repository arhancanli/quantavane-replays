#!/usr/bin/env python3
"""Exact replay of the 78a lower-bound row [FP2011] (teorth/optimizationproblems, main 87b4ebb):
C_78 >= 1, "any cycle of length at least five can be drawn as a thrackle, and such a cycle has |E|/|V| = 1".

C_78 = sup |E(G)|/|V(G)| over graphs with a thrackle drawing (every pair of edges meets exactly once: at a common
endpoint or at a proper crossing), so one thrackled cycle gives C_78 >= 1. We check the classical straight-line
drawings of the odd cycles: for odd n = 5..101, points P_i = (i, i^2), i = 0..n-1 (strictly convex position, hull order
0, 1, ..., n-1), and the cycle P_0, P_k, P_2k, ... with k = (n-1)/2 (the star polygon {n/k}). With integer orientation
tests we check, for every pair of edges:
  1. it is a Hamiltonian cycle (n distinct vertices, n edges);
  2. adjacent edges share exactly their common endpoint (no three points collinear on a parabola);
  3. every pair of non-adjacent edges crosses properly (strict opposite orientations both ways), so meets exactly once.
Even cycles are not checked (they have thrackle drawings, but no straight-line ones); they are not needed for C_78 >= 1.
Red controls (--controls): the convex n-gon (k = 1) is refused (non-adjacent edges do not cross); for n = 6 every
star order k = 1, 2 is refused; a pentagram with two vertices swapped is refused.
usage: python3 -I thrackle_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys

CONTROLS = "--controls" in sys.argv


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def cross_properly(p, q, r, s):
    o1, o2, o3, o4 = orient(p, q, r), orient(p, q, s), orient(r, s, p), orient(r, s, q)
    return o1 * o2 < 0 and o3 * o4 < 0


def meet_only_at_shared(p, q, r):
    """Segments pq and pr (shared endpoint p): they meet only at p iff p, q, r are not collinear."""
    return orient(p, q, r) != 0


def is_thrackle(order, pts):
    n = len(order)
    if len(set(order)) != n:
        return False, "not a cycle on distinct vertices"
    E = [(pts[order[i]], pts[order[(i + 1) % n]], order[i], order[(i + 1) % n]) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b, ia, ib = E[i]
            c, d, ic, id_ = E[j]
            shared = {ia, ib} & {ic, id_}
            if shared:
                (s,) = shared
                p = pts[s]
                q = b if ia == s else a
                r = d if ic == s else c
                if not meet_only_at_shared(p, q, r):
                    return False, "adjacent edges %d, %d overlap" % (i, j)
            elif not cross_properly(a, b, c, d):
                return False, "edges %d and %d do not cross" % (i, j)
    return True, "every pair of edges meets exactly once"


def star(n, k):
    return [(i * k) % n for i in range(n)]


def main():
    ok = True
    print("78a row [FP2011]: C_78 >= 1 from thrackled cycles (straight-line star polygons of the odd cycles)")
    for n in range(5, 102, 2):
        pts = [(i, i * i) for i in range(n)]
        good, why = is_thrackle(star(n, (n - 1) // 2), pts)
        if n <= 9 or n % 24 == 5 or n == 101 or not good:
            print("  n = %3d: star {%d/%d} on (i, i^2): %s (%s); |E|/|V| = %d/%d" % (n, n, (n - 1) // 2, good, why, n, n))
        ok &= good
    print("  (all odd n from 5 to 101 checked; lines shown for a sample)")
    print("RESULT [FP2011]: each odd cycle C_n, 5 <= n <= 101, has a straight-line thrackle drawing with |E|/|V| = 1, "
          "so C_78 >= 1: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        pts = [(i, i * i) for i in range(7)]
        g, why = is_thrackle(star(7, 1), pts)
        print("  the convex 7-gon (k = 1): %s (%s)" % ("thrackle (UNEXPECTED)" if g else "refused, as expected", why))
        c_ok &= not g
        pts6 = [(i, i * i) for i in range(6)]
        for k in (1, 2):
            order = [(i * k) % 6 for i in range(6)] if k == 1 else [0, 2, 4, 1, 3, 5]
            g, why = is_thrackle(order, pts6)
            print("  n = 6, order %s: %s (%s)" % (order, "thrackle (UNEXPECTED)" if g else "refused, as expected", why))
            c_ok &= not g
        pts5 = [(i, i * i) for i in range(5)]
        order = star(5, 2)
        order[1], order[2] = order[2], order[1]
        g, why = is_thrackle(order, pts5)
        print("  pentagram with two vertices swapped %s: %s (%s)" % (order, "thrackle (UNEXPECTED)" if g else
                                                                     "refused, as expected", why))
        c_ok &= not g
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
