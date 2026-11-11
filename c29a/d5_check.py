#!/usr/bin/env python3
"""Exact replay of the 29a lower-bound row [KZ1873] (teorth/optimizationproblems, main 87b4ebb): the kissing number in
dimension 5 is at least 40, "achieved by the D_5 root system (40 points), giving a kissing configuration of size 40".

A set of nonzero vectors v_1, ..., v_N of equal norm is a kissing configuration (unit spheres centred at 2 v_i/|v_i|
all touch the unit sphere at 0 and have disjoint interiors) iff every two distinct ones make an angle of at least
60 degrees, i.e. <v_i, v_j> <= |v|^2 / 2.  D_5's roots are the 40 integer vectors +-e_i +- e_j (i < j) of R^5, all
of squared norm 2, so the condition is <v_i, v_j> <= 1 for i != j.  Integers only; no download.

Checked: the 40 roots are distinct, all of squared norm 2, and all 780 pairs have inner product at most 1 (the
distribution of inner products is printed); a second route recounts the roots as all vectors of Z^5 with squared norm
2 (by enumerating {-1,0,1}^5 and {-2..2}^5 with squared norm 2). The page's other row, 10 from the cross polytope
+-e_i, is checked the same way as an anchor.
Red controls (--controls): adding (1/2, 1/2, 1/2, 1/2, 1) (squared norm 2) makes a pair with inner product 3/2 > 1
(refused); replacing e_1 + e_2 by 2 e_1 (norm 4, not 2) is refused; a claimed 41 from D_5 alone is refused (there are
exactly 40 vectors of squared norm 2 in Z^5).
usage: python3 -I d5_check.py [--controls]
"""
import itertools
import sys
from fractions import Fraction

CONTROLS = "--controls" in sys.argv


def d5_roots():
    out = []
    for i, j in itertools.combinations(range(5), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * 5
                v[i], v[j] = si, sj
                out.append(tuple(v))
    return out


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def kissing(vs):
    """(ok, reason, histogram of inner products)."""
    n2 = {dot(v, v) for v in vs}
    if len(set(vs)) != len(vs):
        return False, "repeated vector", {}
    if len(n2) != 1:
        return False, f"unequal squared norms {sorted(n2)}", {}
    r2 = n2.pop()
    hist = {}
    worst = None
    for u, v in itertools.combinations(vs, 2):
        p = dot(u, v)
        hist[p] = hist.get(p, 0) + 1
        if 2 * p > r2 and worst is None:
            worst = (u, v, p)
    if worst:
        fmt = lambda v: "(" + ", ".join(str(x) for x in v) + ")"
        return False, f"pair {fmt(worst[0])}, {fmt(worst[1])} with inner product {worst[2]} > {Fraction(r2, 2)}", hist
    return True, f"squared norm {r2}, max inner product {max(hist)} <= {Fraction(r2, 2)}", hist


def main():
    ok = True
    R = d5_roots()
    good, why, hist = kissing(R)
    second = sorted(v for v in itertools.product(range(-2, 3), repeat=5) if dot(v, v) == 2)
    route2 = second == sorted(R)
    print(f"D_5 roots +-e_i +- e_j: {len(R)} vectors; {why}")
    print(f"  inner products over the {len(R) * (len(R) - 1) // 2} pairs: "
          + ", ".join(f"{k}: {hist[k]}" for k in sorted(hist)))
    print(f"  second route: the vectors of Z^5 with squared norm 2 (all of {{-2..2}}^5 scanned): {len(second)}, "
          f"the same set: {route2}")
    good = good and route2 and len(R) == 40
    print(f"RESULT [KZ1873]: 40 points of R^5 pairwise at angle >= 60 degrees (D_5), so the kissing number in "
          f"dimension 5 is at least 40: {'HOLDS' if good else 'FAILS'}")
    ok &= good
    cross = [tuple((s if k == i else 0) for k in range(5)) for i in range(5) for s in (1, -1)]
    g2, why2, _ = kissing(cross)
    print(f"anchor (the page's trivial row): cross polytope +-e_i, {len(cross)} vectors; {why2}")
    ok &= g2
    if CONTROLS:
        print("red controls (each must be refused):")
        res = []
        half = [tuple(Fraction(x) for x in v) for v in R] + [(Fraction(1, 2),) * 4 + (Fraction(1),)]
        g, w, _ = kissing(half)
        print(f"  D_5 plus (1/2, 1/2, 1/2, 1/2, 1): {w}")
        res.append(("an added 41st vector", not g))
        R2 = [v for v in R if v != (1, 1, 0, 0, 0)] + [(2, 0, 0, 0, 0)]
        g, w, _ = kissing(R2)
        print(f"  2 e_1 for e_1 + e_2: {w}")
        res.append(("2 e_1 for e_1 + e_2", not g))
        res.append(("a claimed 41 from D_5", len(second) != 41))
        for label, r in res:
            print(f"  {label}: {'refused, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
