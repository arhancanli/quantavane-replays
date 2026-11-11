#!/usr/bin/env python3
"""Exact replay of the 27a lower-bound row [MM1961] (teorth/optimizationproblems, main 87b4ebb): chi(R^2) >= 4 from
the Moser spindle, "a 7-vertex unit-distance graph with chromatic number 4".

The spindle, built here from its standard description: two rhombi O B C D made of two unit equilateral triangles
(O B D and O C D) with the common vertex O, the second rotated about O by the angle theta with cos theta = 5/6, so
that the far tips D and D' are at distance 1 (|D| = sqrt 3 and 2 * 3 (1 - cos theta) = 1).  Coordinates lie in
Q(sqrt 3, sqrt 11); squared distances are computed exactly in that field (rational coordinates on the basis
1, sqrt 3, sqrt 11, sqrt 33).

Checked: the 7 points are distinct; the 11 spindle edges have squared length exactly 1; of all 21 pairs exactly these
11 are at distance 1 (so the unit-distance graph on the 7 points is the spindle); no proper 3-colouring exists (all
assignments with O coloured 0, which loses nothing by symmetry) and a proper 4-colouring exists.  A colouring of
the plane without monochromatic unit pairs restricts to a proper colouring of this graph, so chi(R^2) >= 4.  No
download.
Red controls (--controls): cos theta = 4/5 puts the tips at squared distance 6/5 (the spindle edge fails); the
spindle without the edge D D' is 3-colourable; the 3-colouring found for that graph violates the edge D D'.
usage: python3 -I spindle_check.py [--controls]
"""
import itertools
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv
# a + b sqrt3 + c sqrt11 + d sqrt33 as a tuple of four Fractions
MUL = {(0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (1, 2), (0, 3): (1, 3),
       (1, 1): (3, 0), (1, 2): (1, 3), (1, 3): (3, 2),
       (2, 2): (11, 0), (2, 3): (11, 1), (3, 3): (33, 0)}


def mul(x, y):
    out = [F(0)] * 4
    for i in range(4):
        if x[i] == 0:
            continue
        for j in range(4):
            if y[j] == 0:
                continue
            k, l = (i, j) if i <= j else (j, i)
            f, b = MUL[(k, l)]
            out[b] += f * x[i] * y[j]
    return tuple(out)


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def num(a=0, b=0, c=0, d=0):
    return (F(a), F(b), F(c), F(d))


def rot(p, cos, sin):
    x, y = p
    return (sub(mul(x, cos), mul(y, sin)), add(mul(x, sin), mul(y, cos)))


def sqdist(p, q):
    dx, dy = sub(p[0], q[0]), sub(p[1], q[1])
    return add(mul(dx, dx), mul(dy, dy))


def spindle(cos, sin):
    O = (num(), num())
    B = (num(b=F(1, 2)), num(F(1, 2)))
    C = (num(b=F(1, 2)), num(F(-1, 2)))
    D = (num(b=1), num())
    pts = {"O": O, "B": B, "C": C, "D": D}
    for k in "BCD":
        pts[k + "'"] = rot(pts[k], cos, sin)
    edges = [("O", "B"), ("O", "C"), ("B", "C"), ("B", "D"), ("C", "D"),
             ("O", "B'"), ("O", "C'"), ("B'", "C'"), ("B'", "D'"), ("C'", "D'"), ("D", "D'")]
    return pts, edges


def colourings(names, edges, k):
    idx = {v: i for i, v in enumerate(names)}
    E = [(idx[a], idx[b]) for a, b in edges]
    for col in itertools.product(range(k), repeat=len(names)):
        if col[0] != 0:                      # colour of O fixed to 0 (symmetry); every colouring has a relabelling
            continue
        if all(col[a] != col[b] for a, b in E):
            yield col


def main():
    ok = True
    one = num(1)
    cos, sin = num(F(5, 6)), num(c=F(1, 6))                      # cos theta = 5/6, sin theta = sqrt(11)/6
    need_unit = add(mul(cos, cos), mul(sin, sin)) == one
    pts, edges = spindle(cos, sin)
    names = list(pts)
    distinct = len({pts[v] for v in names}) == 7
    unit_pairs = [(a, b) for a, b in itertools.combinations(names, 2) if sqdist(pts[a], pts[b]) == one]
    edges_unit = all(sqdist(pts[a], pts[b]) == one for a, b in edges)
    same = sorted(map(sorted, unit_pairs)) == sorted(map(sorted, edges))
    c3 = next(colourings(names, edges, 3), None)
    c4 = next(colourings(names, edges, 4), None)
    good = need_unit and distinct and edges_unit and same and c3 is None and c4 is not None
    def fmt(x):
        parts = [f"{v}{s}" for v, s in zip(x, ("", " sqrt3", " sqrt11", " sqrt33")) if v]
        return " + ".join(parts).replace("+ -", "- ") or "0"
    print("Moser spindle in Q(sqrt 3, sqrt 11), cos theta = 5/6, sin theta = sqrt(11)/6 (cos^2 + sin^2 = 1: "
          f"{need_unit}):")
    for v in names:
        print(f"  {v:3} = ({fmt(pts[v][0])}, {fmt(pts[v][1])})")
    print(f"  7 distinct points: {distinct}; the 11 spindle edges have squared length 1: {edges_unit}; pairs at distance "
          f"1 among all 21: {len(unit_pairs)}, exactly the spindle's edges: {same}")
    print(f"  proper 3-colourings (the 729 assignments with O coloured 0, which loses nothing by symmetry): "
          f"{'none' if c3 is None else c3}; a proper "
          f"4-colouring: {dict(zip(names, c4)) if c4 else 'NONE'}")
    print(f"RESULT [MM1961]: the Moser spindle is a 7-vertex unit-distance graph with chromatic number 4, so "
          f"chi(R^2) >= 4: {'HOLDS' if good else 'FAILS'}")
    ok &= good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        p2, e2 = spindle(num(F(4, 5)), num(F(3, 5)))
        d = sqdist(p2["D"], p2["D'"])
        print(f"  cos theta = 4/5: |D D'|^2 = {fmt(d)}")
        res.append(("cos theta = 4/5", d != one))
        e3 = [e for e in edges if e != ("D", "D'")]
        c3b = next(colourings(names, e3, 3), None)
        print(f"  without the edge D D': a 3-colouring {dict(zip(names, c3b)) if c3b else 'NONE'}")
        res.append(("without D D', not 3-colourable", c3b is not None))
        res.append(("that 3-colouring on the spindle", c3b is not None and c3b[names.index("D")] == c3b[names.index("D'")]))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
