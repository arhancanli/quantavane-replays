#!/usr/bin/env python3
"""Exact replay of the 4a lower-bound row [P1970] (teorth/optimizationproblems, main 87b4ebb): the cap set constant
satisfies C_4 >= 2.1146, from a cap of 20 points in F_3^4 (Pellegrino 1970: the largest cap in four dimensions has 20
points).  The page: "Any cap set A in F_3^n gives a lower bound C_4 >= |A|^(1/n) by taking Cartesian powers", so a
20-point cap in F_3^4 gives C_4 >= 20^(1/4) = 2.11474252688...

The paper is not accessible to us, so the cap is found here by our own search: depth-first search in F_3^4 from the
affine frame {0, e_1, e_2, e_3, e_4} (every cap spanning F_3^4 contains five affinely independent points, and an
affine map moves them there), each new point the least admissible one in lexicographic order, with forbidden points
-(x + y) for every pair.  Then, with no reference to the search:
  1. route a: all C(20,3) = 1140 triples of distinct points, none with x + y + z = 0 (the lines of AG(n,3) are exactly
     the triples of distinct points with sum 0);
  2. route b: for each of the 190 pairs the third point of its line, -(x + y), is not in the set;
  3. the Cartesian square: A x A, 400 points of F_3^8, is again a cap (route b on its 79,800 pairs);
  4. completeness, as an anchor: each of the other 61 points lies on a line through two points of the cap;
  5. 20^(1/4) > 2.1146 exactly (2.1146^4 < 20) and 2.1147^4 < 20 < 2.1148^4, so 2.1147 would be the four-decimal
     truncation of 20^(1/4); the printed 2.1146 is a valid, slightly lower, bound.
Red controls (--controls): the cap with one point replaced by a point outside it has a line; the cap with a 21st point
has a line (completeness); the product of the cap with a set that is not a cap is not a cap; 2.1148 is refused.
usage: python3 -I cap_check.py [--controls]
"""
import itertools
import sys
from fractions import Fraction

CONTROLS = "--controls" in sys.argv
N = 4


def neg_sum(x, y):
    return tuple((-a - b) % 3 for a, b in zip(x, y))


def search(n, target):
    pts = list(itertools.product(range(3), repeat=n))
    frame = [tuple(0 for _ in range(n))] + [tuple(1 if i == j else 0 for i in range(n)) for j in range(n)]
    best = []

    def rec(cap, forb, start):
        if len(cap) == target:
            best.append(list(cap))
            return True
        for i in range(start, len(pts)):
            p = pts[i]
            if p in forb or p in cap:
                continue
            new = {neg_sum(p, q) for q in cap}
            if p in new:
                continue
            cap.append(p)
            if rec(cap, forb | new, i + 1):
                return True
            cap.pop()
        return False

    forb = {neg_sum(x, y) for x, y in itertools.combinations(frame, 2)}
    rec(list(frame), forb, 0)
    return best[0] if best else None


def has_line_triples(A):
    for x, y, z in itertools.combinations(A, 3):
        if all((a + b + c) % 3 == 0 for a, b, c in zip(x, y, z)):
            return (x, y, z)
    return None


def has_line_pairs(A):
    S = set(A)
    for x, y in itertools.combinations(A, 2):
        z = neg_sum(x, y)
        if z in S and z != x and z != y:
            return (x, y, z)
    return None


def main():
    ok = True
    A = search(N, 20)
    need = A is not None and len(set(A)) == 20
    print(f"a 20-point cap in F_3^4 found by our search: {A is not None}")
    print("  " + " ".join("".join(map(str, p)) for p in A))
    la, lb = has_line_triples(A), has_line_pairs(A)
    prod = [x + y for x in A for y in A]
    lp = has_line_pairs(prod)
    outside = [p for p in itertools.product(range(3), repeat=N) if p not in set(A)]
    complete = all(any(neg_sum(x, y) == p for x, y in itertools.combinations(A, 2)) for p in outside)
    exact = Fraction(21146, 10000) ** 4 < 20 and Fraction(21147, 10000) ** 4 < 20 < Fraction(21148, 10000) ** 4
    print(f"  route a (1140 triples): {'no line' if la is None else 'LINE ' + str(la)}; route b (190 pairs): "
          f"{'no line' if lb is None else 'LINE ' + str(lb)}")
    print(f"  A x A, 400 points of F_3^8 (79,800 pairs): {'a cap' if lp is None else 'LINE'}; complete (each of the 61 "
          f"other points on a line through two cap points): {complete}")
    print(f"  20^(1/4) = {20 ** 0.25:.12f}; 2.1146^4 < 20: {Fraction(21146, 10000) ** 4 < 20}; 2.1147^4 < 20 < "
          f"2.1148^4: {Fraction(21147, 10000) ** 4 < 20 < Fraction(21148, 10000) ** 4}")
    good = need and la is None and lb is None and lp is None and exact
    print(f"RESULT [P1970]: a 20-point cap in F_3^4 exists, so C_4 >= 20^(1/4) = 2.114742... > 2.1146: "
          f"{'HOLDS' if good else 'FAILS'}")
    ok &= good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        B = A[:-1] + [outside[0]]
        res.append(("one point replaced", has_line_pairs(B) is not None))
        C = A + [outside[-1]]
        res.append(("a 21st point", has_line_pairs(C) is not None and has_line_triples(C) is not None))
        D = A[:2] + [neg_sum(A[0], A[1])]
        res.append(("product with a line", has_line_pairs([x + y for x in A for y in D]) is not None))
        res.append(("2.1148 claimed", not (Fraction(21148, 10000) ** 4 < 20)))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
