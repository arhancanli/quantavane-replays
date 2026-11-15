#!/usr/bin/env python3
"""52a, program A: the first-moment upper bound [FP1983] 5.191 for the random 3-SAT threshold, from the page's model.

The page's model: a formula of density r = m/n is m clauses chosen uniformly at random, with replacement, from the
8 C(n,3) clauses on three distinct variables of n. For a fixed assignment, exactly one of the 8 sign patterns on each
triple is violated, so a random clause is satisfied with probability 7/8 and E[number of satisfying assignments] =
2^n (7/8)^m. Hence P(satisfiable) <= (2 (7/8)^r)^n -> 0 for r > r1 = log 2 / log(8/7), so r_3 <= r1.

This program (route A):
  * exhaustive: over ALL (8 C(n,3))^m formulas, the total number of (formula, satisfying assignment) pairs equals
    (8 C(n,3))^m 2^n (7/8)^m exactly, for (n, m) = (3, 1..3), (4, 1..3), (5, 1..2);
  * r1 < 5.191 in integers: 2 (7/8)^(5191/1000) < 1  <=>  7^5191 < 2^14573; and r1 > 5.190: 7^5190 > 2^14570.
A value from its definition (the row's own first-moment argument). Standard library only.
Usage: c52a_a.py [--claim R] [--variant two-literal]
"""
import sys
from fractions import Fraction as Q
from itertools import combinations, product

def clauses(n, variant):
    if variant == "two-literal":             # control: clauses on two distinct variables
        return [tuple(zip(pr, signs)) for pr in combinations(range(n), 2) for signs in product((0, 1), repeat=2)]
    return [tuple(zip(tri, signs)) for tri in combinations(range(n), 3) for signs in product((0, 1), repeat=3)]

def satisfied(cl, a):                         # literal (v, s) is true when a[v] == s
    return any(a[v] == s for v, s in cl)

def exhaustive(n, m, variant):
    C = clauses(n, variant)
    A = list(product((0, 1), repeat=n))
    sat = [[satisfied(c, a) for a in A] for c in C]
    total = 0
    for f in product(range(len(C)), repeat=m):
        for ai in range(len(A)):
            if all(sat[ci][ai] for ci in f):
                total += 1
    return total, len(C)

def main():
    a = sys.argv[1:]
    opt = dict(zip(a[0::2], a[1::2]))
    claim = Q(opt.get("--claim", "5.191"))
    variant = opt.get("--variant")
    if variant:
        print("CONTROL: clause model", variant)
    print("52a program A: [FP1983] r_3 <= %s by the first moment" % claim)
    ok = True
    for n, m in ((3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (4, 3), (5, 1), (5, 2)):
        total, nc = exhaustive(n, m, variant)
        want = Q(nc) ** m * 2 ** n * Q(7, 8) ** m
        good = total == want
        ok &= good
        print("n = %d, m = %d: %d formulas, %d (formula, satisfying assignment) pairs; (8 C(n,3))^m 2^n (7/8)^m = %s: %s"
              % (n, m, nc ** m, total, want, good))
        if not good:
            break
    if not ok:
        print("VERDICT REFUSED: E[#satisfying assignments] is not 2^n (7/8)^m in this model"); return
    # 2 (7/8)^claim < 1  <=>  (8/7)^claim > 2  <=>  8^p > 2^q 7^p  for claim = p/q
    p, q = claim.numerator, claim.denominator
    below = 8 ** p > 2 ** q * 7 ** p
    print("2 (7/8)^%s < 1 (exact integers: 8^%d > 2^%d 7^%d): %s" % (claim, p, q, p, below))
    lo_ok = 8 ** 519 < 2 ** 100 * 7 ** 519            # r1 > 5.19
    print("information: r1 = log 2/log(8/7) > 5.19 (8^519 < 2^100 7^519): %s; r1 < 5.191 iff 7^5191 < 2^14573: %s"
          % (lo_ok, 7 ** 5191 < 2 ** 14573))
    if below:
        print("row [FP1983] r_3 <= %s: True" % claim)
        print("VERDICT PASS")
    else:
        print("VERDICT REFUSED: 2 (7/8)^%s >= 1" % claim)

if __name__ == "__main__":
    main()
