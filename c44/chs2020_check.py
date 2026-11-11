#!/usr/bin/env python3
"""Exact replay of the two lower-bound rows of the 44a page (teorth/optimizationproblems, main 87b4ebb), where
C_44 = sup_d R_d / 2^d and R_d is the largest number of relevant variables of a Boolean function of degree d:

  [CHS2020] 1     a complete read-once decision tree of depth d has degree d and 2^d - 1 relevant variables;
  [CHS2020] 3/2   the functions Xi_d of Chiarelli, Hatami and Saks (arXiv 1801.08564, section 5; also Shinkar and Tal):
                  Xi_1(x) = x and Xi_d(s, t, x, y) = (s+t)/2 Xi_{d-1}(x) + (s-t)/2 Xi_{d-1}(y) on l(d) = 2 l(d-1) + 2
                  variables, so l(d) = 3 * 2^(d-1) - 2, of degree d.

Over {-1,1} (degree and relevance are the same as over {0,1}).  For d = 1, ..., DMAX the multilinear polynomial is
built exactly by the recursion (products of polynomials on disjoint variables stay multilinear; rational
coefficients), its degree and the set of variables that occur are read off, and the function is evaluated on random
points against the definition (T: follow the tree; Xi_d: s Xi_{d-1}(x) if s = t, else s Xi_{d-1}(y)), which is
+-1-valued by construction.  For the small cases (Xi_d, d <= 3; trees, d <= 4) the whole truth table is transformed
(fast Walsh-Hadamard, integers) and must give exactly the same polynomial.  Then R_d / 2^d >= (2^d - 1)/2^d and
(3 * 2^(d-1) - 2)/2^d, increasing to 1 and 3/2.  Written from the page's rows and the paper's recursion; no download.
usage: python3 -I chs2020_check.py [--controls]
"""
import random
import sys
from fractions import Fraction

CONTROLS = "--controls" in sys.argv
DMAX = 8
rng = random.Random(44)


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def mul(p, q):
    out = {}
    for s1, c1 in p.items():
        for s2, c2 in q.items():
            need(not (s1 & s2), "products must be on disjoint variables")
            k = s1 | s2
            out[k] = out.get(k, 0) + c1 * c2
    return {k: c for k, c in out.items() if c}


def add(p, q):
    out = dict(p)
    for k, c in q.items():
        out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}


def half(sign_s, s, sign_t, t):
    """(sign_s * x_s + sign_t * x_t) / 2 as a polynomial."""
    return {frozenset([s]): Fraction(sign_s, 2), frozenset([t]): Fraction(sign_t, 2)}


class Counter:
    def __init__(self):
        self.n = 0

    def new(self):
        self.n += 1
        return self.n - 1


def tree(d, ctr):
    """(polynomial, evaluator) of the complete read-once decision tree of depth d: T_1(r) = r, and
    T_d(r, L, R) = (1+r)/2 T_{d-1}(L) + (1-r)/2 T_{d-1}(R)."""
    r = ctr.new()
    if d == 1:
        return {frozenset([r]): Fraction(1)}, (lambda x, r=r: x[r])
    pl, fl = tree(d - 1, ctr)
    pr, fr = tree(d - 1, ctr)
    poly = add(mul({frozenset(): Fraction(1, 2), frozenset([r]): Fraction(1, 2)}, pl),
               mul({frozenset(): Fraction(1, 2), frozenset([r]): Fraction(-1, 2)}, pr))
    return poly, (lambda x, r=r, fl=fl, fr=fr: fl(x) if x[r] == 1 else fr(x))


def xi(d, ctr, bad=False):
    """(polynomial, evaluator) of Xi_d.  bad: the control, (s+t)/2 used for both halves."""
    if d == 1:
        v = ctr.new()
        return {frozenset([v]): Fraction(1)}, (lambda x, v=v: x[v])
    s, t = ctr.new(), ctr.new()
    px, fx = xi(d - 1, ctr, bad)
    py, fy = xi(d - 1, ctr, bad)
    poly = add(mul(half(1, s, 1, t), px), mul(half(1, s, 1 if bad else -1, t), py))
    return poly, (lambda x, s=s, t=t, fx=fx, fy=fy: x[s] * fx(x) if x[s] == x[t] else x[s] * fy(x))


def wht_poly(f, n):
    """The multilinear polynomial of f from its truth table (bit i set means x_i = -1)."""
    vals = []
    for m in range(1 << n):
        x = [(-1 if (m >> i) & 1 else 1) for i in range(n)]
        vals.append(f(x))
    v = list(vals)
    h = 1
    while h < len(v):
        for i in range(0, len(v), 2 * h):
            for j in range(i, i + h):
                a, b = v[j], v[j + h]
                v[j], v[j + h] = a + b, a - b
        h *= 2
    return {frozenset(i for i in range(n) if (S >> i) & 1): Fraction(c, 1 << n) for S, c in enumerate(v) if c}, vals


def evaluate(poly, x):
    total = Fraction(0)
    for k, c in poly.items():
        term = c
        for i in k:
            term *= x[i]
        total += term
    return total


def check(name, build, d, small):
    ctr = Counter()
    poly, f = build(d, ctr)
    n = ctr.n
    deg = max(len(k) for k in poly)
    relevant = set().union(*poly.keys())
    for _ in range(200):
        x = [rng.choice((-1, 1)) for _ in range(n)]
        need(evaluate(poly, x) == f(x), f"{name} d = {d}: the polynomial differs from the definition")
    note = ""
    if n <= small:
        p2, vals = wht_poly(f, n)
        need(all(v in (-1, 1) for v in vals), f"{name} d = {d}: not +-1-valued")
        need(p2 == poly, f"{name} d = {d}: the truth table's transform differs from the recursion's polynomial")
        note = f", truth table of 2^{n} points: +-1-valued, its transform equals the polynomial"
    return n, deg, len(relevant), len(poly), note


def main():
    ok = True
    print("complete read-once decision trees (the [CHS2020] row 1):")
    for d in range(1, DMAX + 1):
        n, deg, rel, mons, note = check("tree", tree, d, 15)
        good = deg == d and rel == n == 2 ** d - 1
        ok &= good
        print(f"  d = {d}: {n} variables, degree {deg}, {rel} relevant, {mons} monomials{note}; R_d/2^d >= "
              f"{Fraction(rel, 2 ** d)}{'' if good else '  FAILS'}")
    print(f"RESULT [CHS2020] row 1: R_d >= 2^d - 1 for d = 1..{DMAX}, so C_44 >= sup_d (2^d - 1)/2^d = 1: "
          f"{'HOLDS' if ok else 'FAILS'}")
    ok2 = True
    print("Xi_d (the [CHS2020] row 3/2):")
    for d in range(1, DMAX + 1):
        n, deg, rel, mons, note = check("Xi", xi, d, 10)
        good = deg == d and rel == n == 3 * 2 ** (d - 1) - 2
        ok2 &= good
        print(f"  d = {d}: {n} variables, degree {deg}, {rel} relevant, {mons} monomials{note}; R_d/2^d >= "
              f"{Fraction(rel, 2 ** d)} = {rel / 2 ** d:.6f}{'' if good else '  FAILS'}")
    print(f"RESULT [CHS2020] row 3/2: R_d >= 3 * 2^(d-1) - 2 for d = 1..{DMAX}, so C_44 >= sup_d (3/2 - 2^(1-d)) = 3/2: "
          f"{'HOLDS' if ok2 else 'FAILS'}")
    ok &= ok2
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        ctr = Counter()
        poly, f = xi(3, ctr, bad=True)
        vals = {evaluate(poly, [rng.choice((-1, 1)) for _ in range(ctr.n)]) for _ in range(300)}
        res.append(("Xi_3 with (s+t)/2 in both halves", not vals <= {-1, 1}))
        print(f"  Xi_3 with (s+t)/2 in both halves: values {sorted(vals)} (not +-1-valued)")
        res.append(("a claim of 3 * 2^(d-1) - 1 variables at d = 4", 3 * 2 ** 3 - 1 != check("Xi", xi, 4, 0)[2]))
        res.append(("a claim of degree d - 1 for the tree at d = 3", check("tree", tree, 3, 0)[1] != 2))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
