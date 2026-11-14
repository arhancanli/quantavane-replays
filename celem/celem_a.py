#!/usr/bin/env python3
"""Program A for five elementary rows on pages that had no replay folder (20b, 25a, 39a, 11a), each from the object
the row names; exact rational arithmetic, standard library only.

  20b, "Classical (cube)", C_20b >= 1/sqrt(12): K = [-1/2, 1/2]^n, Cov(K) from the one-dimensional moments
       (int x dx = 0, int x^2 dx = 1/12 over [-1/2, 1/2], coordinates independent), Vol = 1, so
       L_K^(2n) = det Cov / Vol^2 = 12^-n for every n (checked for n <= 12).
  20b, "Classical (simplices)", C_20b >= 1/e: the standard simplex Delta_n = {x >= 0, sum x <= 1} of R^n, its moments by
       the Dirichlet formula  int x^a dx = prod a_i! / (n + |a|)!, Cov exactly, det by exact elimination, Vol = 1/n!;
       L^(2n) equals (n!)^2 / ((n+1)^(n+1) (n+2)^n) for n <= 24, is increasing in n, and stays below e^(-2n). That
       the sequence tends to 1/e (Stirling) is the limit step, not computed.
  25a, the cube, C_25 <= 4: K = [-1, 1]^n, K° = {y : sum |y_i| <= 1} (the support function of K is sum |y_i|), the
       cross-polytope; Vol(K°) = 2^n / n! from its 2^n orthant simplices, so n! Vol(K) Vol(K°) = 4^n (n <= 12).
  39a, the cube, H_3 >= 8: a translate v + int([-1, 1]^3) holds at most one vertex of the cube, since two vertices
       differ by 2 in some coordinate and |p_i - v_i| < 1, |q_i - v_i| < 1 would give |p_i - q_i| < 2.
  11a, n = 1, f(x) = x: E|f - E f| = 1 and |grad f| = |D_1 f| = |x| = 1, so the ratio is 1 and C_11a >= 1.
Red controls (--control NAME): cube20b (side 1 replaced by side 2 with the volume left at 1), simplex (one second
moment off by 1/1000), polar (the cube's polar taken as the cube itself), vertex (a fourth coordinate equal for
all vertices: two vertices of a 4-cube restricted to three coordinates), poincare (D_j f without its factor 1/2).
"""
import itertools
import math
import sys
from fractions import Fraction

ok = True


def check(desc, cond, detail=""):
    global ok
    ok &= bool(cond)
    print("A: check %s: %s%s" % (desc, "PASS" if cond else "FAIL", (" (%s)" % detail) if detail else ""))
    return cond


def det(M):
    M = [row[:] for row in M]
    n, d = len(M), Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= f * M[i][c]
    return d


def e_bounds(k):
    """(lo, hi) for e^k, k a positive integer, from the series of e with its tail bound."""
    s, t, m = Fraction(0), Fraction(1), 0
    while t > Fraction(1, 10 ** 40):
        s += t
        m += 1
        t /= m
    return s ** k, (s + 2 * t) ** k


def main():
    ctl = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
    if ctl:
        print("A: CONTROL %s" % ctl)
    # 20b, cube
    side = Fraction(2) if ctl == "cube20b" else Fraction(1)
    a, b = -side / 2, side / 2
    m1 = (b ** 2 - a ** 2) / 2 / (b - a)             # mean of the uniform law on [a, b]
    m2 = (b ** 3 - a ** 3) / 3 / (b - a)
    var = m2 - m1 ** 2
    good = all(det([[var if i == j else Fraction(0) for j in range(n)] for i in range(n)]) / Fraction(1) ** 2
               == Fraction(1, 12 ** n) for n in range(1, 13))
    check("20b cube: L^(2n) = det Cov / Vol^2 = 12^-n for n = 1..12, so L = 1/sqrt(12) in every dimension", good,
          "variance %s" % var)
    # 20b, simplex
    vals = []
    good = True
    for n in range(1, 25):
        F = math.factorial

        def mom(alpha):                               # int over Delta_n of x^alpha, divided by Vol = 1/n!
            return Fraction(F(n) * math.prod(F(x) for x in alpha), F(n + sum(alpha)))
        mean = [mom([1 if j == i else 0 for j in range(n)]) for i in range(n)]
        C = [[mom([(1 if k == i else 0) + (1 if k == j else 0) for k in range(n)]) - mean[i] * mean[j]
              for j in range(n)] for i in range(n)]
        if ctl == "simplex" and n == 5:
            C[0][0] += Fraction(1, 1000)
        L2n = det(C) / Fraction(1, F(n)) ** 2
        good &= L2n == Fraction(F(n) ** 2, (n + 1) ** (n + 1) * (n + 2) ** n)
        vals.append(L2n)
    inc = all(vals[n] ** n > vals[n - 1] ** (n + 1) for n in range(1, len(vals)))   # L_(n+1) > L_n (2n(n+1)-th powers)
    below = all(vals[n - 1] * e_bounds(2 * n)[1] < 1 for n in range(1, len(vals) + 1))
    check("20b simplex: det Cov / Vol^2 = (n!)^2 / ((n+1)^(n+1) (n+2)^n) for n = 1..24", good)
    check("20b simplex: L_n increasing and below 1/e for n = 1..24 (L_24 = %.9f, 1/e = 0.367879441)" % (
          float(vals[-1]) ** (1 / 48)), inc and below)
    # 25a, cube and cross-polytope
    good = True
    for n in range(1, 13):
        verts = list(itertools.product((-1, 1), repeat=n))
        # support function of K at y: max over vertices of <x, y>; at the 2n points +-e_i it is 1, so +-e_i lie on the
        # boundary of K°; the facets of conv(+-e_i) are sum s_i y_i <= 1, one for each vertex s of K
        supp_ok = all(max(sum(x[i] * (1 if i == k else 0) * sg for i in range(n)) for x in verts) == 1
                      for k in range(n) for sg in (1, -1))
        polar_is_cube = ctl == "polar"
        vol_polar = Fraction(2 ** n, 1) if polar_is_cube else Fraction(2 ** n, math.factorial(n))
        good &= supp_ok and math.factorial(n) * 2 ** n * vol_polar == 4 ** n
    check("25a cube: n! Vol(K) Vol(K°) = 4^n for n = 1..12 (K° the cross-polytope, 2^n simplices of volume 1/n!)",
          good)
    # 39a, cube
    verts = list(itertools.product((-1, 1), repeat=3))
    if ctl == "vertex":
        verts = [v[:3] for v in itertools.product((-1, 1), repeat=4)]   # 16 points, repeated vertices
    sep = all(any(abs(p[i] - q[i]) >= 2 for i in range(3)) for p, q in itertools.combinations(verts, 2))
    check("39a cube: any two vertices differ by 2 in some coordinate, so each translate of int(K) holds at most one of "
          "the %d vertices and H_3 >= 8" % len(verts), sep and len(verts) == 8)
    # 11a, n = 1, f(x) = x
    pts = [-1, 1]
    f = lambda x: x
    Ef = sum(Fraction(f(x)) for x in pts) / 2
    half = Fraction(1) if ctl == "poincare" else Fraction(1, 2)   # the control drops the 1/2 in D_j f
    lhs = sum(abs(f(x) - Ef) for x in pts) / Fraction(2)
    rhs = sum(abs((f(x) - f(-x)) * half) for x in pts) / Fraction(2)
    check("11a: n = 1, f(x) = x: E|f - Ef| = %s, E|grad f| = %s, ratio 1" % (lhs, rhs), lhs == 1 and rhs == 1)
    print("A: VERDICT %s" % ("PASS" if ok else "REFUSED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
