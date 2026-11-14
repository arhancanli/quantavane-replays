#!/usr/bin/env python3
"""Program B for the five elementary rows of celem/ (20b cube and simplices, 25a cube, 39a cube, 11a n = 1), written
separately from celem_a.py and, where there is a second route, by it; exact rational arithmetic, standard library.

  20b cube: the moments of [-1/2, 1/2]^n by iterated exact integration of monomials over a box, Cov built entry by
       entry, det by Bareiss's fraction-free elimination: 12^-n for n <= 12.
  20b simplices: the moments of Delta_n by peeling one coordinate at a time, int_{Delta_n} x_1^a x_2^b =
       B(a+1, n+b) int_{Delta_(n-1)} x_2^b, each Beta integral by expanding (1 - x)^m binomially and integrating term by
       term (not the factorial formula); L^(2n) against (n!)^2 / ((n+1)^(n+1) (n+2)^n) for n <= 24.
  25a cube: Vol of the cross-polytope from its Ehrhart polynomial: the lattice-point counts of t B_1^n for
       t = 0..n by a recursion over the coordinates, the polynomial by Lagrange interpolation, its leading
       coefficient the volume: 2^n/n!, so n! 2^n Vol(K°) = 4^n for n <= 12.
  39a cube: for every pair of vertices p, q, the open boxes p + (-1, 1)^3 and q + (-1, 1)^3 (the possible positions v
       of a translate v + int K holding p, resp. q) are disjoint: their projections on some axis do not overlap.
  11a: E|f - E f| and E|grad f| by enumeration of {-1, 1}^n, for n = 1 and f(x) = x_1.
Red controls (--control NAME) as in celem_a.py: cube20b, simplex, polar, vertex, poincare.
"""
import itertools
import math
import sys
from fractions import Fraction

ok = True


def check(desc, cond, detail=""):
    global ok
    ok &= bool(cond)
    print("B: check %s: %s%s" % (desc, "PASS" if cond else "FAIL", (" (%s)" % detail) if detail else ""))
    return cond


def bareiss_det(M):
    """Determinant of a matrix of Fractions: scale rows to integers, then Bareiss's fraction-free elimination."""
    n = len(M)
    den = 1
    A = []
    for row in M:
        l = 1
        for x in row:
            l = l * x.denominator // math.gcd(l, x.denominator)
        den *= l
        A.append([int(x * l) for x in row])
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((r for r in range(k + 1, n) if A[r][k] != 0), None)
            if sw is None:
                return Fraction(0)
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return Fraction(sign * A[n - 1][n - 1], den)


def box_moment(exps, lo, hi):
    """int over the box [lo, hi]^n of prod x_i^e_i, iterated: one coordinate at a time."""
    r = Fraction(1)
    for e in exps:
        r *= (hi ** (e + 1) - lo ** (e + 1)) / (e + 1)
    return r


def beta_int(a, m):
    """int_0^1 x^a (1 - x)^m dx by expanding (1 - x)^m."""
    return sum(Fraction((-1) ** k * math.comb(m, k), a + k + 1) for k in range(m + 1))


def simplex_moment(n, a, b):
    """int over Delta_n (n >= 1) of x_1^a x_2^b (b = 0 allowed; for n = 1 only x_1^a)."""
    if n == 1:
        return beta_int(a + b, 0) if b == 0 else None
    # x_1 = s, the rest in (1 - s) Delta_(n-1): x_2^b scales by (1 - s)^b, the volume element by (1 - s)^(n-1)
    inner = simplex_moment(n - 1, b, 0)
    return beta_int(a, n - 1 + b) * inner


def ehrhart_count(n, t):
    """Number of integer points with sum |x_i| <= t in Z^n, by a recursion over the coordinates."""
    cnt = [1] + [0] * t                                    # cnt[s] = points of the first k coordinates with sum s
    for _ in range(n):
        new = [0] * (t + 1)
        for s, c in enumerate(cnt):
            if c:
                for x in range(-(t - s), t - s + 1):
                    new[s + abs(x)] += c
        cnt = new
    return sum(cnt)


def lead_coeff(xs, ys):
    """Leading coefficient of the interpolating polynomial of degree len(xs) - 1."""
    return sum(Fraction(y) / math.prod(Fraction(x - z) for z in xs if z != x) for x, y in zip(xs, ys))


def main():
    ctl = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
    if ctl:
        print("B: CONTROL %s" % ctl)
    # 20b cube
    side = Fraction(2) if ctl == "cube20b" else Fraction(1)
    lo, hi = -side / 2, side / 2
    good = True
    for n in range(1, 13):
        vol = Fraction(1)                                  # the row's volume 1 (the control changes the side only)
        def E(exps):
            return box_moment(exps, lo, hi) / box_moment([0] * n, lo, hi)
        mean = [E([1 if k == i else 0 for k in range(n)]) for i in range(n)]
        C = [[E([(k == i) + (k == j) for k in range(n)]) - mean[i] * mean[j] for j in range(n)] for i in range(n)]
        good &= bareiss_det(C) / vol ** 2 == Fraction(1, 12 ** n)
    check("20b cube: det Cov / Vol^2 = 12^-n for n = 1..12 (iterated integrals, Bareiss)", good)
    # 20b simplices
    good = True
    for n in range(2, 25):
        vol = simplex_moment(n, 0, 0)
        m1 = simplex_moment(n, 1, 0) / vol
        v11 = simplex_moment(n, 2, 0) / vol - m1 * m1
        v12 = simplex_moment(n, 1, 1) / vol - m1 * m1
        if ctl == "simplex" and n == 5:
            v11 += Fraction(1, 1000)
        C = [[v11 if i == j else v12 for j in range(n)] for i in range(n)]
        good &= vol == Fraction(1, math.factorial(n))
        good &= bareiss_det(C) / vol ** 2 == Fraction(math.factorial(n) ** 2, (n + 1) ** (n + 1) * (n + 2) ** n)
    check("20b simplex: Vol = 1/n! and det Cov / Vol^2 = (n!)^2 / ((n+1)^(n+1) (n+2)^n) for n = 2..24 (Beta "
          "integrals by binomial expansion)", good)
    # 25a cube
    good = True
    for n in range(1, 13):
        ts = list(range(n + 1))
        vol = lead_coeff(ts, [ehrhart_count(n, t) for t in ts])
        if ctl == "polar":
            vol = Fraction(2 ** n)                         # the cube taken as its own polar
        good &= vol == Fraction(2 ** n, math.factorial(n)) and math.factorial(n) * 2 ** n * vol == 4 ** n
    check("25a cube: Vol(B_1^n) from the Ehrhart polynomial is 2^n/n!, and n! 2^n Vol = 4^n for n = 1..12", good)
    # 39a cube
    verts = list(itertools.product((-1, 1), repeat=3))
    if ctl == "vertex":
        verts = [v[:3] for v in itertools.product((-1, 1), repeat=4)]
    disjoint = all(any(max(p[i], q[i]) - 1 >= min(p[i], q[i]) + 1 for i in range(3))
                   for p, q in itertools.combinations(verts, 2))
    check("39a cube: for every pair of the %d points the sets of translates v with p, resp. q in v + int K are "
          "disjoint (a projection does not overlap)" % len(verts), disjoint and len(set(verts)) == len(verts) == 8)
    # 11a
    n = 1
    pts = list(itertools.product((-1, 1), repeat=n))
    f = {x: Fraction(x[0]) for x in pts}
    Ef = sum(f.values()) / len(pts)
    dj = 1 if ctl == "poincare" else 2                       # the control drops the 1/2 in D_j f
    lhs = sum(abs(f[x] - Ef) for x in pts) / len(pts)
    grad = []
    for x in pts:
        sq = Fraction(0)
        for j in range(n):
            xj = tuple(-v if i == j else v for i, v in enumerate(x))
            sq += ((f[x] - f[xj]) / dj) ** 2
        grad.append(sq)                                    # |grad f|^2; here a square of a rational
    rg = [Fraction(math.isqrt(g.numerator), math.isqrt(g.denominator)) for g in grad]
    exact = all(r * r == g for r, g in zip(rg, grad))
    rhs = sum(rg) / len(pts)
    check("11a: n = 1, f(x) = x: E|f - Ef| = %s and E|grad f| = %s by enumeration" % (lhs, rhs),
          exact and lhs == 1 and rhs == 1)
    print("B: VERDICT %s" % ("PASS" if ok else "REFUSED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
