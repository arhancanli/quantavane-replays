#!/usr/bin/env python3
"""Exact replay of two 9a upper-bound rows (teorth/optimizationproblems, main 87b4ebb) on the Shannon capacity of C_7:

  [L1979]  Theta(C_7) <= theta(C_7) ~ 3.3177 (Lovasz).  Lovasz's bound for a regular graph G on n vertices with adjacency
           eigenvalues lambda_1 >= ... >= lambda_n: theta(G) <= -n lambda_n / (lambda_1 - lambda_n), the value of the
           feasible matrix B = J - t A, t = n/(lambda_1 - lambda_n) (B_ij = 1 on the diagonal and on non-edges), whose
           largest eigenvalue is -n lambda_n/(lambda_1 - lambda_n).  For C_7 this is 7 cos(pi/7)/(1 + cos(pi/7)).
  [S1956]  Theta(C_7) <= 7/2 (Shannon): the fractional clique cover (weight 1/2 on each of the 7 edges) has value 7/2;
           the fractional packing (1/2 on each vertex) shows the LP optimum alpha*(C_7) is exactly 7/2.

Computed (integers and rationals; no download):
  1. the characteristic polynomial of A(C_7) by Faddeev-LeVerrier in exact rationals: (x - 2)(x^3 + x^2 - 2x - 1)^2;
  2. the least root lambda_7 of x^3 + x^2 - 2x - 1, isolated by bisection over the rationals to width 10^-40 (one sign
     change; the other two roots are bracketed too, so lambda_7 is the least eigenvalue), and -2 cos(pi/7) inside it;
  3. the bound -7 lambda_7/(2 - lambda_7) as an interval (it decreases in lambda_7) and the printed 3.3177 against it;
  4. the certificate: the eigenvalues of vI - B (v the bound) are v - 7 + 2t = 0 on the all-ones vector and v + t mu
     >= v + t lambda_7 = 0 on the others (mu the other eigenvalues of A), so vI - B is positive semidefinite; also
     checked in floating point on the circulant's eigenvalues 2 cos(2 pi k/7);
  5. [S1956]: the cover and the packing in exact rationals (every clique of C_7 is a vertex or an edge).
Red controls (--controls): t = 1 instead of 7/(2 - lambda_7) gives a larger lambda_max(B) (so a worse bound); the
claimed 3.3176 is below the exact value (refused as an upper bound); weights 1/3 on the edges do not cover; the
eigenvalue route with lambda_7 replaced by the second root gives a matrix vI - B with a negative eigenvalue.
usage: python3 -I theta_check.py [--controls]
"""
import math
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv
N = 7


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def matmul(X, Y):
    n = len(X)
    return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def charpoly(A):
    """Faddeev-LeVerrier: coefficients c_n = 1, c_(n-1), ..., c_0 of det(xI - A)."""
    n = len(A)
    M = [[F(0)] * n for _ in range(n)]
    I = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    coeffs = [F(1)]
    c = F(1)
    for k in range(1, n + 1):
        AM = matmul(A, M)
        M = [[AM[i][j] + c * I[i][j] for j in range(n)] for i in range(n)]
        AM = matmul(A, M)
        c = -sum(AM[i][i] for i in range(n)) / k
        coeffs.append(c)
    return coeffs                     # highest degree first


def polymul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def ev(p, x):
    v = 0
    for a in p:
        v = v * x + a
    return v


def bisect(p, lo, hi, width):
    flo = ev(p, lo)
    need(flo * ev(p, hi) < 0, "no sign change")
    while hi - lo > width:
        mid = (lo + hi) / 2
        fm = ev(p, mid)
        if fm == 0:
            return mid, mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return lo, hi


def main():
    ok = True
    A = [[F(int((j - i) % N in (1, N - 1))) for j in range(N)] for i in range(N)]
    cp = charpoly(A)
    cubic = [1, 1, -2, -1]
    expect = polymul([1, -2], polymul(cubic, cubic))
    cp_ok = [int(c) for c in cp] == expect and all(c.denominator == 1 for c in cp)
    W = F(1, 10 ** 40)
    r1 = bisect(cubic, F(-2), F(-1), W)          # least root
    r2 = bisect(cubic, F(-1), F(0), W)
    r3 = bisect(cubic, F(1), F(2), W)
    lam_lo, lam_hi = r1
    cospi7 = math.cos(math.pi / 7)
    inside = float(lam_lo) - 1e-12 <= -2 * cospi7 <= float(lam_hi) + 1e-12
    f = lambda lam: -7 * lam / (2 - lam)
    v_lo, v_hi = f(lam_hi), f(lam_lo)            # decreasing in lambda
    printed = F(33177, 10000)
    holds = printed >= v_hi
    # certificate in floating point on the circulant spectrum
    lam = float(lam_lo)
    t = 7 / (2 - lam)
    v = float(v_lo)
    spec_A = [2 * math.cos(2 * math.pi * k / 7) for k in range(7)]
    spec = [v - 7 + 2 * t] + [v + t * mu for mu in spec_A[1:]]
    psd = all(s > -1e-12 for s in spec)
    good = cp_ok and r1[1] < r2[0] and inside and holds and psd
    print(f"[L1979] characteristic polynomial of A(C_7) (Faddeev-LeVerrier, exact): {[int(c) for c in cp]} = "
          f"(x - 2)(x^3 + x^2 - 2x - 1)^2: {cp_ok}")
    print(f"  roots of x^3 + x^2 - 2x - 1 isolated in (-2,-1), (-1,0), (1,2); least eigenvalue lambda_7 in "
          f"[{float(lam_lo):.15f}, {float(lam_hi):.15f}] (width 10^-40), contains -2 cos(pi/7): {inside}")
    print(f"  theta(C_7) <= -7 lambda_7/(2 - lambda_7) in [{float(v_lo):.15f}, {float(v_hi):.15f}]; 7 cos(pi/7)/(1 + "
          f"cos(pi/7)) = {7 * cospi7 / (1 + cospi7):.15f}; printed 3.3177 >= upper end: {holds} (rounded up)")
    print(f"  certificate B = J - tA, t = 7/(2 - lambda_7) = {t:.12f}: spectrum of vI - B (float) min "
          f"{min(spec):.3e}, positive semidefinite: {psd}")
    print(f"RESULT [L1979]: Theta(C_7) <= theta(C_7) <= {float(v_hi):.13f}... <= 3.3177: {'HOLDS' if good else 'FAILS'}")
    ok &= good
    edges = [(i, (i + 1) % N) for i in range(N)]
    y = {e: F(1, 2) for e in edges}
    cover = all(sum(y[e] for e in edges if v_ in e) >= 1 for v_ in range(N))
    w = [F(1, 2)] * N
    pack = all(w[a] + w[b] <= 1 for a, b in edges) and all(x <= 1 for x in w)
    good2 = cover and pack and sum(y.values()) == F(7, 2) == sum(w)
    print(f"[S1956] cliques of C_7: 7 vertices and 7 edges; cover 1/2 on each edge covers every vertex: {cover}, value "
          f"{sum(y.values())}; packing 1/2 on each vertex: {pack}, value {sum(w)}; so alpha*(C_7) = 7/2")
    print(f"RESULT [S1956]: Theta(C_7) <= alpha*(C_7) = 7/2: {'HOLDS' if good2 else 'FAILS'}")
    ok &= good2
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        lmax_t1 = max([7 - 2] + [-mu for mu in spec_A[1:]])        # B = J - A: 7 - 2 on 1, -mu elsewhere
        res.append(("t = 1", lmax_t1 > float(v_hi)))
        print(f"  t = 1: lambda_max(J - A) = {lmax_t1:.6f} > {float(v_hi):.6f}")
        res.append(("a claimed 3.3176", not (F(33176, 10000) >= v_lo)))
        y2 = {e: F(1, 3) for e in edges}
        res.append(("1/3 on the edges", not all(sum(y2[e] for e in edges if v_ in e) >= 1 for v_ in range(N))))
        lam2 = float(r2[0])
        t2 = 7 / (2 - lam2)
        v2 = -7 * lam2 / (2 - lam2)
        spec2 = [v2 - 7 + 2 * t2] + [v2 + t2 * mu for mu in spec_A[1:]]
        res.append(("the second root for lambda_7", min(spec2) < -1e-9))
        print(f"  the second root: min eigenvalue of vI - B = {min(spec2):.4f}")
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
