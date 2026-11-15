#!/usr/bin/env python3
"""83a, the search step (not part of any proof): an approximate positive eigenfunction psi of Wirsing's operator V.

Wirsing's operator (from the Gauss-Kuzmin operator T phi(x) = sum_k p_k(x) phi(1/(x+k)), p_k = (1+x)/((x+k)(x+k+1)),
by Abel summation; (T phi)' = -V(phi')):
    V psi(x) = sum_{k>=1} [ k/(x+k+1)^2 * int_{1/(x+k+1)}^{1/(x+k)} psi(t) dt + (1+x)/((x+k)^3 (x+k+1)) * psi(1/(x+k)) ].
This script looks for a polynomial psi of degree n with V psi ~ lambda psi on [0, 1] by collocation at n+1 Chebyshev
points and power iteration, in 60-digit decimal arithmetic WITHOUT error control, and writes psi.json (the
coefficients rounded to multiples of 2^-64, and the bounds to be proved). Both checking programs read psi.json and
prove the bounds from scratch; nothing computed here is trusted. Standard library only.
Usage: c83a_search.py [--degree N] [--out psi.json]
"""
import sys, json
from decimal import Decimal as D, getcontext
from fractions import Fraction as Q
from math import comb, cos, pi

getcontext().prec = 60

def bernoulli(n):
    B = [Q(0)] * (n + 1); B[0] = Q(1)
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, j) * B[j] for j in range(m)) / (m + 1)
    return B
BER = bernoulli(60)

def hurwitz(s, a, N=25, M=25):
    """zeta(s, a), Euler-Maclaurin, no error control (search only)"""
    tot = sum(D(1) / (a + k) ** s for k in range(N))
    aN = a + N
    tot += aN ** (1 - s) / (s - 1) + aN ** (-s) / 2
    poch = D(s)
    for j in range(1, M + 1):
        if j > 1:
            poch *= (s + 2 * j - 3) * (s + 2 * j - 2)
        b = BER[2 * j]
        tot += D(b.numerator) / D(b.denominator) / D(_fact(2 * j)) * poch * aN ** (-s - 2 * j + 1)
    return tot

def _fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r

def pf(a, b):
    """1/(y^a (y+1)^b) = sum A_i y^-i + sum B_i (y+1)^-i"""
    out = []
    if b == 0:
        return [(("y", a), 1)]
    if a == 0:
        return [(("z", b), 1)]
    for i in range(1, a + 1):
        out.append((("y", i), (-1) ** (a - i) * comb(a + b - i - 1, b - 1)))
    for i in range(1, b + 1):
        out.append((("z", i), (-1) ** a * comb(a + b - i - 1, a - 1)))
    return out

def terms_of_monomial(j):
    """V(u^j) as {('y'|'z', i): [alpha, beta]} meaning sum_k (alpha + beta X) y^-i, X = 1 + x, y = x + k, z = y + 1"""
    acc = {}
    def add(a, b, al, be):
        for key, co in pf(a, b):
            v = acc.setdefault(key, [Q(0), Q(0)])
            v[0] += al * co; v[1] += be * co
    d = Q(1, j + 1)
    add(j + 1, 1, d, 0)
    add(0, j + 2, -d, 0)
    add(j + 1, 2, 0, -d)
    add(0, j + 3, 0, d)
    add(j + 3, 1, 0, Q(1))
    return acc

def V_monomial_at(j, x, zcache):
    X = 1 + x
    acc = terms_of_monomial(j)
    y1 = acc.get(("y", 1), [0, 0]); z1 = acc.get(("z", 1), [0, 0])
    assert y1[0] + z1[0] == 0 and y1[1] + z1[1] == 0
    tot = (D(y1[0].numerator) / y1[0].denominator + D(y1[1].numerator) / y1[1].denominator * X) / X if y1 != [0, 0] else D(0)
    for (kind, i), (al, be) in acc.items():
        if i == 1:
            continue
        c = D(al.numerator) / al.denominator + D(be.numerator) / be.denominator * X
        zi = zcache.setdefault((i, x), hurwitz(i, X))
        if kind == "y":
            tot += c * zi
        else:
            tot += c * (zi - X ** (-i))
    return tot

def solve(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]

def main():
    n = 18; out = "psi.json"
    a = sys.argv[1:]
    for i, t in enumerate(a):
        if t == "--degree": n = int(a[i + 1])
        if t == "--out": out = a[i + 1]
    xs = [D(str((1 - cos(pi * (2 * i + 1) / (2 * (n + 1)))) / 2)) for i in range(n + 1)]
    zc = {}
    A = [[V_monomial_at(j, x, zc) for j in range(n + 1)] for x in xs]
    B = [[x ** j if j else D(1) for j in range(n + 1)] for x in xs]
    # C = B^-1 A column by column
    cols = [solve(B, [A[i][j] for i in range(n + 1)]) for j in range(n + 1)]
    C = [[cols[j][i] for j in range(n + 1)] for i in range(n + 1)]
    c = [D(1)] + [D(0)] * n
    lam = D(0)
    for it in range(200):
        nc = [sum(C[i][j] * c[j] for j in range(n + 1)) for i in range(n + 1)]
        s = sum(nc)            # psi(1)
        nl = s / sum(c)
        c = [v / s for v in nc]
        if it > 10 and abs(nl - lam) < D(10) ** -40:
            lam = nl; break
        lam = nl
    print("collocation eigenvalue estimate:", +lam)
    co = [Q(round(v * (1 << 64)), 1 << 64) for v in c]
    # non-rigorous look at the ratio on a grid
    def psi(t):
        r = D(0)
        for v in reversed(co):
            r = r * t + D(v.numerator) / v.denominator
        return r
    rs = []
    for g in range(0, 41):
        x = D(g) / 40
        zc2 = {}
        Vx = sum(D(v.numerator) / v.denominator * V_monomial_at(j, x, zc2) for j, v in enumerate(co))
        rs.append(Vx / psi(x))
    print("grid ratio min", +min(rs), "max", +max(rs))
    print("psi(0), psi(1/2), psi(1):", +psi(D(0)), +psi(D("0.5")), +psi(D(1)))
    print("max |coefficient|:", float(max(abs(v) for v in co)))
    from math import floor, ceil
    lo = Q(floor((Q(str(min(rs))) - Q(1, 10 ** 16)) * 10 ** 15), 10 ** 15)
    hi = Q(ceil((Q(str(max(rs))) + Q(1, 10 ** 16)) * 10 ** 15), 10 ** 15)
    data = {
        "page": "83a",
        "operator": "Wirsing's V on C[0,1]",
        "psi_monomial_coefficients": [str(v) for v in co],
        "claimed_lower": "%.15f" % lo,
        "claimed_upper": "%.15f" % hi,
        "note": "found by c83a_search.py (collocation, power iteration, no error control); proved by c83a_a.py and c83a_b.py",
    }
    with open(out, "w") as f:
        json.dump(data, f, indent=1)
    print("written", out.split("/")[-1], "claims", "%.15f" % lo, "%.15f" % hi)

if __name__ == "__main__":
    main()
