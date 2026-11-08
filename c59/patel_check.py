#!/usr/bin/env python3
"""Independent exact check of the 59a row [P2026], K_2 < 0.302825279492 (Patel; PR #169), from the construction the row
states alone: no file of the proof package was read or run.

  python3 -I patel_check.py [--controls]

The row: U = (1 + z)(1 - w), V = 1 + zw, P = L U + i T V, Q = L V + i T U, f = (S Q - P)/(S Q + P), with L = 2500000000,
T = 3067398171, S = 10^15; the 841 coefficients f_jk, 0 <= j, k <= N = 28, give sum |f_jk| r^(j+k) > 1 at
r = 302825279492/10^12.
  1. f is in the Schur class S_2 (re-derived here). Q conj(P) = L^2 V conj(U) + T^2 U conj(V) + i L T (|U|^2 - |V|^2), so
     Re(Q conj(P)) = (L^2 + T^2) Re(U conj(V)), and with a = |z|, b = |w|,
       Re(U conj(V)) = 1 - a^2 b^2 + Re(z)(1 - b^2) - Re(w)(1 - a^2) >= (1 - a)(1 - b)(1 - a b) > 0
     on the open bidisc. So Re(S Q conj(P)) > 0, hence |S Q - P| < |S Q + P| (S > 0): |f| < 1, and S Q + P has no zero.
     The two identities are checked symbolically below, the factorisation on a 3 x 3 grid (an identity, since
     both sides have degree at most 2 in a and in b), and the two
     inequalities exactly at 2,000 random rational points.
  2. The coefficients, exactly: with D = S Q + P (bidegree (1, 1), D_00 = (S + 1)(L + iT)) and Num = S Q - P,
     f_jk = G_jk / D_00^(j+k+1) where G_jk = Num_jk D_00^(j+k) - D_10 G_(j-1,k) - D_01 G_(j,k-1) - D_11 D_00 G_(j-1,k-1),
     all Gaussian integers. Checked against f at two rational points, exactly (the truncation error is below 10^-25 there).
  3. sum_{j,k <= 28} |f_jk| r^(j+k) > 1 from below: each |G_jk| rounded down and |D_00| up at 10^-40 (integer square roots),
     every term floored at 10^-80. All terms of the full sum are nonnegative, so ||f_r||_1 > 1, and K_2 < r (the
     dilation norm increases with r and is continuous, so it also exceeds 1 slightly below r).
Python 3.9+, standard library only."""
import random
import sys
import time
from fractions import Fraction
from math import isqrt

L0, T0, S0, N0 = 2500000000, 3067398171, 10 ** 15, 28
R0 = Fraction(302825279492, 10 ** 12)
SCALE = 10 ** 80


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gsc(c, a):
    return (c * a[0], c * a[1])


def poly_UV():
    """U and V as {(j, k): Gaussian integer}"""
    U = {(0, 0): (1, 0), (1, 0): (1, 0), (0, 1): (-1, 0), (1, 1): (-1, 0)}
    V = {(0, 0): (1, 0), (1, 1): (1, 0)}
    return U, V


def lin(a, A, b, B):
    """a A + b B for Gaussian scalars a, b and polynomials A, B"""
    out = {}
    for P_, c in ((A, a), (B, b)):
        for k, v in P_.items():
            out[k] = gadd(out.get(k, (0, 0)), gmul(c, v))
    return out


def coefficients(L, T, S, N):
    U, V = poly_UV()
    P = lin((L, 0), U, (0, T), V)
    Q = lin((L, 0), V, (0, T), U)
    Num = lin((S, 0), Q, (-1, 0), P)
    D = lin((S, 0), Q, (1, 0), P)
    assert set(D) <= {(0, 0), (1, 0), (0, 1), (1, 1)} and set(Num) <= set(D) | {(0, 0)}
    d00, d10, d01, d11 = (D.get(k, (0, 0)) for k in ((0, 0), (1, 0), (0, 1), (1, 1)))
    assert d00 == ((S + 1) * L, (S + 1) * T)
    pw = [(1, 0)]
    for _ in range(2 * N + 2):
        pw.append(gmul(pw[-1], d00))
    G = {}
    for j in range(N + 1):
        for k in range(N + 1):
            g = gmul(Num.get((j, k), (0, 0)), pw[j + k])
            if j > 0:
                g = gadd(g, gsc(-1, gmul(d10, G[j - 1, k])))
            if k > 0:
                g = gadd(g, gsc(-1, gmul(d01, G[j, k - 1])))
            if j > 0 and k > 0:
                g = gadd(g, gsc(-1, gmul(gmul(d11, d00), G[j - 1, k - 1])))
            G[j, k] = g
    return G, d00, (U, V, P, Q, Num, D)


def series_value(G, d00, N, z, w):
    """the truncated series at a rational point (z, w as pairs of Fractions), exactly"""
    D2 = d00[0] ** 2 + d00[1] ** 2
    cj = (d00[0], -d00[1])
    tot = (Fraction(0), Fraction(0))
    zp = [(Fraction(1), Fraction(0))]
    wp = [(Fraction(1), Fraction(0))]
    for _ in range(N):
        zp.append(gmul(zp[-1], z))
        wp.append(gmul(wp[-1], w))
    cpow = [(1, 0)]
    for _ in range(2 * N + 2):
        cpow.append(gmul(cpow[-1], cj))
    for (j, k), g in G.items():
        n = j + k + 1
        num = gmul(g, cpow[n])
        c = (Fraction(num[0], D2 ** n), Fraction(num[1], D2 ** n))
        tot = gadd(tot, gmul(c, gmul(zp[j], wp[k])))
    return tot


def f_exact(z, w, L, T, S):
    U = gmul(gadd((1, 0), z), gadd((1, 0), gsc(-1, w)))
    V = gadd((1, 0), gmul(z, w))
    P = gadd(gsc(L, U), gmul((0, T), V))
    Q = gadd(gsc(L, V), gmul((0, T), U))
    num = gadd(gsc(S, Q), gsc(-1, P))
    den = gadd(gsc(S, Q), P)
    d2 = den[0] ** 2 + den[1] ** 2
    q = gmul(num, (den[0], -den[1]))
    return (q[0] / d2, q[1] / d2)


def bohr_lower(G, d00, r, N):
    """rigorous lower bound, as an integer at SCALE, for sum_{j,k <= N} |f_jk| r^(j+k)"""
    E = 40
    A_up = Fraction(isqrt((d00[0] ** 2 + d00[1] ** 2) * 10 ** (2 * E)) + 1, 10 ** E)   # >= |D_00|
    tot = 0
    for (j, k), g in G.items():
        if j > N or k > N:
            continue
        mag = Fraction(isqrt((g[0] ** 2 + g[1] ** 2) * 10 ** (2 * E)), 10 ** E)   # <= |G_jk|
        term = mag * r ** (j + k) / A_up ** (j + k + 1)
        tot += term.numerator * SCALE // term.denominator           # floored
    return tot


def schur_identity_checks(L, T):
    """the two identities behind step 1, symbolically in the real coordinates (polynomials as dicts of exponents)."""
    # variables: x1, y1 (z = x1 + i y1), x2, y2 (w = x2 + i y2); a polynomial is {(e1, e2, e3, e4): Fraction}
    def mul(p, q):
        out = {}
        for a, u in p.items():
            for b, v in q.items():
                k = tuple(i + j for i, j in zip(a, b))
                out[k] = out.get(k, 0) + u * v
        return {k: v for k, v in out.items() if v}

    def add(*ps):
        out = {}
        for p in ps:
            for k, v in p.items():
                out[k] = out.get(k, 0) + v
        return {k: v for k, v in out.items() if v}

    def c(v):
        return {(0, 0, 0, 0): Fraction(v)}

    x1, y1, x2, y2 = ({tuple(int(i == j) for i in range(4)): Fraction(1)} for j in range(4))
    neg = c(-1)
    # complex numbers as (re, im) pairs of polynomials
    def cm(a, b):
        return (add(mul(a[0], b[0]), mul(neg, mul(a[1], b[1]))), add(mul(a[0], b[1]), mul(a[1], b[0])))

    def ca(a, b):
        return (add(a[0], b[0]), add(a[1], b[1]))

    def cs(s, a):
        return (mul(c(s), a[0]), mul(c(s), a[1]))

    def conj(a):
        return (a[0], mul(neg, a[1]))

    one = (c(1), {})
    z, w = (x1, y1), (x2, y2)
    U = cm(ca(one, z), ca(one, cs(-1, w)))
    V = ca(one, cm(z, w))
    UV = cm(U, conj(V))
    a2 = add(mul(x1, x1), mul(y1, y1))
    b2 = add(mul(x2, x2), mul(y2, y2))
    target = add(c(1), mul(neg, mul(a2, b2)), mul(x1, add(c(1), mul(neg, b2))), mul(neg, mul(x2, add(c(1), mul(neg, a2)))))
    ok1 = add(UV[0], mul(neg, target)) == {}
    P = ca(cs(L, U), cm((c(0), c(T)), V))
    Q = ca(cs(L, V), cm((c(0), c(T)), U))
    QP = cm(Q, conj(P))
    ok2 = add(QP[0], mul(c(-(L * L + T * T)), UV[0])) == {}
    # (1 - a)(1 - b)(1 - ab) = 1 - a - b + a^2 b + a b^2 - a^2 b^2, the value of the bound at Re z = -a, Re w = b
    ok3 = True   # both sides have degree at most 2 in a and in b, so agreement on a 3 x 3 grid is an identity
    for a in (Fraction(0), Fraction(1, 3), Fraction(2, 3)):
        for b in (Fraction(0), Fraction(1, 3), Fraction(2, 3)):
            ok3 &= (1 - a) * (1 - b) * (1 - a * b) == 1 - a * a * b * b - a * (1 - b * b) - b * (1 - a * a)
    return ok1, ok2, ok3


def main():
    t0 = time.time()
    controls = "--controls" in sys.argv
    ok = True
    print("59a row [P2026]: K_2 < 0.302825279492 (U = (1+z)(1-w), V = 1+zw, P = LU + iTV, Q = LV + iTU, "
          "f = (SQ - P)/(SQ + P); L = 2500000000, T = 3067398171, S = 10^15, N = 28)")
    ok1, ok2, ok3 = schur_identity_checks(L0, T0)
    print(f"  Re(U conj V) = 1 - |z|^2|w|^2 + Re z (1 - |w|^2) - Re w (1 - |z|^2), as polynomials: {ok1}")
    print(f"  Re(Q conj P) = (L^2 + T^2) Re(U conj V), as polynomials: {ok2}")
    print(f"  1 - a^2 b^2 - a(1 - b^2) - b(1 - a^2) = (1 - a)(1 - b)(1 - ab) (degree <= 2 in each, on a 3 x 3 grid): {ok3}")
    rnd = random.Random(1)
    good = True
    for _ in range(2000):   # exact, at random rational points of the bidisc: Re(U conj V) > 0 and |SQ - P| < |SQ + P|
        pts = []
        for _ in range(2):
            while True:
                x, y = Fraction(rnd.randint(-999, 999), 1000), Fraction(rnd.randint(-999, 999), 1000)
                if x * x + y * y < 1:
                    pts.append((x, y))
                    break
        (x1, y1), (x2, y2) = pts
        U = ((1 + x1) * (1 - x2) + y1 * y2, y1 * (1 - x2) - (1 + x1) * y2)      # (1 + z)(1 - w)
        V = (1 + x1 * x2 - y1 * y2, x1 * y2 + y1 * x2)                          # 1 + zw
        re_uv = U[0] * V[0] + U[1] * V[1]
        P = (L0 * U[0] - T0 * V[1], L0 * U[1] + T0 * V[0])
        Q = (L0 * V[0] - T0 * U[1], L0 * V[1] + T0 * U[0])
        num = (S0 * Q[0] - P[0], S0 * Q[1] - P[1])
        den = (S0 * Q[0] + P[0], S0 * Q[1] + P[1])
        good &= re_uv > 0 and num[0] ** 2 + num[1] ** 2 < den[0] ** 2 + den[1] ** 2
    print(f"  at 2,000 random rational points of the bidisc, exactly: Re(U conj V) > 0 and |SQ - P| < |SQ + P|: {good}")
    ok &= ok1 and ok2 and ok3 and good
    G, d00, _ = coefficients(L0, T0, S0, N0)
    print(f"  {len(G)} coefficients f_jk = G_jk / D_00^(j+k+1), D_00 = (S + 1)(L + iT), exact Gaussian integers "
          f"(largest |G| has {max(len(str(abs(g[0]))) for g in G.values())} digits)")
    good = True
    for z, w in (((Fraction(1, 20), Fraction(1, 50)), (Fraction(-3, 100), Fraction(1, 25))),
                 ((Fraction(0), Fraction(-3, 50)), (Fraction(7, 100), Fraction(0)))):
        a, b = series_value(G, d00, N0, z, w), f_exact(z, w, L0, T0, S0)
        good &= abs(a[0] - b[0]) < Fraction(1, 10 ** 25) and abs(a[1] - b[1]) < Fraction(1, 10 ** 25)
    print(f"  the truncated series agrees with f at two rational points, exactly, to 10^-25 (the tail): {good}")
    ok &= good
    lo = bohr_lower(G, d00, R0, N0)
    margin = Fraction(lo - SCALE, SCALE)
    good = lo > SCALE
    print(f"  sum_(j,k <= 28) |f_jk| r^(j+k) at r = 0.302825279492 >= 1 + {float(margin):.6e} (rigorous): {good}")
    good2 = margin > Fraction(1, 10 ** 26)
    print(f"  the row's 'greater than 1 + 10^-26': {good2}")
    ok &= good and good2
    if controls:
        print("Red controls (each must leave the finite sum below 1):")
        for what, (L, T, S, N, r) in (("r = 0.3028252794 (10^-10 lower)", (L0, T0, S0, N0, Fraction(3028252794, 10 ** 10))),
                                      ("r = 0.3006, Knese's lower bound for K_2", (L0, T0, S0, N0, Fraction(3006, 10 ** 4))),
                                      ("T = 0 (no complex phase)", (L0, 0, S0, N0, R0)),
                                      ("L = 0 (no real part)", (0, T0, S0, N0, R0))):
            Gc, dc, _ = coefficients(L, T, S, N)
            v, vhi = bohr_lower(Gc, dc, r, N), bohr_upper(Gc, dc, r, N)
            refused = vhi < SCALE
            print(f"  {what}: the finite sum is 1 - {float(Fraction(SCALE - vhi, SCALE)):.4e} at most: "
                  f"{'below 1, as expected' if refused else 'NOT BELOW 1 (control failed)'}")
            ok &= refused
        Gc, dc, _ = coefficients(L0, T0, S0, N0)
        r1 = R0 - Fraction(1, 10 ** 12)
        v = bohr_lower(Gc, dc, r1, N0)
        print(f"  (not a control: at r = 0.302825279491, one unit lower, the same sum is still >= 1 + "
              f"{float(Fraction(v - SCALE, SCALE)):.3e})")
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


def bohr_upper(G, d00, r, N):
    E = 40
    A_lo = Fraction(isqrt((d00[0] ** 2 + d00[1] ** 2) * 10 ** (2 * E)), 10 ** E)       # <= |D_00|
    tot = 0
    for (j, k), g in G.items():
        if j > N or k > N:
            continue
        mag = Fraction(isqrt((g[0] ** 2 + g[1] ** 2) * 10 ** (2 * E)) + 1, 10 ** E)   # >= |G_jk|
        term = mag * r ** (j + k) / A_lo ** (j + k + 1)
        tot += -((-term.numerator * SCALE) // term.denominator)     # ceiled
    return tot


if __name__ == "__main__":
    sys.exit(main())
