#!/usr/bin/env python3
"""Exact replay of the two linear-programming rows of the 29a page (kissing number in dimension 5, teorth/optimization-
problems, main 87b4ebb): tau_5 <= 46.345 [OS1979] and tau_5 <= 48 [Cox1963] (48 is L_5(5, 1/2), the bound the cited
survey [BDM2012] states). We could not read the 1979 paper or Coxeter's; both certificates below are ours.

The bound (Delsarte; Kabatiansky-Levenshtein, Odlyzko-Sloane), re-derived here. Let G_k be the Gegenbauer polynomials
of S^4, normalised by G_k(1) = 1: G_0 = 1, G_1 = t, (k + 3) G_{k+1} = (2k + 3) t G_k - k G_{k-1}. For any finite
X in S^4, sum_{x,y in X} G_k(<x,y>) >= 0 (Schoenberg; the addition formula), and in a kissing configuration, read as
unit vectors, distinct points have <x,y> <= 1/2. If f = sum_k f_k G_k with f_0 > 0, f_k >= 0 (k >= 1) and
f(t) <= 0 on [-1, 1/2], then
    f_0 N^2 <= sum_{x,y} f(<x,y>) = N f(1) + sum_{x != y} f(<x,y>) <= N f(1),  so  tau_5 <= f(1)/f_0.

Checked here, all in exact rational arithmetic:
  * G_k by the recurrence and, as a second route, by the explicit sum for C_k^(3/2)(t)/C_k^(3/2)(1), k <= 40.
  * [OS1979]: our degree-10 polynomial F = (1 - eps) + sum f_k G_k (f_6 = f_7 = f_8 = 0; the f_k are rationals from our
    own floating-point LP search, lp_search.py, and eps = 10^-6): every f_k >= 0; the constant and every other
    coefficient also by orthogonal projection (weight (1 - t^2) on [-1, 1], exact integrals); F < 0 on [-1, 1/2] by a
    Sturm sequence over Q (no root, F(-1) < 0) and, as a second route, by Bernstein coefficients on a subdivision of
    [-1, 1/2] (all negative); the bound F(1)/F_0 as an exact fraction, compared with 46.345.
  * [Cox1963] value 48: Levenshtein's degree-5 polynomial L = (t - 1/2)(t + 5/7)^2 (t + 1/7)^2 (our LP search finds
    its double roots near -0.714 and -0.143): its Gegenbauer coefficients exactly (projection, and a triangular solve
    as the second route), all positive; L <= 0 on [-1, 1/2] from the factorisation; L(1)/L_0 = 48 exactly.
  * Dual certificates (what no polynomial can do): points t_i in [-1, 1/2] and weights y_i >= 0 with
    sum_i y_i G_k(t_i) >= -1 for k = 1..K give f(1)/f_0 >= 1 + sum_i y_i for every admissible f of degree <= K
    (0 >= sum_i y_i f(t_i)/f_0 >= sum_i y_i - (f(1)/f_0 - 1)). With K = 9 this shows degree 10 is needed to reach
    46.345; with K = 40, that the LP bound for degree <= 40 is within 5e-5 of our certificate.
Red controls (--controls): the same F with s = 0.51 in place of 1/2, with f_10 = 0, and with the Gegenbauer
polynomials of S^3 in place of those of S^4 are refused; a claim of 46.3375 from F is refused; the degree-9 dual refuses
46.345 for every polynomial of degree <= 9; L with 5/7 replaced by 2/3 is refused as a certificate of 48; L is refused
as a certificate of 47.9.
Usage: python3 -I os1979_check.py [--controls]   (standard library only; a few seconds)
"""
import sys
from fractions import Fraction as Fr
from math import comb, factorial

CONTROLS = "--controls" in sys.argv
KMAX = 40

# ---------------------------------------------------------------- the certificates (ours)
# [OS1979] row: degree-10 Gegenbauer coefficients f_1..f_10 (from lp_search.py 10 6001, limit_denominator(10^6))
F_COEF = ["3414134/751195", "6599613/699224", "11310459/864874", "8681257/812868", "1684978/285159",
          "0", "0", "0", "1313483/855749", "19699/128763"]
EPS = Fr(1, 10 ** 6)
PRINTED_OS = Fr(46345, 1000)
# Levenshtein's polynomial for L_5(5, 1/2): roots 1/2, -5/7 (double), -1/7 (double)
LEV_ROOTS = [(Fr(1, 2), 1), (Fr(-5, 7), 2), (Fr(-1, 7), 2)]
# dual certificates: points t (grid points of the LP search) and weights y from the LP's dual solutions
# (lp_search.py 9 6001 and 10 6001), the weights scaled by (1 - 10^-9) after rounding
DUAL9 = [("-0.78025", "0.7333583251195"), ("-0.78", "3.849063142348"), ("-0.4615", "7.195802995762"),
         ("-0.092", "11.38498102542"), ("-0.09175", "7.010732871727"), ("0.5", "15.17197632414")]
DUAL10 = [("-0.7735", "4.507040707515"), ("-0.77325", "0.3835041966617"), ("-0.43975", "3.978342164210"),
          ("-0.4395", "3.316303960857"), ("-0.0885", "17.72306900275"), ("-0.08825", "0.2659282103352"),
          ("0.5", "15.16338387229")]


# ---------------------------------------------------------------- polynomials over Q (coefficient lists, low first)
def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(p, q):
    n = max(len(p), len(q))
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def pscale(p, c):
    return trim([c * x for x in p])


def pmul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return trim(out)


def peval(p, t):
    v = Fr(0)
    for c in reversed(p):
        v = v * t + c
    return v


def pderiv(p):
    return trim([i * p[i] for i in range(1, len(p))]) if len(p) > 1 else [Fr(0)]


def prem(a, b):
    a = list(a)
    while len(a) >= len(b) and any(a):
        c = a[-1] / b[-1]
        s = len(a) - len(b)
        for i, x in enumerate(b):
            a[s + i] -= c * x
        a = trim(a[:-1]) if a[-1] == 0 else trim(a)
        if len(a) < len(b):
            break
    return trim(a)


def sturm_roots(p, a, b):
    """number of distinct real roots of p in (a, b] (p squarefree part handled by the classical sequence)"""
    seq = [trim(p), pderiv(p)]
    while len(seq[-1]) > 1 or seq[-1][0] != 0:
        r = prem(seq[-2], seq[-1])
        if len(r) == 1 and r[0] == 0:
            break
        seq.append(pscale(r, -1))
        if len(seq[-1]) == 1:
            break

    def changes(x):
        vals = [peval(q, x) for q in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for u, v in zip(vals, vals[1:]) if (u < 0) != (v < 0))
    return changes(a) - changes(b)


def integral_w(p):
    """exact integral of p(t)(1 - t^2) over [-1, 1]"""
    q = pmul(p, [Fr(1), Fr(0), Fr(-1)])
    return sum(c * (Fr(2, i + 1) if i % 2 == 0 else 0) for i, c in enumerate(q))


def bernstein_negative(p, a, b, depth=0, maxdepth=30, stats=None):
    """True if p < 0 on [a, b], decided by Bernstein coefficients on a bisection of [a, b] (exact)."""
    # p(a + (b - a) u) in the monomial basis of u
    n = len(p) - 1
    q = [Fr(0)] * (n + 1)
    h = b - a
    for i, c in enumerate(p):
        # c (a + h u)^i
        for j in range(i + 1):
            q[j] += c * comb(i, j) * a ** (i - j) * h ** j
    bern = [sum(Fr(comb(k, j), comb(n, j)) * q[j] for j in range(k + 1)) for k in range(n + 1)]
    if stats is not None:
        stats[0] += 1
    if max(bern) < 0:
        return True
    if bern[0] >= 0 or bern[-1] >= 0 or depth >= maxdepth:
        return False
    m = (a + b) / 2
    return (bernstein_negative(p, a, m, depth + 1, maxdepth, stats)
            and bernstein_negative(p, m, b, depth + 1, maxdepth, stats))


# ---------------------------------------------------------------- Gegenbauer polynomials of S^{n-1}
def gegen_rec(K, n=5):
    G = [[Fr(1)], [Fr(0), Fr(1)]]
    for k in range(1, K):
        G.append(trim(padd(pscale([Fr(0)] + G[k], Fr(2 * k + n - 2, k + n - 2)),
                           pscale(G[k - 1], Fr(-k, k + n - 2)))))
    return G[:K + 1]


def gegen_explicit(K):
    """C_k^(3/2)(t) / C_k^(3/2)(1) from C_k^(lam)(t) = sum_j (-1)^j (lam)_{k-j} / (j! (k-2j)!) (2t)^{k-2j}"""
    def poch(x, m):
        out = Fr(1)
        for i in range(m):
            out *= x + i
        return out
    lam = Fr(3, 2)
    G = []
    for k in range(K + 1):
        p = [Fr(0)] * (k + 1)
        for j in range(k // 2 + 1):
            p[k - 2 * j] += (-1) ** j * poch(lam, k - j) / (factorial(j) * factorial(k - 2 * j)) * 2 ** (k - 2 * j)
        norm = Fr((k + 1) * (k + 2), 2)
        G.append(trim([c / norm for c in p]))
    return G


def combine(coefs, G):
    p = [Fr(0)]
    for k, c in enumerate(coefs):
        if c:
            p = padd(p, pscale(G[k], c))
    return p


def project(p, G):
    """Gegenbauer coefficients of p by orthogonality under the weight (1 - t^2)"""
    return [integral_w(pmul(p, G[k])) / integral_w(pmul(G[k], G[k])) for k in range(len(p))]


def triangular(p, G):
    """Gegenbauer coefficients of p by peeling off the leading terms (second route)"""
    p = list(p)
    out = [Fr(0)] * len(p)
    for k in range(len(p) - 1, -1, -1):
        c = (p[k] if k < len(p) else 0) / G[k][k]
        out[k] = c
        p = padd(p, pscale(G[k], -c)) + [Fr(0)] * len(p)
        p = p[:len(out)]
    return out


def check_certificate(name, coefs, G, s, claim, out, verbose=True):
    """coefs[0] is f_0; returns (ok, bound)"""
    ok = True
    f0 = coefs[0]
    if f0 <= 0 or any(c < 0 for c in coefs[1:]):
        out.append("  %s: a coefficient has the wrong sign" % name)
        return False, None
    P = combine(coefs, G)
    bound = peval(P, Fr(1)) / f0
    nroots = sturm_roots(P, Fr(-1), s)
    left = peval(P, Fr(-1))
    right = peval(P, s)
    neg_sturm = nroots == 0 and left < 0 and right < 0
    stats = [0]
    neg_bern = bernstein_negative(P, Fr(-1), s, stats=stats)
    if verbose:
        out.append("  %s: F(-1) = %.6g, F(s) = %.6g, Sturm: %d distinct roots in (-1, s]; Bernstein: %s on %d pieces"
                   % (name, float(left), float(right), nroots, "all negative" if neg_bern else "not all negative",
                      stats[0]))
    ok = neg_sturm and neg_bern
    if not ok:
        out.append("  %s: F is not negative on [-1, %s]" % (name, s))
    if claim is not None and not bound <= claim:
        out.append("  %s: the bound %.10f is above the claim %s" % (name, float(bound), float(claim)))
        ok = False
    return ok, bound


def dual_lower(pairs, K, G, scale=Fr(1) - Fr(1, 10 ** 9)):
    """returns (feasible for k = 1..K, lower bound 1 + sum y) for the dual points"""
    pts = [(Fr(t), Fr(y) * scale) for t, y in pairs]
    feas = True
    worst = None
    for k in range(1, K + 1):
        v = sum(y * peval(G[k], t) for t, y in pts)
        if v < -1:
            feas = False
        worst = v if worst is None or v < worst else worst
    assert all(-1 <= t <= Fr(1, 2) and y >= 0 for t, y in pts)
    return feas, 1 + sum(y for _, y in pts), worst


def main():
    out = []
    ok = True
    G = gegen_rec(KMAX)
    G2 = gegen_explicit(KMAX)
    same = all(trim(a) == trim(b) for a, b in zip(G, G2)) and all(peval(g, Fr(1)) == 1 for g in G)
    out.append("29a, kissing number in dimension 5: the LP rows [OS1979] 46.345 and [Cox1963] 48")
    out.append("Gegenbauer polynomials of S^4, k <= %d: recurrence = explicit C_k^(3/2) sum, all G_k(1) = 1: %s"
               % (KMAX, same))
    ok &= same

    # [OS1979]
    fk = [Fr(1) - EPS] + [Fr(x) for x in F_COEF]
    out.append("[OS1979] tau_5 <= 46.345: our degree-10 polynomial F = (1 - 10^-6) + sum f_k G_k, f_k = "
               + ", ".join(F_COEF))
    P = combine(fk, G)
    proj = project(P, G)
    same_coef = proj[:len(fk)] == fk and all(c == 0 for c in proj[len(fk):])
    out.append("  coefficients recovered by orthogonal projection (weight 1 - t^2, exact integrals): %s" % same_coef)
    ok &= same_coef
    good, bound = check_certificate("F", fk, G, Fr(1, 2), PRINTED_OS, out)
    ok &= good
    out.append("  bound F(1)/F_0 = %s / %s" % (bound.numerator, bound.denominator))
    out.append("        = %.12f <= 46.345: %s (margin %.6f)" % (float(bound), bound <= PRINTED_OS,
                                                                   float(PRINTED_OS - bound)))

    # duals
    feas9, lb9, w9 = dual_lower(DUAL9, 9, G)
    out.append("  dual, degree <= 9: sum_i y_i G_k(t_i) >= -1 for k = 1..9: %s (least %.12f); every admissible "
               "polynomial of degree <= 9 gives at least %.10f > 46.345: %s" % (feas9, float(w9), float(lb9),
                                                                            lb9 > PRINTED_OS))
    ok &= feas9 and lb9 > PRINTED_OS
    feas40, lb40, w40 = dual_lower(DUAL10, KMAX, G)
    out.append("  dual, degree <= %d: feasible for k = 1..%d: %s (least %.12f); every admissible polynomial of degree "
               "<= %d gives at least %.10f" % (KMAX, KMAX, feas40, float(w40), KMAX, float(lb40)))
    ok &= feas40 and lb40 <= bound
    out.append("  so the LP bound with polynomials of degree <= %d lies in [%.9f, %.9f]; the printed 46.345 holds, "
               "and degree <= 9 cannot give it" % (KMAX, float(lb40), float(bound)))

    # Levenshtein / [Cox1963] 48
    L = [Fr(1)]
    for r, m in LEV_ROOTS:
        for _ in range(m):
            L = pmul(L, [-r, Fr(1)])
    lc = project(L, G)
    lc2 = triangular(L, G)
    out.append("[Cox1963] tau_5 <= 48 (L_5(5, 1/2) = 48 as [BDM2012] states): L = (t - 1/2)(t + 5/7)^2 (t + 1/7)^2")
    out.append("  Gegenbauer coefficients (projection = triangular solve: %s): %s" % (lc == lc2[:len(lc)],
                                                                                 ", ".join(str(c) for c in lc)))
    ok &= lc == lc2[:len(lc)]
    pos = lc[0] > 0 and all(c >= 0 for c in lc[1:])
    ratio = peval(L, Fr(1)) / lc[0]
    out.append("  all coefficients >= 0: %s; L <= 0 on [-1, 1/2] (t - 1/2 <= 0 times squares); L(1)/L_0 = %s"
               % (pos, ratio))
    ok &= pos and ratio == 48

    if CONTROLS:
        out.append("")
        out.append("red controls (each must be refused):")
        res = []
        sub = []
        g, _ = check_certificate("F with s = 0.51", fk, G, Fr(51, 100), None, sub, verbose=False)
        res.append(("F on [-1, 0.51] (s = 0.51)", not g))
        f_no10 = fk[:10] + [Fr(0)]
        g, _ = check_certificate("F with f_10 = 0", f_no10, G, Fr(1, 2), None, sub, verbose=False)
        res.append(("F with f_10 = 0", not g))
        G4 = gegen_rec(KMAX, n=4)
        g, b4 = check_certificate("F with S^3 polynomials", fk, G4, Fr(1, 2), PRINTED_OS, sub, verbose=False)
        res.append(("the same coefficients on the Gegenbauer polynomials of S^3", not g))
        res.append(("a claim of 46.3375 from F", not bound <= Fr(463375, 10000)))
        res.append(("46.345 from any polynomial of degree <= 9 (the degree-9 dual)", lb9 > PRINTED_OS))
        L2 = pmul(pmul(pmul(pmul([-Fr(1, 2), Fr(1)], [Fr(2, 3), Fr(1)]), [Fr(2, 3), Fr(1)]), [Fr(1, 7), Fr(1)]),
                  [Fr(1, 7), Fr(1)])
        c2 = project(L2, G)
        ok2 = c2[0] > 0 and all(c >= 0 for c in c2[1:]) and peval(L2, Fr(1)) / c2[0] <= 48
        res.append(("L with 5/7 replaced by 2/3 as a certificate of 48 (ratio %.6f, least coefficient %.4g)"
                    % (float(peval(L2, Fr(1)) / c2[0]), float(min(c2[1:]))), not ok2))
        res.append(("L as a certificate of 47.9", not ratio <= Fr(479, 10)))
        for what, refused in res:
            out.append("  %s: %s" % (what, "refused, as expected" if refused else "NOT REFUSED"))
        allc = all(r for _, r in res)
        out.append("controls all as expected: %s" % allc)
        ok &= allc
    out.append("RESULT: %s" % ("PASS (tau_5 <= %.10f < 46.345, and tau_5 <= 48 from L)" % float(bound) if ok
                               else "FAIL"))
    print("\n".join(out))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
