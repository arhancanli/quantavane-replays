#!/usr/bin/env python3
"""Replay of two 3a lower-bound rows (teorth/optimizationproblems, main 87b4ebb) through the lemma the page states,
C_3a >= 1 + log(|U-U|/|U+U|) / log(2 max U + 1) for a finite set U of non-negative integers containing 0
(Gyarmati-Hennecart-Ruzsa 2007), and its limit form:

  [Ru96]   U = {0, 1, 3} as printed on the page: |U+U| = 6, |U-U| = 7, q = 7, so C_3a >= 2 - log 6/log 7 = 1.0792...;
           both sets are counted exactly by two routes (Python sets; sorted lists of all pair values).
  [Z2025]  Zheng, arXiv 2506.01896v1, Theorem 1: for the sets U(m, L, B) of Gerbicz with B fixed and m, L -> infinity,
           theta >= T(B, r, a) / log(2B + 1) for every r > 0 and a in (0, min(1, 1/r)), where
             T = log 2 + ar log B + (1 - ar) log(B + 1) - I(ar, 1) - ar I((1-a)/a, B-1) - (1 - ar) I(r/(1 - ar), B)
                 + I(2r, 2B),
             I(c, B) = sup_t ( tc - log((1 + e^t + ... + e^{Bt})/(B + 1)) ) for c < B/2, and 0 for c >= B/2.
           The row prints 1.173077* ("not certified by a finite-depth computation").

What is computed for [Z2025] (the paper's LaTeX read as text; its MATLAB code was not run):
  1. the restricted supremum the paper tabulates, for B = 3..7, by our own floating-point search (Newton for each I,
     golden section in a, then in r): it must reproduce the paper's table (epsilon = 1e-10 column) to 1e-13;
  2. the certificate: B = 5, r0 = 8050530/10^7, a0 = 4668230/10^7 (rationals near the optimum), every term enclosed
     in 60-digit decimal interval arithmetic. For each I with a minus sign we need an upper bound: phi(t) = tc - log M(t)
     is concave, so with t1 < t2 around the maximiser (phi'(t1) > 0 > phi'(t2), checked on enclosures),
     sup phi <= max(phi(t1), phi(t2)) + max(phi'(t1), -phi'(t2)) (t2 - t1). For I(2r, 2B), with a plus sign, the
     value phi(t) at one t is a lower bound. Logarithms and exponentials are correctly rounded at 60 digits and moved
     outward by 10^-50; sums and products are rounded outward.
  3. the conditions: a0 in (0, min(1, 1/r0)) and each I-argument below its B/2 (so the sup formula applies).
What rests on mathematics (read, not checked): the GHR lemma, Gerbicz's formulas for |U+U| and |U-U|, and Zheng's
Lemmas 1-2 (Cramer's theorem and semicontinuity), which turn the expression into a bound on theta.

Red controls (--controls): the same certificate at B = 4 (best value 1.17214 < 1.173077) and at r0 = 0.6 fail; the
term I(2r, 2B) dropped fails; a claimed 1.1730773 (above the float supremum 1.17307727978...) is not supported; for
[Ru96], U = {0, 1, 2} gives 1, and U = {0, 1, 4} gives 1.0702 < 1.0792.
usage: python3 -I z2025_check.py CERTS_ROOT [--controls]     (CERTS_ROOT holds z2025/sumdif.tex)
"""
import math
import os
import re
import sys
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction

getcontext().prec = 60
EPS = Decimal(10) ** -50
ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ---------------------------------------------------------------- [Ru96]
def counts(U):
    U = sorted(set(U))
    s = {x + y for x in U for y in U}
    d = {x - y for x in U for y in U}
    pairs_s = sorted(x + y for i, x in enumerate(U) for y in U[i:])
    pairs_d = sorted(x - y for x in U for y in U)
    nd = lambda xs: sum(1 for i, x in enumerate(xs) if i == 0 or x != xs[i - 1])
    need(len(s) == nd(pairs_s) and len(d) == nd(pairs_d), "the two counting routes disagree")
    return len(s), len(d), 2 * max(U) + 1


def ghr(U):
    s, d, q = counts(U)
    return s, d, q, 1 + (Decimal(d) / s).ln() / Decimal(q).ln()


# ---------------------------------------------------------------- interval helpers
def add(x, y):
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        lo = x[0] + y[0]
        c.rounding = ROUND_CEILING
        hi = x[1] + y[1]
    return (lo, hi)


def neg(x):
    return (-x[1], -x[0])


def mul(x, y):
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        lo = min(x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
        c.rounding = ROUND_CEILING
        hi = max(x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
    return (lo, hi)


def div(x, y):
    need(y[0] > 0, "division by an interval that is not positive")
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        lo = min(x[0] / y[0], x[0] / y[1])
        c.rounding = ROUND_CEILING
        hi = max(x[1] / y[0], x[1] / y[1])
    return (lo, hi)


def pt(v):
    return (Decimal(v), Decimal(v))


def ln_i(x):
    need(x[0] > 0, "log of a non-positive interval")
    return (x[0].ln() - EPS, x[1].ln() + EPS)


def exp_i(x):
    return (x[0].exp() * (1 - EPS), x[1].exp() * (1 + EPS))


def frac(f):
    return div(pt(f.numerator), pt(f.denominator))


def phi_and_slope(t, c, B):
    """Enclosures of phi(t) = t c - log M(t) and phi'(t) = c - mean_t, M(t) = (1/(B+1)) sum_k e^{kt}."""
    T = pt(t)
    S = pt(0)
    S1 = pt(0)
    for k in range(B + 1):
        e = exp_i(mul(pt(k), T))
        S = add(S, e)
        S1 = add(S1, mul(pt(k), e))
    logM = add(ln_i(S), neg(ln_i(pt(B + 1))))
    phi = add(mul(T, c), neg(logM))
    slope = add(c, neg(div(S1, S)))
    return phi, slope


def t_star(c, B):
    """Newton for mean_t = c in plain Decimal (only a good point is needed; the bounds come from enclosures)."""
    t = Decimal(0)
    cc = Decimal(c.numerator) / Decimal(c.denominator)
    for _ in range(200):
        w = [(k * t).exp() for k in range(B + 1)]
        s = sum(w)
        mu = sum(k * x for k, x in enumerate(w)) / s
        var = sum((k - mu) ** 2 * x for k, x in enumerate(w)) / s
        step = (mu - cc) / var
        t -= step
        if abs(step) < Decimal(10) ** -45:
            break
    return t


def I_upper(c, B):
    if c >= Fraction(B, 2):
        return pt(0)
    cI = frac(c)
    ts = t_star(c, B)
    d = Decimal(10) ** -25
    p1, s1 = phi_and_slope(ts - d, cI, B)
    p2, s2 = phi_and_slope(ts + d, cI, B)
    need(s1[0] > 0 and s2[1] < 0, "the slopes do not bracket the maximiser")
    hi = max(p1[1], p2[1]) + max(s1[1], -s2[0]) * (2 * d)
    return (hi - EPS, hi + EPS)          # an upper bound is hi (+ EPS); the interval is used only as an upper bound


def I_lower(c, B):
    if c >= Fraction(B, 2):
        return pt(0)
    p, _ = phi_and_slope(t_star(c, B), frac(c), B)
    return p


def T_lower(B, r, a, drop_last=False):
    """A rigorous lower bound for Zheng's expression divided by log(2B + 1)."""
    ar = a * r
    need(0 < a < min(1, 1 / r), "a is outside (0, min(1, 1/r))")
    need(ar < 1, "ar >= 1")
    terms = [ln_i(pt(2)), mul(frac(ar), ln_i(pt(B))), mul(frac(1 - ar), ln_i(pt(B + 1)))]
    up = [I_upper(ar, 1), mul(frac(ar), I_upper((1 - a) / a, B - 1)), mul(frac(1 - ar), I_upper(r / (1 - ar), B))]
    total = pt(0)
    for x in terms:
        total = add(total, x)
    for x in up:
        total = add(total, (-x[1], -x[1]))     # minus an upper bound
    if not drop_last:
        lo = I_lower(2 * r, 2 * B)[0]
        total = add(total, (lo, lo))
    val = div(total, ln_i(pt(2 * B + 1)))
    return val[0]


# ---------------------------------------------------------------- floating-point search (location only)
def f_stats(t, B):
    m = max(0.0, B * t)
    w = [math.exp(k * t - m) for k in range(B + 1)]
    s = sum(w)
    mu = sum(k * x for k, x in enumerate(w)) / s
    var = sum((k - mu) ** 2 * x for k, x in enumerate(w)) / s
    return m + math.log(s / (B + 1)), mu, var


def f_I(c, B):
    if c >= B / 2:
        return 0.0
    t = 0.0
    for _ in range(100):
        lm, mu, var = f_stats(t, B)
        step = (mu - c) / var
        t -= step
        if abs(step) < 1e-15:
            break
    lm, mu, var = f_stats(t, B)
    return t * c - lm


def f_theta(B, r, a):
    ar = a * r
    num = (math.log(2) + ar * math.log(B) + (1 - ar) * math.log(B + 1) - f_I(ar, 1) - ar * f_I((1 - a) / a, B - 1)
           - (1 - ar) * f_I(r / (1 - ar), B) + f_I(2 * r, 2 * B))
    return num / math.log(2 * B + 1)


def golden(f, lo, hi, it=90):
    for _ in range(it):
        m1 = lo + (hi - lo) * 0.381966011250105
        m2 = lo + (hi - lo) * 0.618033988749895
        if f(m1) < f(m2):
            lo = m1
        else:
            hi = m2
    return (lo + hi) / 2


def f_sup(B):
    best_a = lambda r: golden(lambda a: f_theta(B, r, a), 1e-9, min(1, 1 / r) - 1e-9)
    r = golden(lambda r: f_theta(B, r, best_a(r)), 0.5, 2.0)
    a = best_a(r)
    return r, a, f_theta(B, r, a)


def paper_table(path):
    text = open(path, encoding="utf-8").read()
    rows = re.findall(r"^\s*(\d+)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*(?:\{\\color\{blue\})?([0-9.]+)\}?\s*&\s*"
                      r"(?:\{\\color\{blue\})?([0-9.]+)\}?\s*\\\\", text, re.M)
    need(len(rows) == 8, f"the paper's table was not found ({len(rows)} rows)")
    need("we have $\\theta \\ge 1.173077$" in text, "the theorem statement was not found")
    return {int(r[0]): r[4] for r in rows}


def main():
    ok = True
    s, d, q, val = ghr([0, 1, 3])
    good = ((s, d, q) == (6, 7, 7) and abs(val - (2 - Decimal(6).ln() / Decimal(7).ln())) < EPS
            and str(val).startswith("1.0792"))
    print(f"[Ru96] U = {{0, 1, 3}}: |U+U| = {s}, |U-U| = {d}, 2 max U + 1 = {q} (the page: 6, 7); "
          f"1 + log(7/6)/log 7 = 2 - log 6/log 7 = {str(val)[:24]}...")
    print(f"RESULT [Ru96]: C_3a >= 2 - log 6/log 7 = 1.0792...: {'HOLDS' if good else 'FAILS'}")
    ok &= good

    table = paper_table(os.path.join(ROOT, "z2025", "sumdif.tex"))
    print("[Z2025] the restricted supremum of theta - 1, our search against the paper's table (eps = 1e-10 column):")
    tab_ok = True
    for B in range(3, 8):
        r, a, th = f_sup(B)
        dev = abs((th - 1) - float(table[B]))
        tab_ok &= dev < 1e-13
        print(f"  B = {B}: r = {r:.10f}, a = {a:.10f}, theta - 1 = {th - 1:.15f}   (paper: {table[B]}; "
              f"difference {dev:.1e})")
    r0, a0, B = Fraction(8050530, 10 ** 7), Fraction(4668230, 10 ** 7), 5
    cert = T_lower(B, r0, a0)
    conds = (r0 * a0 < Fraction(1, 2) and (1 - a0) / a0 < Fraction(B - 1, 2) and r0 / (1 - a0 * r0) < Fraction(B, 2)
             and 2 * r0 < B)
    good = tab_ok and conds and cert > Decimal("1.173077")
    print(f"  certificate B = 5, r0 = {r0}, a0 = {a0}: a0 r0 = {float(a0 * r0):.6f}, (1 - a0)/a0 = "
          f"{float((1 - a0) / a0):.6f}, r0/(1 - a0 r0) = {float(r0 / (1 - a0 * r0)):.6f}, 2 r0 = {float(2 * r0):.6f}, "
          f"each below its B/2: {conds}")
    print(f"  theta >= {str(cert)[:26]}... (60-digit interval arithmetic; margin over 1.173077: "
          f"{float(cert - Decimal('1.173077')):.3e})")
    print(f"RESULT [Z2025]: theta >= T(5, r0, a0)/log 11 > 1.173077: {'HOLDS' if good else 'FAILS'}")
    ok &= good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        c1 = T_lower(4, Fraction(6722244, 10 ** 7), Fraction(5314125, 10 ** 7))
        res.append(("B = 4 at its own optimum", c1 < Decimal("1.173077")))
        print(f"  B = 4 at r = 0.6722244, a = 0.5314125: theta >= {str(c1)[:14]}")
        c2 = T_lower(5, Fraction(6, 10), a0)
        res.append(("B = 5 at r = 0.6", c2 < Decimal("1.173077")))
        print(f"  B = 5 at r = 0.6, a = a0: {str(c2)[:14]}")
        c3 = T_lower(5, r0, a0, drop_last=True)
        res.append(("the term I(2r, 2B) dropped", c3 < Decimal("1.173077")))
        print(f"  without I(2r, 2B): {str(c3)[:14]}")
        _, _, th5 = f_sup(5)
        res.append(("a claimed 1.1730773", th5 < 1.1730773))
        for U, label in (([0, 1, 2], "U = {0, 1, 2}"), ([0, 1, 4], "U = {0, 1, 4}")):
            s2, d2, q2, v2 = ghr(U)
            print(f"  [Ru96] {label}: |U+U| = {s2}, |U-U| = {d2}, q = {q2}: {str(v2)[:10]}")
            res.append((f"[Ru96] {label}", v2 < Decimal("1.0792")))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
