#!/usr/bin/env python3
"""Exact replay of the 3a row [K2026] C_3a > 1.19102809* (kleinwaks, PR #134): the masked-digit limit construction
with M = <1518, 1524, 1587, 2024, 2032, 2116> n [0, 17032] (B = 17032, base 2B + 1 = 34065), from the statement in the
pull request alone (no file of its verification package was read or run).  The PR's bound, for 0 < x < 1:
    C_3a >= 1 + (log P_-(x) - log P_+(x)) / log(2B + 1),
    P_+(x) = sum_{s in M+M} x^s,   P_-(x) = sum_{d in M-M} x^kappa(d),   kappa(d) = min{a + b : a, b in M, a - b = d}.

1. M by reachability in the semigroup (integers); M + M and every kappa(d) by bit sets (kappa(d) = d + 2 min{b in M :
   b + d in M} for d >= 0, and symmetrically); |M|, |M+M|, |M-M| against the PR's 3121, 18730, 32369.  A second route
   for kappa: the minimum over all |M|^2 pairs, on the restriction to M n [0, 4000].
2. The best x located in floating point, then the bound certified at the rational x0 = 9996789/10^7: every power x0^k
   enclosed in 60-digit decimal interval arithmetic (directed rounding), P_- summed from below and P_+ from above,
   logarithms correctly rounded and moved one unit outward, log(34065) from above.
3. Anchor: with the bounded digits {0, ..., 5} (P_+ = 1 + x + ... + x^10, P_- = 1 + 2x + ... + 2x^5, base 11) the same
   formula gives Zheng's 1.173077 ([Z2025]), as the pull request says.
The lemma itself (the PR's proof paper) is read, not checked.
usage: python3 -I k2026_check.py [--controls]
"""
import math
import sys
from decimal import Decimal, getcontext, ROUND_FLOOR, ROUND_CEILING, localcontext

CONTROLS = "--controls" in sys.argv
GENS = [1518, 1524, 1587, 2024, 2032, 2116]
B = 17032
X0 = Decimal(9996789) / Decimal(10 ** 7)
getcontext().prec = 60


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def mask(gens, cap):
    reach = [False] * (cap + 1)
    reach[0] = True
    for v in range(1, cap + 1):
        reach[v] = any(v >= g and reach[v - g] for g in gens)
    return [v for v in range(cap + 1) if reach[v]]


def data(M, cap):
    bits = 0
    for a in M:
        bits |= 1 << a
    S = 0
    for a in M:
        S |= bits << a
    sums = [s for s in range(2 * cap + 1) if (S >> s) & 1]
    kappa = {}
    for d in range(-cap, cap + 1):
        inter = bits & (bits >> d) if d >= 0 else bits & (bits << (-d))
        if inter:
            p = (inter & -inter).bit_length() - 1      # least p in M with p + d in M (d >= 0) or p + d in M (d < 0)
            kappa[d] = 2 * p + d
    return sums, kappa


def kappa_pairs(M):
    out = {}
    for a in M:
        for b in M:
            d = a - b
            if d not in out or a + b < out[d]:
                out[d] = a + b
    return out


def F_float(sums, kv, base, x):
    lx = math.log(x)
    pp = sum(math.exp(s * lx) for s in sums)
    pm = sum(math.exp(k * lx) for k in kv)
    return 1 + (math.log(pm) - math.log(pp)) / math.log(base)


def powers(x, kmax, rounding):
    out = [Decimal(1)]
    with localcontext() as ctx:
        ctx.prec = 60
        ctx.rounding = rounding
        for _ in range(kmax):
            out.append(out[-1] * x)
    return out


def certify(sums, kv, base, x):
    """A rigorous lower bound for F(x)."""
    kmax = max(max(sums), max(kv))
    lo, hi = powers(x, kmax, ROUND_FLOOR), powers(x, kmax, ROUND_CEILING)
    with localcontext() as ctx:
        ctx.prec = 60
        ctx.rounding = ROUND_FLOOR
        pm_lo = sum((lo[k] for k in kv), Decimal(0))
        ctx.rounding = ROUND_CEILING
        pp_hi = sum((hi[s] for s in sums), Decimal(0))
    ulp = Decimal(10) ** -55
    ln_pm_lo = pm_lo.ln() - ulp            # correctly rounded at 60 digits, moved outward
    ln_pp_hi = pp_hi.ln() + ulp
    ln_base_hi = Decimal(base).ln() + ulp
    num = ln_pm_lo - ln_pp_hi
    need(num > 0, "the numerator is not positive")
    return 1 + num / ln_base_hi - ulp


def main():
    ok = True
    M = mask(GENS, B)
    sums, kappa = data(M, B)
    sub = [a for a in M if a <= 4000]
    _, k_sub = data(sub, 4000)
    need(kappa_pairs(sub) == k_sub, "the two routes for kappa disagree on M n [0, 4000]")
    counts = (len(M), len(sums), len(kappa))
    good = counts == (3121, 18730, 32369)
    kv = list(kappa.values())
    base = 2 * B + 1
    xs = [9990 + i / 10 for i in range(100)]
    fbest, xbest = max((F_float(sums, kv, base, x / 10000), x / 10000) for x in xs)
    lo, hi = xbest - 1e-4, xbest + 1e-4
    for _ in range(80):
        m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if F_float(sums, kv, base, m1) < F_float(sums, kv, base, m2):
            lo = m1
        else:
            hi = m2
    fmax = F_float(sums, kv, base, lo)
    cert = certify(sums, kv, base, X0)
    good &= cert > Decimal("1.19102809")
    ok &= good
    print(f"M = <{', '.join(map(str, GENS))}> n [0, {B}]: |M| = {counts[0]}, |M+M| = {counts[1]}, |M-M| = {counts[2]} "
          f"(the pull request: 3121, 18730, 32369); kappa by bit sets equals the pairwise minimum on M n [0, 4000]")
    print(f"  the bound as a function of x peaks near x = {lo:.10f} at {fmax:.12f} (floating point)")
    print(f"  certified at x0 = 9996789/10^7: F(x0) >= {str(cert)[:22]}... (60-digit interval arithmetic)")
    print(f"RESULT [K2026]: C_3a >= F(x0) > 1.19102809 (margin {float(cert - Decimal('1.19102809')):.2e}): "
          f"{'HOLDS' if good else 'FAILS'}")
    zs = list(range(11))
    zk = [0] + [k for k in range(1, 6) for _ in (0, 1)]
    zmax = max(F_float(zs, zk, 11, x / 100000) for x in range(1, 100000))
    print(f"anchor [Z2025]: bounded digits 0..5, base 11: sup_x F = {zmax:.9f} (Zheng's 1.173077)")
    need(abs(zmax - 1.173077) < 1e-6, "Zheng anchor failed")
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        M2 = mask(GENS, B - 1)
        s2, k2 = data(M2, B - 1)
        res.append(("cap 17031", (len(M2), len(s2), len(k2)) != (3121, 18730, 32369)))
        print(f"  cap 17031: |M| = {len(M2)}, |M+M| = {len(s2)}, |M-M| = {len(k2)}")
        M3 = mask([1518, 1524, 1587, 2024, 2032, 2117], B)
        s3, k3 = data(M3, B)
        f3 = max(F_float(s3, list(k3.values()), base, x / 100000) for x in range(99900, 99999))
        res.append(("generator 2117 for 2116", f3 < 1.19102809))
        print(f"  generator 2117 for 2116: the best value near x = 0.999 is {f3:.9f}")
        res.append(("a claimed 1.19102810", fmax < 1.19102810))
        res.append(("x0 = 0.99 (far from the peak)", certify(sums, kv, base, Decimal("0.99")) < Decimal("1.19102809")))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
