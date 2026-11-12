#!/usr/bin/env python3
"""Replay of the printed value of the 87a lower-bound row [Odl1990] under GRH (teorth/optimizationproblems, main
87b4ebb): "215.333...", the GRH-conditional Odlyzko-Serre bound for root discriminants of totally real fields.

The closed form is the one Hajir, Maire and Ramakrishna state (arXiv 1901.04354v1, the page's [HMR2019], e-print sha256
19ba78bc...d929, lines 295-297 of Cutting_towers_arxiv.tex, read as text): 8 pi e^(gamma + pi/2) ~ 215.33 for totally
real fields and 8 pi e^gamma ~ 44.763 for totally complex ones. Odlyzko's 1990 survey itself was not machine-readable
to us. The bounds themselves are theorems (read); we check the printed decimals:
  1. Euler's constant gamma enclosed exactly: gamma = H_n - ln n - 1/(2n) + sum_{k=1..m} B_2k/(2k n^2k) + R with
     |R| <= |B_(2m+2)|/((2m+2) n^(2m+2)) (Euler-Maclaurin; n = 1000, m = 12; the Bernoulli numbers computed exactly),
     ln 1000 correctly rounded at 70 digits and widened;
  2. pi from Machin's formula with alternating-series bounds (exact rationals); exponentials correctly rounded at 70
     digits and widened by 1e-60;
  3. 8 pi e^(gamma + pi/2) = 215.33251595...: the row's "215.333..." is not a truncation of it; as printed it lies
     4.8e-4 ABOVE the bound (a lower bound rounded up). Also printed for reference: 4 pi e^(1 + gamma) = 60.8395...
     (the standard form of the unconditional bound; the row's 60.8395... is its truncation), and the two values the
     page's notes quote for totally complex fields: 8 pi e^gamma = 44.7632... (44.763, a truncation) and
     4 pi e^gamma = 22.3816... (22.382, rounded up).
Red controls (--controls): gamma with the Euler-Maclaurin remainder dropped and n = 10 misses the 30-digit value by more
than 1e-15 (refused as an enclosure); 215.333 as a lower bound for 8 pi e^(gamma + pi/2) is refused; 215.332 is accepted
as its truncation.
usage: python3 -I odl1990_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from math import comb

CONTROLS = "--controls" in sys.argv
getcontext().prec = 70
W = Decimal(10) ** -60


def bernoulli(nmax):
    B = [F(1)]
    for m in range(1, nmax + 1):
        B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    return B                      # B[1] = -1/2 convention


def dec_lo(x):
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        return Decimal(x.numerator) / Decimal(x.denominator)


def dec_hi(x):
    with localcontext() as c:
        c.rounding = ROUND_CEILING
        return Decimal(x.numerator) / Decimal(x.denominator)


def gamma_interval(n=1000, m=12, remainder=True):
    B = bernoulli(2 * m + 2)
    H = sum(F(1, k) for k in range(1, n + 1))
    s = H - F(1, 2 * n) + sum(B[2 * k] / (2 * k * F(n) ** (2 * k)) for k in range(1, m + 1))
    R = abs(B[2 * m + 2]) / ((2 * m + 2) * F(n) ** (2 * m + 2)) if remainder else F(0)
    ln_n = Decimal(n).ln()
    return dec_lo(s - R) - ln_n - W, dec_hi(s + R) - ln_n + W


def arctan_bounds(x, terms):
    s, lo, hi = F(0), None, None
    for k in range(terms):
        s += F((-1) ** k) * x ** (2 * k + 1) / (2 * k + 1)
        if k % 2 == 0:
            hi = s
        else:
            lo = s
    return lo, hi


def pi_interval():
    l5, h5 = arctan_bounds(F(1, 5), 50)
    l239, h239 = arctan_bounds(F(1, 239), 16)
    return dec_lo(16 * l5 - 4 * h239), dec_hi(16 * h5 - 4 * l239)


def expi(lo, hi):
    return lo.exp() - W, hi.exp() + W


def mul(a, b):
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        lo = a[0] * b[0]
        c.rounding = ROUND_CEILING
        hi = a[1] * b[1]
    return lo, hi           # all factors positive here


def main():
    ok = True
    g = gamma_interval()
    p = pi_interval()
    print("87a row [Odl1990] (GRH): the printed 215.333... against 8 pi e^(gamma + pi/2)")
    print("  gamma in [%s, %s]" % (str(g[0])[:42], str(g[1])[:42]))
    print("  pi    in [%s, %s]" % (str(p[0])[:42], str(p[1])[:42]))
    ok &= str(g[0])[:32] == str(g[1])[:32]
    e1 = expi(g[0] + p[0] / 2, g[1] + p[1] / 2)
    grh_real = mul((8 * p[0], 8 * p[1]), e1)
    eg = expi(*g)
    grh_cplx = mul((8 * p[0], 8 * p[1]), eg)
    unc_real = mul((4 * p[0], 4 * p[1]), expi(1 + g[0], 1 + g[1]))
    unc_cplx = mul((4 * p[0], 4 * p[1]), eg)
    for name, v in (("8 pi e^(gamma + pi/2)", grh_real), ("4 pi e^(1 + gamma)", unc_real),
                    ("8 pi e^gamma", grh_cplx), ("4 pi e^gamma", unc_cplx)):
        print("  %-22s in [%s, %s]" % (name, str(v[0])[:34], str(v[1])[:34]))
    trunc_ok = Decimal("215.332") <= grh_real[0] and grh_real[1] < Decimal("215.333")
    row_as_lower = Decimal("215.333") <= grh_real[0]
    print("  the row's 215.333...: a valid lower bound as printed: %s; above the bound by %.4e; the truncation is "
          "215.3325...: %s" % (row_as_lower, Decimal("215.333") - grh_real[1], trunc_ok))
    print("  [for reference] 60.8395... is the truncation of 4 pi e^(1 + gamma): %s; the notes' 44.763 truncates "
          "8 pi e^gamma: %s; the notes' 22.382 is above 4 pi e^gamma (rounded up): %s"
          % (str(unc_real[0]).startswith("60.8395") and str(unc_real[1]).startswith("60.8395"),
             str(grh_cplx[0]).startswith("44.763"), Decimal("22.382") > unc_cplx[1]))
    ok &= trunc_ok and not row_as_lower
    print("RESULT [Odl1990] GRH: 8 pi e^(gamma + pi/2) = %s...; the printed 215.333... is rounded up past the bound "
          "(it should read 215.3325...): %s" % (str(grh_real[0])[:13], "FINDING CONFIRMED" if ok else "CHECK FAILED"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g2 = gamma_interval(n=10, m=0, remainder=False)
        miss = abs(g2[0] - g[0]) > Decimal(10) ** -15
        print("  gamma from n = 10 with no correction terms: %s (off by %.2e): %s"
              % (str(g2[0])[:20], abs(g2[0] - g[0]), "refused as an enclosure, as expected" if miss else "UNEXPECTED"))
        c_ok &= miss
        print("  215.333 as a lower bound for 8 pi e^(gamma + pi/2): %s" % ("holds (UNEXPECTED)" if row_as_lower else
                                                                         "refused, as expected"))
        c_ok &= not row_as_lower
        acc = Decimal("215.332") <= grh_real[0]
        print("  215.332 as a lower bound: %s" % ("accepted, as expected" if acc else "refused (UNEXPECTED)"))
        c_ok &= acc
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
