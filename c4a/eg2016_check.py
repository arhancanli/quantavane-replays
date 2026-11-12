#!/usr/bin/env python3
"""Replay of the finite counting step of the 4a upper-bound row [EG2016] (teorth/optimizationproblems, main 87b4ebb):
C_4a <= 2.756 (Ellenberg-Gijswijt).

Ellenberg and Gijswijt bound a cap A in F_3^n by |A| <= 3 M_n, where M_n is the number of monomials prod x_i^(a_i) with
every a_i in {0, 1, 2} and sum a_i <= 2n/3 (the polynomial-method lemma; read, not checked). The exponential base of
M_n is the finite content we check, exactly in rational arithmetic:
  1. for every 0 < x <= 1, M_n <= sum_a x^(sum a - 2n/3) = (1 + x + x^2)^n x^(-2n/3); at x = 5931/10000 this gives
     M_n^3 x^(2n) <= (1 + x + x^2)^(3n), and (1 + x + x^2)^3 / x^2 < 2.756^3, so M_n <= 2.756^n;
  2. M_n computed exactly (dynamic programming over the coordinate sum) for n = 3, 6, ..., 600, and both
     inequalities checked as integer comparisons; (3 M_n)^(1/n) is printed, increasing towards the base;
  3. the optimal base Gamma = min_{0<x<1} (1 + x + x^2)/x^(2/3), attained at x0 = (sqrt 33 - 1)/8 (the root of
     4x^2 + x - 2), enclosed: Gamma = 2.755104613..., below the printed 2.756 (the row's value is Gamma rounded up,
     valid as an upper bound).
So r_3(F_3^n)^(1/n) <= (3 M_n)^(1/n) <= 3^(1/n) * 2.756 and C_4a <= 2.756.
Red controls (--controls): x = 1/2 gives (1 + x + x^2)^3/x^2 = 343/16 > 2.756^3 (refused at that x); a claimed 2.755
is below Gamma and is refused; M_n with the degree cap 2n/3 replaced by n (all monomials of degree <= n) exceeds
2.756^n at n = 600 (refused).
usage: python3 -I eg2016_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv
getcontext().prec = 50


def M(n, cap_num=2, cap_den=3):
    """Number of a in {0,1,2}^n with sum a <= (cap_num/cap_den) n."""
    cnt = [1]                      # cnt[s] = number of vectors with coordinate sum s
    for _ in range(n):
        new = [0] * (len(cnt) + 2)
        for s, c in enumerate(cnt):
            new[s] += c
            new[s + 1] += c
            new[s + 2] += c
        cnt = new
    lim = (cap_num * n) // cap_den
    return sum(cnt[: lim + 1])


def base_ok(x, B=F(2756, 1000)):
    return (1 + x + x * x) ** 3 / (x * x) < B ** 3


def main():
    ok = True
    x = F(5931, 10000)
    b = base_ok(x)
    print("4a row [EG2016]: the base of M_n (monomials of degree <= 2n/3 with exponents <= 2)")
    print("  1. x = 5931/10000: (1 + x + x^2)^3/x^2 = %.9f < 2.756^3 = %.9f: %s"
          % ((1 + x + x * x) ** 3 / (x * x), 2.756 ** 3, b))
    ok &= b
    bad = []
    for n in range(3, 601, 3):
        m = M(n)
        chern = m ** 3 * x ** (2 * n) <= (1 + x + x * x) ** (3 * n)
        within = m <= F(2756, 1000) ** n
        if not (chern and within):
            bad.append(n)
        if n in (3, 30, 300, 600):
            r = (Decimal(3 * m).ln() / n).exp()
            print("  2. n = %3d: M_n = %s%s; M_n^3 x^2n <= (1+x+x^2)^3n: %s; M_n <= 2.756^n: %s; (3 M_n)^(1/n) = %s"
                  % (n, str(m)[:12], "..." if len(str(m)) > 12 else "", chern, within, str(r)[:10]))
    print("     all n = 3, 6, ..., 600: %s" % (not bad))
    ok &= not bad
    # Gamma at x0 = (sqrt33 - 1)/8, in high precision (x0 is the critical point: 4x^2 + x - 2 = 0)
    s33 = Decimal(33).sqrt()
    x0 = (s33 - 1) / 8
    crit = 4 * x0 * x0 + x0 - 2
    gam = (1 + x0 + x0 * x0) / ((x0 * x0).ln() / 3).exp()          # x0^(2/3)
    print("  3. x0 = %s (4 x0^2 + x0 - 2 = %.1e); Gamma = %s...; Gamma < 2.756: %s"
          % (str(x0)[:20], crit, str(gam)[:20], gam < Decimal("2.756")))
    ok &= gam < Decimal("2.756") and gam > Decimal("2.7551")
    print("RESULT [EG2016]: M_n <= 2.756^n for every n (the exact inequality at x = 5931/10000), so with |A| <= 3 M_n, "
          "C_4a <= 2.756 (Gamma = %s...): %s" % (str(gam)[:12], "HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g = base_ok(F(1, 2))
        print("  x = 1/2: (1 + x + x^2)^3/x^2 = %s: %s" % ((1 + F(1, 2) + F(1, 4)) ** 3 / F(1, 4),
                                                         "below 2.756^3 (UNEXPECTED)" if g else "refused, as expected"))
        c_ok &= not g
        claim = gam <= Decimal("2.755")
        print("  a claimed C_4a <= 2.755 from this bound: %s" % ("holds (UNEXPECTED)" if claim else "refused, as expected"))
        c_ok &= not claim
        m_all = M(600, 1, 1)
        g = m_all <= F(2756, 1000) ** 600
        print("  degree cap n instead of 2n/3, n = 600: M = %s... <= 2.756^600: %s" % (str(m_all)[:8],
              "UNEXPECTED" if g else "no, refused as expected"))
        c_ok &= not g
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
