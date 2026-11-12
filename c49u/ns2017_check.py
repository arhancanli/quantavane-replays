#!/usr/bin/env python3
"""Replay of the counting step of the 49a upper-bound row [NS2017] (teorth/optimizationproblems, main 87b4ebb):
C_49 = lim f(n)^(1/n) <= 3/2^(2/3) ~ 1.8898815748 (Naslund-Sawin).

Naslund and Sawin prove |F| <= 3(n + 1) sum_{k <= n/3} C(n, k) for a sunflower-free family F of subsets of [n] (the
slice-rank argument; read, not checked). The base of that sum is the finite content we check, exactly:
  1. sum_{k <= n/3} C(n, k) <= (1 + x)^n x^(-n/3) for every 0 < x <= 1; at x = 1/2 this is (3/2)^n 2^(n/3) =
     (3/2^(2/3))^n, i.e. (sum)^3 <= 27^n / 4^n as an integer inequality, checked for n = 1..900 with exact binomials;
  2. hence f(n) <= 3(n + 1) (3/2^(2/3))^n and C_49 <= 3/2^(2/3); the base is optimal for this sum (x = 1/2 is the
     minimiser of (1 + x) x^(-1/3)), and 3/2^(2/3) = 1.88988157484230974... (cube of 3/2^(2/3) is 27/4 exactly); the
     printed 1.8898815748 is its truncation (and rounding) to 10 decimals.
Red controls (--controls): the cut n/2 in place of n/3 exceeds (3/2^(2/3))^n at n = 900 (refused); x = 2/3 gives the base
(5/3)/(2/3)^(1/3) > 3/2^(2/3) (not the optimum: its cube 375/54 > 27/4); a claimed 1.8898 (below the base) is refused.
usage: python3 -I ns2017_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F
from math import comb

CONTROLS = "--controls" in sys.argv
getcontext().prec = 40


def partial(n, num=1, den=3):
    return sum(comb(n, k) for k in range((num * n) // den + 1))


def main():
    ok = True
    print("49a row [NS2017]: the base of sum_{k <= n/3} C(n, k)")
    bad = [n for n in range(1, 901) if partial(n) ** 3 * 4 ** n > 27 ** n]
    print("  1. (sum_{k<=n/3} C(n,k))^3 <= (27/4)^n for n = 1..900: %s" % (not bad))
    ok &= not bad
    base = (Decimal(27) / 4).ln() / 3
    b = base.exp()
    print("  2. 3/2^(2/3) = (27/4)^(1/3) = %s...; printed 1.8898815748 is its truncation: %s"
          % (str(b)[:20], str(b)[:12] == "1.8898815748"))
    ok &= str(b)[:12] == "1.8898815748"
    x_opt_cube = (1 + F(1, 2)) ** 3 / F(1, 2)
    print("     at x = 1/2: ((1 + x) x^(-1/3))^3 = %s = 27/4: %s" % (x_opt_cube, x_opt_cube == F(27, 4)))
    ok &= x_opt_cube == F(27, 4)
    for n in (30, 300, 900):
        v = (Decimal(3 * (n + 1) * partial(n)).ln() / n).exp()
        print("     n = %d: (3(n+1) sum)^(1/n) = %s" % (n, str(v)[:10]))
    print("RESULT [NS2017]: with |F| <= 3(n+1) sum_{k<=n/3} C(n,k) <= 3(n+1)(27/4)^(n/3), C_49 <= 3/2^(2/3) = "
          "1.8898815748...: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g = partial(900, 1, 2) ** 3 * 4 ** 900 <= 27 ** 900
        print("  cut at n/2, n = 900: %s" % ("within (27/4)^(n/3) (UNEXPECTED)" if g else "exceeds, refused as expected"))
        c_ok &= not g
        c23 = (1 + F(2, 3)) ** 3 / F(2, 3)
        print("  x = 2/3: cube of the base %s > 27/4: %s" % (c23, "refused as the optimum, as expected" if c23 > F(27, 4)
                                                             else "UNEXPECTED"))
        c_ok &= c23 > F(27, 4)
        claim = Decimal("1.8898") >= b
        print("  a claimed C_49 <= 1.8898 from this count: %s" % ("holds (UNEXPECTED)" if claim else "refused, as expected"))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
