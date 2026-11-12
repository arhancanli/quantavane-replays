#!/usr/bin/env python3
"""Exact replay of the 69a lower-bound row [Tao2022] (teorth/optimizationproblems, main 87b4ebb): C_69 >= 1, from
f(z) = z^n - 1 and lambda_0 = 1, "the critical points are at 0 and lie on the boundary of D(lambda_0, 1)".

For n = 2..60, exactly (integer coefficient lists): f(1) = 0; every zero of f has |z| = 1 (z^n = 1 gives |z|^n = 1);
f' = n z^(n-1), so the critical points are 0 only (multiplicity n - 1); |0 - lambda_0|^2 = 1. Hence for every R < 1 the
disk D(lambda_0, R) holds no critical point, and C_69 >= 1 (the definition takes the infimum of admissible R).
Red controls (--controls): R = 0.999 is not admissible for this f (refused); f(z) = z^n - 2 has a zero outside the
closed unit disk and is refused as an example; lambda_0 = 0 is not a zero of z^n - 1.
usage: python3 -I sendov_example_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


def poly(n, c=1):
    """z^n - c as a coefficient list, constant term first."""
    return [-c] + [0] * (n - 1) + [1]


def ev(p, z):
    return sum(a * z ** k for k, a in enumerate(p))


def deriv(p):
    return [k * a for k, a in enumerate(p)][1:]


def example_ok(n, c=1, lam=1, R=None):
    f = poly(n, c)
    if ev(f, lam) != 0:
        return False, "lambda_0 is not a zero"
    # zeros of z^n - c have |z|^n = |c|: inside the closed unit disk iff |c| <= 1
    if abs(c) > 1:
        return False, "zeros of modulus |c|^(1/n) > 1"
    d = deriv(f)
    crit_only_zero = all(a == 0 for a in d[:-1]) and d[-1] == n      # f' = n z^(n-1)
    dist2 = abs(0 - lam) ** 2
    if R is not None:                                                  # is some critical point within R?
        return dist2 <= R * R, "distance^2 %s vs R^2 %s" % (dist2, R * R)
    return crit_only_zero and dist2 == 1, "f' = %d z^%d, |0 - 1|^2 = %s" % (n, n - 1, dist2)


def main():
    ok = True
    print("69a row [Tao2022]: C_69 >= 1 from f(z) = z^n - 1, lambda_0 = 1")
    for n in range(2, 61):
        g, why = example_ok(n)
        if n <= 4 or n == 60 or not g:
            print("  n = %2d: %s (%s)" % (n, g, why))
        ok &= g
    print("  (n = 2..60 all checked)")
    print("RESULT [Tao2022]: the only critical point 0 is at distance exactly 1 from lambda_0 = 1, so no R < 1 is "
          "admissible and C_69 >= 1: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g, why = example_ok(5, R=F(999, 1000))
        print("  R = 0.999 for z^5 - 1: %s (%s)" % ("admissible (UNEXPECTED)" if g else "not admissible, as expected", why))
        c_ok &= not g
        g, why = example_ok(5, c=2, lam=1)
        print("  f = z^5 - 2: %s (%s)" % ("accepted (UNEXPECTED)" if g else "refused, as expected", why))
        c_ok &= not g
        g, why = example_ok(5, lam=0)
        print("  lambda_0 = 0: %s (%s)" % ("accepted (UNEXPECTED)" if g else "refused, as expected", why))
        c_ok &= not g
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
