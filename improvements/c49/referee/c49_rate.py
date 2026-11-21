#!/usr/bin/env python3
"""Exact integer rate checks: rho = a + sqrt(b*c), bound rho^(1/52); decimal brackets by integer powers."""
from math import isqrt
from decimal import Decimal, getcontext

getcontext().prec = 80


def ge_bound(a, D, num, den, k=52):
    """True iff (num/den)^k <= a + sqrt(D)  (exact)."""
    X, S = num ** k, den ** k
    lhs = X - a * S          # (x - a) * S
    if lhs <= 0:
        return True
    return lhs * lhs <= D * S * S


def lt_bound(a, D, num, den, k=52):
    """True iff a + sqrt(D) < (num/den)^k  (exact; D not a square)."""
    Y, S = num ** k, den ** k
    rhs = Y - a * S
    return rhs > 0 and D * S * S < rhs * rhs


def report(name, b, a, c, lo, hi, digits):
    D = b * c
    assert isqrt(D) ** 2 != D
    den = 10 ** digits
    r = (Decimal(a) + Decimal(D).sqrt()) ** (Decimal(1) / Decimal(52))
    print("%s: N00=%d N01=N10=%d N11=%d total=%d" % (name, b, a, c, b + 2 * a + c))
    print("   b*c = %d" % D)
    print("   rho^(1/52) ~ %s" % str(r)[:40])
    print("   exact: %s <= rate : %s ;  rate < %s : %s" % (
        lo, ge_bound(a, D, int(lo.replace('.', '')), 10 ** (len(lo) - 2)),
        hi, lt_bound(a, D, int(hi.replace('.', '')), 10 ** (len(hi) - 2))))
    return a, D


if __name__ == "__main__":
    t = report("Tilley", 4397311488, 7211663100, 11504366352,
               "1.56786638820497592522", "1.56786638820497592523", 20)
    print("   b*c matches note's 50588282321810251776:", t[1] == 50588282321810251776)
    m = report("X* (one set)", 4397440944, 7211751967, 11504366352,
               "1.5678667956347", "1.5678667956348", 13)
    print("   b*c matches claimed 50589771631060716288:", m[1] == 50589771631060716288)
    # tighter bracket of the claimed decimal 1.5678667956347499661...
    print("   tighter: ", ge_bound(m[0], m[1], 15678667956347499661, 10 ** 19),
          lt_bound(m[0], m[1], 15678667956347499662, 10 ** 19))
    m2 = report("X* + X' (two sets)", 4397570400, 7211840834, 11504366352,
                "1.567867203055880727", "1.567867203055880728", 18)
    print("   b*c matches claimed 50591260940311180800:", m2[1] == 50591260940311180800)
    # comparisons
    print("one-set rate > Tilley rate (a'>a and D'>D):", m[0] > t[0] and m[1] > t[1])
    print("two-set rate > one-set rate:", m2[0] > m[0] and m2[1] > m[1])
    # previous public rates per note page 3
    for nm, a, D in (("Koshta", 64728, 3018723123), ("Pinho", 64728, 3019700660)):
        r = (Decimal(a) + Decimal(D).sqrt()) ** (Decimal(1) / Decimal(26))
        print("%s rate %s" % (nm, str(r)[:22]))
