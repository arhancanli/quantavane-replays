#!/usr/bin/env python3
"""19a [E1956]: Esseen's lower bound c_E = (sqrt 10 + 3)/(6 sqrt(2 pi)) ~ 0.4097321837, program A.

Definition used (Esseen 1956, as the row describes it: attained asymptotically by normalized Bernoulli sums). For a
two-point law with P(X = 1) = p, q = 1 - p, standardized, the lattice span is h = 1/sqrt(pq), the third moment
mu_3 = (q - p)/sqrt(pq) and beta_3 = E|X|^3 = (p^2 + q^2)/sqrt(pq). Esseen's lattice expansion gives
sqrt(n) sup_x |F_n(x) - Phi(x)| -> (|mu_3|/6 + h/2) phi(0), so the ratio to beta_3 tends to
g(p) = (|1 - 2p|/6 + 1/2) / ((p^2 + q^2) sqrt(2 pi)), and c_E = max_p g(p).

This program: with u = |1 - 2p| in [0, 1], p^2 + q^2 = (1 + u^2)/2 and g = f(u)/(3 sqrt(2 pi)), f(u) = (u + 3)/(1 + u^2).
f'(u) = (1 - 6u - u^2)/(1 + u^2)^2 changes sign once, at u* = sqrt 10 - 3; f(u*) is computed exactly in Q(sqrt 10)
and equals (sqrt 10 + 3)/2. The decimal value uses integer fixed point (scale 10^-60) with isqrt and pi by Machin's
formula with explicit error terms. Standard library only.
"""
import sys
from fractions import Fraction as Fr
from math import isqrt

S = 10 ** 60


class Q10:
    """a + b sqrt(10), a and b rational"""
    def __init__(self, a, b=0): self.a, self.b = Fr(a), Fr(b)
    def __add__(s, o): o = o if isinstance(o, Q10) else Q10(o); return Q10(s.a + o.a, s.b + o.b)
    def __sub__(s, o): o = o if isinstance(o, Q10) else Q10(o); return Q10(s.a - o.a, s.b - o.b)
    def __mul__(s, o):
        o = o if isinstance(o, Q10) else Q10(o)
        return Q10(s.a * o.a + 10 * s.b * o.b, s.a * o.b + s.b * o.a)
    def inv(s):
        n = s.a * s.a - 10 * s.b * s.b
        return Q10(s.a / n, -s.b / n)
    def __truediv__(s, o): o = o if isinstance(o, Q10) else Q10(o); return s * o.inv()
    def __eq__(s, o): o = o if isinstance(o, Q10) else Q10(o); return s.a == o.a and s.b == o.b
    def sign(s):
        # sign of a + b sqrt 10, exactly
        a, b = s.a, s.b
        if b == 0: return (a > 0) - (a < 0)
        if a == 0: return (b > 0) - (b < 0)
        if (a > 0) == (b > 0): return 1 if a > 0 else -1
        # opposite signs: compare a^2 with 10 b^2
        return (1 if a > 0 else -1) * (1 if a * a > 10 * b * b else -1)


def atan_inv(n):
    """S * atan(1/n) as an integer interval. Terms floor((S g / n^(2k+1)) / (2k+1)) at a guard scale g = 10^5: the
    repeated floors make each term low by less than 3 units; the alternating tail is below the first omitted term."""
    g = 10 ** 5
    p = S * g // n
    acc = 0
    k = 0
    while p:
        t = p // (2 * k + 1)
        acc += t if k % 2 == 0 else -t
        k += 1
        p //= n * n
    return (acc - 3 * k - 3) // g, (acc + 3 * k + 3) // g + 1


def pi_iv():
    a = atan_inv(5); b = atan_inv(239)
    return 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]


def main():
    a = sys.argv[1:]
    ctrl = a[a.index("--control") + 1] if "--control" in a else None
    printed = "0.4097321838" if ctrl == "print" else "0.4097321837"
    print("program A (exact algebra in Q(sqrt 10), fixed point at 10^-60)%s" % ((" CONTROL " + ctrl) if ctrl else ""))
    # f(u) = (u + w)/(1 + u^2) with w = 3 (h/2 term) or, in a control, the h/2 term doubled
    lat = Fr(1, 2) if ctrl != "lattice" else Fr(1)
    # g(u) = (u/6 + lat) / ((1 + u^2)/2) / sqrt(2 pi) = (u + 6 lat)/(3 (1 + u^2)) / sqrt(2 pi)
    w = 6 * lat
    # f'(u) numerator: (1 + u^2) - 2u(u + w) = 1 - 2 w u - u^2 ; root in (0, 1): u* = -w + sqrt(w^2 + 1)
    disc = w * w + 1
    if ctrl == "beta":
        # control: the factor p^2 + q^2 of beta_3 dropped: g = (u/6 + lat) phi(0), increasing in u, max at u = 1
        best = Q10(Fr(1, 6) + lat)
        val_over = best          # times phi(0) = 1/sqrt(2 pi)
        num_den = None
    else:
        assert disc == 10 or ctrl == "lattice"
        if disc == 10:
            ustar = Q10(-w, 1)
            # self-test: u* in (0, 1) and f' changes sign there (exact)
            assert ustar.sign() > 0 and (Q10(1) - ustar).sign() > 0
            fstar = (ustar + w) / (Q10(1) + ustar * ustar)
            # f is increasing on [0, u*] and decreasing on [u*, 1]: compare with the end points
            f0 = Q10(w); f1 = Q10((1 + w) / 2)
            assert (fstar - f0).sign() > 0 and (fstar - f1).sign() > 0
            print("u* = sqrt 10 - 3, p* = (1 - u*)/2 = (4 - sqrt 10)/2; f(u*) = %s + %s sqrt 10" % (fstar.a, fstar.b))
            assert fstar == Q10(Fr(3, 2), Fr(1, 2)), "f(u*) = (sqrt 10 + 3)/2"
            print("f(u*) = (sqrt 10 + 3)/2 exactly; f(0) = 3 and f(1) = 2 are smaller")
            val_over = fstar / 3     # c_E = val_over / sqrt(2 pi)
        else:
            # control: w = 6, disc = 37; u* = sqrt 37 - 6; evaluate numerically through Fractions
            from decimal import Decimal, getcontext
            getcontext().prec = 50
            us = Decimal(37).sqrt() - 6
            fs = (us + 6) / (1 + us * us)
            val_over = Q10(Fr(str(fs / 3)))
    # numeric: c = (a + b sqrt 10) / sqrt(2 pi)
    s10 = isqrt(10 * S * S)                       # floor(S sqrt 10)
    num_lo = val_over.a * S + val_over.b * s10     # b >= 0 here
    num_hi = val_over.a * S + val_over.b * (s10 + 1)
    p0, p1 = pi_iv()
    r_lo = isqrt(2 * p0 * S)                      # S sqrt(2 pi) lower
    r_hi = isqrt(2 * p1 * S) + 1
    c_lo = Fr(num_lo) / r_hi
    c_hi = Fr(num_hi) / r_lo
    k = 40
    print("c_E in [%s, %s]" % (dec(c_lo, k), dec(c_hi, k)))
    d = Fr(printed)
    ulp = Fr(1, 10 ** 10)
    rounding = abs(c_lo - d) <= ulp / 2 and abs(c_hi - d) <= ulp / 2
    trunc = d <= c_lo and c_hi < d + ulp
    print("row 19a [E1956]: c_E = %s...; page %s: %s" % (dec(c_lo, 30), printed,
          ", ".join([x for x, y in (("truncation", trunc), ("correct rounding", rounding)) if y]) or "NONE"))
    good = rounding
    print("VERDICT " + (("REFUSED" if not good else "NOT REFUSED") if ctrl else ("PASS" if good else "FAILED")))


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


if __name__ == "__main__":
    main()
