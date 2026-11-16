#!/usr/bin/env python3
"""41b, lower bound [Romik2016] 1.64495521: the area of Romik's ambidextrous sofa, program A (exact rationals).

Romik (arXiv 1606.08111, abstract and Theorem 5, read through its HTML rendering): the area is X + arctan Y, X and Y
the real solutions of x^2 (x + 3) = 8 and x (4 x^2 + 3) = 1; equivalently (his eq. (5))
    cbrt(3 + 2 sqrt 2) + cbrt(3 - 2 sqrt 2) - 1 + arctan[(cbrt(sqrt 2 + 1) - cbrt(sqrt 2 - 1))/2].
That the sofa has this area is Romik's theorem (read, not checked). Here: the cubics x^3 + 3x^2 - 8 and 4x^3 + 3x - 1
each have exactly one real root (the first is increasing for x > 0 and negative for x <= 0; the second is increasing
everywhere), enclosed by bisection in exact rationals; arctan Y by halving and the alternating series; the value
compared with the printed 1.64495521 (and the numerical 1.64495 of [Gibbs2014]). The radical form is program B's.
Red controls (--control): the claim 1.64495522 (above the value); Y from 4x^3 + 3x + 1 (a sign slip); X from
x^2 (x + 3) = 9. Standard library only.
"""
import sys
from fractions import Fraction as F

CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None


def root(p, lo, hi, steps=150):
    assert p(lo) < 0 < p(hi)
    for _ in range(steps):
        m = (lo + hi) / 2
        if p(m) < 0: lo = m
        else: hi = m
    return lo, hi


def atan_iv(x):
    """arctan of a rational in (0, 1), enclosed: halving to below 1/8 (with rational square-root bounds), series"""
    from math import isqrt
    def sq(q, up):
        n = q.numerator * 10 ** 80; d = q.denominator
        r = isqrt(n // d) if not up else isqrt(-(-n // d)) + 1
        return F(r, 10 ** 40)
    out = []
    for up in (0, 1):
        y = x; m = 1
        while y > F(1, 8):
            y = y / (1 + sq(1 + y * y, not up)); m *= 2      # upper end: smaller denominator
            y = F(round(y * 10 ** 60) + (1 if up else -1), 10 ** 60)
        s = F(0); p = y; k = 0
        while p > F(1, 10 ** 60):
            s += (-1) ** k * p / (2 * k + 1); p *= y * y; k += 1
        tail = p / (2 * k + 1)
        out.append(m * ((s + tail) if up else (s - tail)))
    return out[0], out[1]


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def main():
    c8 = 9 if CONTROL == "X" else 8
    py = (lambda x: 4 * x ** 3 + 3 * x + 1) if CONTROL == "Y" else (lambda x: 4 * x ** 3 + 3 * x - 1)
    px = lambda x: x ** 3 + 3 * x ** 2 - c8
    X = root(px, F(0), F(3))
    try:
        Y = root(py, F(0), F(1))
    except AssertionError:
        Y = root(py, F(-1), F(0))
    print("X (x^2 (x + 3) = %d) in [%s, %s]; Y (4x^3 + 3x %s 1 = 0) in [%s, %s]"
          % (c8, dec(X[0], 30), dec(X[1], 30), "+" if CONTROL == "Y" else "-", dec(Y[0], 30), dec(Y[1], 30)))
    if Y[0] > 0:
        at = (atan_iv(Y[0])[0], atan_iv(Y[1])[1])
    else:
        at = (-atan_iv(-Y[0])[1], -atan_iv(-Y[1])[0])
    lo, hi = X[0] + at[0], X[1] + at[1]
    print("area X + arctan Y in [%s, %s]" % (dec(lo, 32), dec(hi, 32)))
    printed = F("1.64495522") if CONTROL == "claim" else F("1.64495521")
    tr = printed <= lo and hi < printed + F(1, 10 ** 8)
    g = F("1.64495") <= lo and hi < F("1.64496")
    print("printed %s: a truncation of the value (so a valid lower bound): %s; [Gibbs2014] 1.64495 likewise: %s" % (printed, tr, g))
    ok = tr and g
    print("key %s" % dec(lo, 25))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
