#!/usr/bin/env python3
"""41b, lower bound [Romik2016] 1.64495521, program B (written separately from program A): Romik's radical form
    cbrt(3 + 2 sqrt 2) + cbrt(3 - 2 sqrt 2) - 1 + arctan[(cbrt(sqrt 2 + 1) - cbrt(sqrt 2 - 1))/2]
(arXiv 1606.08111, eq. (5)), in 60-digit decimals: square roots from the decimal module, cube roots by Newton's
method with a final bracket check (r^3 <= x <= (r + 10^-55)^3), arctan by the Euler series
arctan y = sum_n 2^(2n) (n!)^2 / (2n + 1)! * y^(2n+1)/(1 + y^2)^(n+1) (all terms positive; the tail is bounded by a
geometric series), every result moved outward by 10^-50. Also checked: the two cube-root combinations satisfy the
cubics of program A's form, u^3 - 3u - 6 = 0 for u = X + 1 and w^3 + 3w - 2 = 0 for w = 2Y, to 10^-45.
Red controls (--control): the claim 1.64495522; 3 + 2 sqrt 2 replaced by 3 + 2 sqrt 3; the arctan term dropped.
Standard library only.
"""
import sys
from decimal import Decimal as D, getcontext

getcontext().prec = 60
PAD = D(10) ** -50
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None


def cbrt(x):
    r = D(float(x) ** (1.0 / 3.0))
    for _ in range(200):
        r2 = r - (r * r * r - x) / (3 * r * r)
        if abs(r2 - r) < D(10) ** -58:
            r = r2; break
        r = r2
    r = r.quantize(D(10) ** -56)
    while r * r * r > x: r -= D(10) ** -56
    while (r + D(10) ** -56) ** 3 < x: r += D(10) ** -56
    assert r ** 3 <= x <= (r + D(10) ** -55) ** 3
    return r, r + D(10) ** -55


def atan_euler(y):
    """arctan y for 0 < y < 1 by Euler's series; returns (lo, hi)"""
    z = y * y / (1 + y * y)
    t = y / (1 + y * y); s = D(0); n = 0
    while t > D(10) ** -58:
        s += t
        n += 1
        t = t * z * (2 * n) / (2 * n + 1)
    tail = t / (1 - z)
    return s - PAD, s + tail + PAD


def main():
    s2 = D(2).sqrt()
    base = 3 + 2 * (D(3).sqrt() if CONTROL == "radical" else s2)
    a = cbrt(base); b = cbrt(3 - 2 * s2)
    c = cbrt(s2 + 1); d = cbrt(s2 - 1)
    u = (a[0] + b[0], a[1] + b[1])
    w = (c[0] - d[1], c[1] - d[0])
    cub1 = abs(u[0] ** 3 - 3 * u[0] - 6) < D(10) ** -45
    cub2 = abs(w[0] ** 3 + 3 * w[0] - 2) < D(10) ** -45
    print("u = X + 1 = cbrt(3 + 2 sqrt 2) + cbrt(3 - 2 sqrt 2) in [%s, %s]; u^3 - 3u - 6 = 0: %s"
          % (format(u[0], ".40f"), format(u[1], ".40f"), cub1))
    print("w = 2Y = cbrt(sqrt 2 + 1) - cbrt(sqrt 2 - 1) in [%s, %s]; w^3 + 3w - 2 = 0: %s"
          % (format(w[0], ".40f"), format(w[1], ".40f"), cub2))
    at = atan_euler(w[0] / 2)[0], atan_euler(w[1] / 2)[1]
    if CONTROL == "atan":
        at = (D(0), D(0))
    lo = u[0] - 1 + at[0] - PAD; hi = u[1] - 1 + at[1] + PAD
    print("area in [%s, %s]" % (format(lo, ".34f"), format(hi, ".34f")))
    printed = D("1.64495522") if CONTROL == "claim" else D("1.64495521")
    tr = printed <= lo and hi < printed + D("0.00000001")
    g = D("1.64495") <= lo and hi < D("1.64496")
    print("printed %s: a truncation of the value: %s; [Gibbs2014] 1.64495 likewise: %s" % (printed, tr, g))
    ok = tr and g and (cub1 and cub2 if CONTROL != "radical" else True)
    print("key %s" % format(lo, ".40f")[:27])
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
