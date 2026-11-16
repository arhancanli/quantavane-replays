#!/usr/bin/env python3
"""11a, upper bound [IS2024] C_11a <= pi/2 - delta, delta ~ 0.00013, program B (written separately from program A).

delta = int_{1/2}^{1} min{4 p^3 (1-p)^3, (2p-1)^4 p (1-p)} / (45 - 3 (p (1-p))^3) dp / sqrt(p (1-p))
(arXiv 2407.04835, Section 2.4). Here in the variable s = arcsin(2p - 1) in [0, pi/2] (dp/sqrt(p(1-p)) = ds,
p(1-p) = cos^2 s/4), where the integrand is cos^2 s sin^4 s/(4 E) below s* = arctan(1/sqrt 2) and cos^6 s/(16 E)
above it, E = 45 - 3 cos^6 s/64. A different exact route from program A: every power is expanded into a cosine
polynomial, cos^n s = 2^-n sum_j C(n, j) cos((n - 2j) s), sin^4 s = (3 - 4 cos 2s + cos 4s)/8, products by
cos a cos b = (cos(a+b) + cos(a-b))/2, then integrated term by term; sin(2j s*) = sin(2 s*) U_(j-1)(cos 2 s*) with
cos 2s* = 1/3 and sin 2s* = 2 sqrt 2/3 (Chebyshev recurrence in exact rationals). 1/E = (1/45) sum_k (cos^6 s/960)^k,
12 terms, tail bounded by 960^-12/(1 - 1/960) times the first bracket. Numbers: pi by the Gauss-Legendre iteration in
60-digit decimals, s* = arctan(1/sqrt 2) by Euler's series, all moved outward by 10^-45. Also, for information, a
plain midpoint sum in the variable p (floats, not a proof) as a check of the substitution. Red controls (--control):
45 -> 44; the crossover at sin^2 = 1/4; the claim delta ~ 0.000133. Standard library only.
"""
import sys
from fractions import Fraction as F
from math import comb, sqrt
from decimal import Decimal as D, getcontext

getcontext().prec = 60
PAD = D(10) ** -45
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
D0 = 44 if CONTROL == "den" else 45


def cospow(n):
    """cos^n s as {frequency: coefficient}"""
    out = {}
    for j in range(n + 1):
        f = abs(n - 2 * j)
        out[f] = out.get(f, F(0)) + F(comb(n, j), 2 ** n)
    return out


def mul(a, b):
    out = {}
    for fa, ca in a.items():
        for fb, cb in b.items():
            for f in (fa + fb, abs(fa - fb)):
                out[f] = out.get(f, F(0)) + ca * cb / 2
    return out


SIN4 = {0: F(3, 8), 2: F(-4, 8), 4: F(1, 8)}


def main():
    # crossover data: sin^2 s* = 1/3 (control: 1/4); cos 2s* = 1 - 2 sin^2; sin 2s* = 2 sin cos = c2 * sqrt(rad)
    s2 = F(1, 4) if CONTROL == "cross" else F(1, 3)
    c2s = 1 - 2 * s2                      # cos 2s*
    # sin 2s* = 2 sqrt(s2 (1 - s2)): 2 sqrt(2/9) = (2/3) sqrt 2 ; control: 2 sqrt(3/16) = (1/2) sqrt 3
    coef2, rad = (F(2, 3), 2) if CONTROL != "cross" else (F(1, 2), 3)
    U = [F(1), 2 * c2s]                   # U_0, U_1 at cos 2s*
    for _ in range(200):
        U.append(2 * c2s * U[-1] - U[-2])
    K = 12
    # delta = (1/D0) sum_k (3/(64 D0))^k [ (1/4) int_0^s* sin^4 cos^(2+6k) + (1/16) int_s*^(pi/2) cos^(6+6k) ]
    # an integral of sum_f c_f cos(f s) over [0, s*] = c_0 s* + sum_{f>0} c_f sin(f s*)/f; over [s*, pi/2] = c_0 (pi/2 - s*)
    # - sum_{f>0} c_f sin(f s*)/f (f even, so sin(f pi/2) = 0). sin(f s*) = sin(2 j s*) = coef2 sqrt(rad) U_(j-1).
    A = {"s*": F(0), "pi": F(0), "rad": F(0)}
    first = None
    for k in range(K):
        w = F(3, 64 * D0) ** k / D0
        p1 = mul(SIN4, cospow(2 + 6 * k)); p2 = cospow(6 + 6 * k)
        br = {"s*": F(0), "pi": F(0), "rad": F(0)}
        for f, c in p1.items():
            c = c / 4
            if f == 0: br["s*"] += c
            else: br["rad"] += c * coef2 * U[f // 2 - 1] / f
        for f, c in p2.items():
            c = c / 16
            if f == 0: br["pi"] += c / 2; br["s*"] -= c
            else: br["rad"] -= c * coef2 * U[f // 2 - 1] / f
        if first is None:
            first = dict(br)
        for key in A: A[key] += w * br[key]
    # numbers
    def gl_pi():
        a = D(1); b = 1 / D(2).sqrt(); t = D(1) / 4; p = D(1)
        for _ in range(8):
            an = (a + b) / 2; b = (a * b).sqrt(); t -= p * (a - an) ** 2; a = an; p *= 2
        return (a + b) ** 2 / (4 * t)
    PI = gl_pi()
    def atan_euler(y):
        z = y * y / (1 + y * y); t = y / (1 + y * y); s = D(0); n = 0
        while t > D(10) ** -58:
            s += t; n += 1; t = t * z * (2 * n) / (2 * n + 1)
        return s + t / (1 - z) / 2
    SR = D(rad).sqrt()
    SS = PI / 6 if CONTROL == "cross" else atan_euler(1 / D(2).sqrt())
    def val(dct):
        return sum((D(c.numerator) / D(c.denominator)) * x for c, x in ((dct["s*"], SS), (dct["pi"], PI), (dct["rad"], SR)))
    delta = val(A)
    tail = abs(val(first)) / D0 * (D(3) / (64 * D0)) ** K / (1 - D(3) / (64 * D0))
    lo, hi = delta - PAD, delta + tail + PAD
    print("delta = a s* + b pi + c sqrt%d with a = %.6e, b = %.6e, c = %.6e (series to K = %d, tail <= %.1e)"
          % (rad, float(A["s*"]), float(A["pi"]), float(A["rad"]), K, float(tail)))
    print("delta in [%s, %s]" % (format(lo, ".30f"), format(hi, ".30f")))
    print("pi/2 - delta in [%s, %s]" % (format(PI / 2 - hi - PAD, ".26f"), format(PI / 2 - lo + PAD, ".26f")))
    # information only: midpoint sum in p (floats), to catch a slip in the substitution
    n = 200000; acc = 0.0
    for i in range(n):
        u = (i + 0.5) / n                 # p = 1/2 + (1/2)(1 - (1-u)^2): dense near p = 1 where 1/sqrt blows up
        p = 0.5 + 0.5 * (1 - (1 - u) ** 2)
        dp = (1 - u) / n
        q = p * (1 - p)
        acc += min(4 * q ** 3, (2 * p - 1) ** 4 * q) / (D0 - 3 * q ** 3) / sqrt(q) * dp
    print("for information: a midpoint sum in p gives %.12f" % acc)
    cs = "0.000133" if CONTROL == "claim" else "0.00013"
    half = D(1) / (2 * 10 ** len(cs.split(".")[1]))
    rd = abs(lo - D(cs)) <= half and abs(hi - D(cs)) <= half
    print("the row's delta ~ %s: correct rounding %s%s" % (cs, rd, (" ; delta >= %s: %s" % (cs, lo >= D(cs))) if CONTROL == "claim" else ""))
    ok = rd and (lo >= D(cs) if CONTROL == "claim" else True) and abs(acc - float(delta)) < 1e-8
    print("key %s" % format(lo, ".40f")[:24])
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
