#!/usr/bin/env python3
"""11a, upper bound [IS2024] C_11a <= pi/2 - delta, delta ~ 0.00013: delta from its definition, program A.

Ivanisvili and Stone (arXiv 2407.04835, Section 2.4, read through its HTML rendering) give
    delta = int_{1/2}^{1} min{4 p^3 (1-p)^3, (2p-1)^4 p (1-p)} / (45 - 3 (p (1-p))^3) dp / sqrt(p (1-p)).
Their theorem (that pi/2 - delta bounds C_11a) is read, not checked; we compute delta. With p = (1 + sin t)/2,
t in [0, pi/2]: p (1-p) = cos^2 t/4 =: q, 2p - 1 = sin t, dp/sqrt(p (1-p)) = dt; the two terms of the min are
cos^6 t/16 and sin^4 t cos^2 t/4, equal where tan^2 t = 1/2, i.e. t* = arctan(1/sqrt 2) (sin t* = 1/sqrt 3,
cos t* = sqrt(2/3)); the second is the smaller on [0, t*]. With 1/(45 - 3q^3) = (1/45) sum_k (cos^6 t/960)^k,
    delta = (1/45) sum_{k>=0} 960^-k [ (1/4) int_0^t* sin^4 cos^(2+6k) + (1/16) int_t*^(pi/2) cos^(6+6k) ],
and sin^4 = (1 - cos^2)^2 reduces everything to J_n(a, b) = int_a^b cos^n, by J_n = [cos^(n-1) sin]/n
+ (n-1)/n J_(n-2), J_0 = b - a. So each truncation is exactly r1 + r2 sqrt 2 + r3 t* + r4 pi with rationals r_i
(computed in Python Fractions, kept in that basis); the tail after K terms is at most 960^-K/(1 - 1/960) times the
k = 0 bracket. Then pi (Machin), sqrt 2 (integer square roots) and t* = arctan(1/sqrt 2) (halving and the
alternating series) enclose delta. Red controls (--control): the denominator 45 replaced by 44; the crossover taken
at tan^2 t = 1/3; the claim delta >= 0.000133 (so pi/2 - 0.000133 claimed as the bound). Standard library only.
"""
import sys
from fractions import Fraction as F
from math import isqrt

CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
D0 = 44 if CONTROL == "den" else 45


class Q2:
    """r1 + r2 sqrt2 + r3 ts + r4 pi with rational r_i (ts = the crossover angle)"""
    def __init__(self, a=0, b=0, c=0, d=0): self.v = [F(a), F(b), F(c), F(d)]
    def __add__(s, o): return Q2(*[x + y for x, y in zip(s.v, o.v)])
    def __sub__(s, o): return Q2(*[x - y for x, y in zip(s.v, o.v)])
    def scale(s, r): return Q2(*[x * r for x in s.v])


# crossover: tan^2 t* = 1/2 gives sin^2 = 1/3, cos^2 = 2/3; sin cos^(n-1) = (2/3)^((n-2)/2) * sqrt(2)/3 for even n
# (control "cross": tan^2 = 1/3, sin^2 = 1/4, cos^2 = 3/4; then sin cos^(n-1) = (3/4)^((n-2)/2) * sqrt(3)/4: we keep the
# sqrt 2 slot for that algebraic number and adjust its value below)
CS2 = F(3, 4) if CONTROL == "cross" else F(2, 3)
SC1 = F(1, 4) if CONTROL == "cross" else F(1, 3)        # (sin t* cos t*)/sqrt(2) or /sqrt(3), squared factor


def boundary(n, at):
    """cos^(n-1) sin at 0, t* or pi/2, for even n >= 2, as a Q2"""
    if at in ("0", "pi/2"):
        return Q2()
    # cos^(n-1) sin = cos^(n-2) (sin cos) = CS2^((n-2)/2) * sin cos; sin cos = sqrt(2)/3 (or sqrt(3)/4 in the control)
    return Q2(0, CS2 ** ((n - 2) // 2) * (F(1, 3) if CONTROL != "cross" else F(1, 4)), 0, 0)


def J(n, a, b, memo):
    """int_a^b cos^n for even n, a, b in {'0', 't*', 'pi/2'}"""
    key = (n, a, b)
    if key in memo:
        return memo[key]
    if n == 0:
        val = {("0", "t*"): Q2(0, 0, 1, 0), ("t*", "pi/2"): Q2(0, 0, -1, F(1, 2))}[(a, b)]
    else:
        val = (boundary(n, b) - boundary(n, a)).scale(F(1, n)) + J(n - 2, a, b, memo).scale(F(n - 1, n))
    memo[key] = val
    return val


def main():
    K = 14
    memo = {}
    tot = Q2()
    first = None
    for k in range(K):
        m = 2 + 6 * k
        s1 = J(m, "0", "t*", memo) - J(m + 2, "0", "t*", memo).scale(2) + J(m + 4, "0", "t*", memo)   # sin^4 cos^m
        s2 = J(6 + 6 * k, "t*", "pi/2", memo)
        br = s1.scale(F(1, 4)) + s2.scale(F(1, 16))
        if first is None:
            first = br
        tot = tot + br.scale(F(3, 64 * D0) ** k)
    tot = tot.scale(F(1, D0))
    print("delta (series to K = %d) = r1 + r2 sqrt2 + r3 t* + r4 pi with r = (%s, %s, %s, %s) (denominators of %d, %d, %d, %d digits)"
          % (K, *[("%.6e" % float(x)) for x in tot.v], *[len(str(x.denominator)) for x in tot.v]))
    # enclosures: pi, sqrt 2 (or sqrt 3 in the control), t* = arctan(1/sqrt 2) (or arctan(1/sqrt 3) = pi/6)
    def atan_series(x, terms):
        s = F(0); p = x
        for kk in range(terms):
            s += (-1) ** kk * p / (2 * kk + 1); p *= x * x
        nxt = p / (2 * terms + 1)
        return (s, s + nxt) if terms % 2 == 0 else (s - nxt, s)
    a5 = atan_series(F(1, 5), 50); a239 = atan_series(F(1, 239), 24)
    PI = (16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0])
    rad = 3 if CONTROL == "cross" else 2
    R = (F(isqrt(rad * 10 ** 120), 10 ** 60), F(isqrt(rad * 10 ** 120) + 1, 10 ** 60))
    # t* = arctan(1/sqrt2) = 2 arctan(x) with x = (1/sqrt2)/(1 + sqrt(3/2)) = 1/(sqrt2 + sqrt3); for the control pi/6
    if CONTROL == "cross":
        TS = (PI[0] / 6, PI[1] / 6)
    else:
        S3 = (F(isqrt(3 * 10 ** 120), 10 ** 60), F(isqrt(3 * 10 ** 120) + 1, 10 ** 60))
        xl = 1 / (R[1] + S3[1]); xh = 1 / (R[0] + S3[0])
        TS = (2 * atan_series(xl, 60)[0], 2 * atan_series(xh, 60)[1])
    def enclose(q):
        lo = q.v[0]; hi = q.v[0]
        for coef, iv in zip(q.v[1:], (R, TS, PI)):
            a, b = coef * iv[0], coef * iv[1]
            lo += min(a, b); hi += max(a, b)
        return lo, hi
    lo, hi = enclose(tot)
    f1 = enclose(first)
    tail = f1[1] / D0 * F(3, 64 * D0) ** K / (1 - F(3, 64 * D0))
    hi += tail
    print("delta in [%s, %s] (series tail <= %.1e)" % (dec(lo, 30), dec(hi, 30), float(tail)))
    bl, bh = PI[0] / 2 - hi, PI[1] / 2 - lo
    print("pi/2 - delta in [%s, %s]" % (dec(bl, 26), dec(bh, 26)))
    cs = "0.000133" if CONTROL == "claim" else "0.00013"
    claim = F(cs); half = F(1, 2 * 10 ** len(cs.split(".")[1]))
    ok_round = abs(lo - claim) <= half and abs(hi - claim) <= half
    ok_side = lo >= claim if CONTROL == "claim" else True
    print("the row's delta ~ %s: correct rounding %s%s" % (cs, ok_round, (" ; delta >= %s: %s" % (claim, lo >= claim)) if CONTROL == "claim" else ""))
    ok = ok_round and ok_side
    print("key %s" % dec(lo, 22))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


def dec(q, k):
    neg = q < 0; q = abs(q)
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return ("-" if neg else "") + s[:-k] + "." + s[-k:]


if __name__ == "__main__":
    main()
