#!/usr/bin/env python3
"""71a, lower bound [Hod2017] C >= beta(1/2) > 6.4547837 (Hod, arXiv 1711.00762, Theorem 4.4), beta from its
definition, program A (exact rational intervals).

As read in the paper (Section 4): h(p) = -p log2 p - (1 - p) log2(1 - p); q_0 = z, q_(m+1) = 1/(1 + q_m);
pi_m = q_1 ... q_m (pi_0 = 1), pi_(-1) = 1/z, pi_(-2) = 1/(1 - z);
    beta(z) = 4 z (1 - z) sum_(k >= -2) h(q_(k+2)) pi_k,
and Theorem 4.4 takes z = 1/2 (the limit-of-limits construction from the identity function; read, not checked).
Here: the recursion in exact rationals for k <= K = 220, with q_k <= 2/3 for k >= 1 checked along the way, so the tail
sum_(k > K) h(q_(k+2)) pi_k <= sum_(k > K) pi_k <= 3 pi_(K+1) (h <= 1); log2 by our own ln (2 atanh series in fixed-point
integers at 2^-256, floor and ceiling at every step, tails bounded). Red controls (--control): the claim 6.4547838; z = 0.499; q_(m+1) = 1/(2 + q_m).
Standard library only.
"""
import sys
from fractions import Fraction as F

CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
P = 256
ONE = 1 << P


def fx_lo(q): return (q.numerator << P) // q.denominator          # floor(q 2^P)
def fx_hi(q): return -((-q.numerator << P) // q.denominator)      # ceil(q 2^P)


def ln_fx(x):
    """ln of a rational x > 0 as integers (lo, hi) at scale 2^-P: x = 2^k y, y in [1, 2), ln y = 2 atanh(u),
    u = (y - 1)/(y + 1) <= 1/3; series in fixed point with floor/ceil at every step and the tail bounded"""
    k = 0
    while x >= 2: x /= 2; k += 1
    while x < 1: x *= 2; k -= 1
    def at2(u):
        ul, uh = fx_lo(u), fx_hi(u)
        u2l, u2h = (ul * ul) >> P, -((-uh * uh) >> P)
        sl = sh = 0; pl, ph = ul, uh; j = 0
        while ph > 0 and j < 400:
            sl += pl // (2 * j + 1); sh += -((-ph) // (2 * j + 1))
            pl = (pl * u2l) >> P; ph = -((-ph * u2h) >> P); j += 1
            if ph < 4:
                break
        # tail: sum_{i >= j} u^(2i+1)/(2i+1) <= ph / (1 - u^2) <= 2 ph (u <= 1/3), plus one unit per step of rounding
        return 2 * sl, 2 * (sh + 2 * ph + j + 2)
    lo, hi = at2((x - 1) / (x + 1))
    l2lo, l2hi = at2(F(1, 3))
    return (lo + k * (l2lo if k >= 0 else l2hi), hi + k * (l2hi if k >= 0 else l2lo))


L2 = ln_fx(F(2))


def h_iv(p):
    """binary entropy of a rational p in (0, 1), enclosed as Fractions: -(p ln p + (1-p) ln(1-p))/ln 2"""
    a = ln_fx(p); b = ln_fx(1 - p)
    num_lo = -(p * F(a[1], ONE) + (1 - p) * F(b[1], ONE)); num_hi = -(p * F(a[0], ONE) + (1 - p) * F(b[0], ONE))
    return num_lo / F(L2[1], ONE), num_hi / F(L2[0], ONE)


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def main():
    z = F(499, 1000) if CONTROL == "z" else F(1, 2)
    add = 2 if CONTROL == "rec" else 1
    K = 220
    q = [z]
    for m in range(K + 3):
        q.append(1 / (add + q[-1]))
    bounded = all(x <= F(2, 3) for x in q[1:]) or CONTROL == "rec"
    pi = {-2: 1 / (1 - z), -1: 1 / z, 0: F(1)}
    for m in range(1, K + 2):
        pi[m] = pi[m - 1] * q[m]
    lo = hi = F(0)
    for k in range(-2, K + 1):
        a, b = h_iv(q[k + 2])
        lo += a * pi[k]; hi += b * pi[k]
    tail = 3 * pi[K + 1] if CONTROL != "rec" else 2 * pi[K + 1]
    c = 4 * z * (1 - z)
    lo, hi = c * lo, c * (hi + tail)
    print("q_k <= 2/3 for 1 <= k <= %d: %s; tail after k = %d at most %.1e" % (K + 2, bounded, K, float(c * tail)))
    print("beta(%s) in [%s, %s]" % (z, dec(lo, 30), dec(hi, 30)))
    claim = F("6.4547838") if CONTROL == "claim" else F("6.4547837")
    ok = bounded and lo > claim
    print("beta > %s (the row): %s (margin %.3e)" % (claim, lo > claim, float(lo - claim)))
    print("key %s" % dec(lo, 25))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
