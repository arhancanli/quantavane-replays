#!/usr/bin/env python3
"""Exact replay of the 51a lower-bound row [CH1964] (teorth/optimizationproblems, main 87b4ebb): B >= 4/7 ~ 0.5714285714
(Clunie and Hayman, "scaling-identity construction").

The construction as He and Tang recall it (arXiv 2602.12217v1, e-print sha256 5134482a...4f7d, section 1, read as
text; Clunie-Hayman's paper was not read): for K > 1, f(z) = sum_{n>=0} (-1)^(n(n-1)/2) K^(-n(n+1)/2) z^n, with the Laurent series k(z) over all n in Z
satisfying k(Kz) = z k(-z); then beta(f) = 1/max_{|z|=1} |k(z)| (Clunie-Hayman p. 156; He-Tang's Theorem for general
phases; read, not checked). Pairing n with -n-1 (as in our [R2026b] replay, c51_check.py) gives, for this phase
(epsilon = -1), max_{|z|=1} |k| = max over x in [-1, 1] of |Q(x)| with
    Q(x) = 2 sum_{n>=0} (-1)^(T_n) K^(-T_n) T_(2n+1)(x),   T_n = n(n+1)/2,  T_m the Chebyshev polynomial.
At K = 6 we check, exactly in rational arithmetic:
  1. Q_N, the sum truncated at N = 8 (degree 15, integer Chebyshev coefficients), and the tail bound
     |Q - Q_N| <= 2 * 6^(-T_8) / (1 - 6^(-9)) on [-1, 1];
  2. with c = 7/4 - tail: Q_N - c and Q_N + c have no root in [-1, 1] (Sturm sequences over Q), Q_N(0) - c < 0 and
     Q_N(0) + c > 0, so max |Q| < 7/4 and beta(f) > 4/7;
  3. a second route, in floating point with an explicit Lipschitz margin: max over a grid of 20001 points of [0, pi/2]
     plus L h / 2 with L = 2 sum (2n+1) 6^(-T_n) (|P'| <= L), as a cross-check, and the value 1/max |Q| = 0.57465...
The printed 0.5714285714 stands next to the exact 4/7.
Red controls (--controls): K = 3 (max |Q| = 1.944... > 7/4) is refused by the Sturm test; the phase epsilon = +1 at
K = 6 (max = Q(1) = 2.33...) is refused; a claimed 0.5750 at K = 6 (above 1/max |Q|) is refused.
usage: python3 -I ch1964_check.py [--controls]      (Python 3.9+, standard library only)
"""
import math
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def chebyshev(m):
    T0, T1 = [F(1)], [F(0), F(1)]
    if m == 0:
        return T0
    for _ in range(m - 1):
        T2 = [F(0)] + [2 * c for c in T1]
        for i, c in enumerate(T0):
            T2[i] -= c
        T0, T1 = T1, trim(T2)
    return T1


def padd(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def pscale(a, c):
    return trim([c * x for x in a])


def ev(p, x):
    r = F(0)
    for c in reversed(p):
        r = r * x + c
    return r


def deriv(p):
    return trim([i * c for i, c in enumerate(p)][1:])


def prem(a, b):
    a = trim(a)
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]
        k = len(a) - len(b)
        for i, y in enumerate(b):
            a[i + k] -= c * y
        a = trim(a)
    return a


def sturm_roots(p, lo, hi):
    """Number of distinct real roots of p in (lo, hi] by Sturm's theorem."""
    seq = [trim(p), deriv(p)]
    while seq[-1]:
        r = prem(seq[-2], seq[-1])
        seq.append([-c for c in r])
    seq = [s for s in seq if s]

    def changes(x):
        vals = [ev(s, x) for s in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for i in range(1, len(vals)) if (vals[i] > 0) != (vals[i - 1] > 0))
    return changes(lo) - changes(hi)


def Q_poly(K, N, eps_sign=-1):
    Q = []
    for n in range(N):
        T = n * (n + 1) // 2
        Q = padd(Q, pscale(chebyshev(2 * n + 1), 2 * F(eps_sign) ** T / F(K) ** T))
    return Q


def tail(K, N):
    TN = N * (N + 1) // 2
    return 2 * F(1) / F(K) ** TN / (1 - F(1) / F(K) ** (N + 1))


def certify(K, target, eps_sign=-1, N=8):
    Q = Q_poly(K, N, eps_sign)
    c = target - tail(K, N)
    up = padd(Q, [-c])
    dn = padd(Q, [c])
    r_up = sturm_roots(up, F(-1), F(1)) + (1 if ev(up, F(-1)) == 0 else 0)
    r_dn = sturm_roots(dn, F(-1), F(1)) + (1 if ev(dn, F(-1)) == 0 else 0)
    ok = r_up == 0 and r_dn == 0 and ev(up, F(0)) < 0 and ev(dn, F(0)) > 0
    return ok, len(Q) - 1, r_up, r_dn, c


def grid_max(K, eps_sign=-1, N=12, M=20000):
    best, arg = 0.0, 0.0
    for j in range(M + 1):
        th = (math.pi / 2) * j / M
        s = sum(2 * (eps_sign ** (n * (n + 1) // 2)) * K ** (-(n * (n + 1) // 2)) * math.cos((2 * n + 1) * th)
                for n in range(N))
        if abs(s) > best:
            best, arg = abs(s), th
    L = sum(2 * (2 * n + 1) * K ** (-(n * (n + 1) // 2)) for n in range(N + 5))
    return best, arg, best + L * (math.pi / 2) / M / 2 + 1e-9


def main():
    ok = True
    print("51a row [CH1964]: B >= 4/7 from the Clunie-Hayman series with K = 6 (phase epsilon = -1)")
    good, deg, r_up, r_dn, c = certify(6, F(7, 4))
    print("  1. Q_8 of degree %d; tail bound %s (< 1e-16: %s)" % (deg, float(tail(6, 8)), tail(6, 8) < F(1, 10 ** 16)))
    print("  2. Sturm: roots of Q_8 - c in [-1, 1]: %d, of Q_8 + c: %d; signs at 0 as needed: %s; "
          "so max |Q| < 7/4: %s" % (r_up, r_dn, good, good))
    ok &= good
    best, arg, upper = grid_max(6)
    print("  3. floating-point grid: max |Q| = %.9f at theta = %.6f; with the Lipschitz margin <= %.9f < 1.75: %s; "
          "1/max = %.9f" % (best, arg, upper, upper < 1.75, 1 / best))
    ok &= upper < 1.75
    print("RESULT [CH1964]: beta(f) = 1/max |k| > 4/7 at K = 6 (1/max |k| = %.6f), so B > 4/7 = 0.5714285714...: %s"
          % (1 / best, "HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g, _, r_up, _, _ = certify(3, F(7, 4))
        print("  K = 3 (max |Q| = %.6f): %s (roots of Q - c: %d)" % (grid_max(3)[0], "certified (UNEXPECTED)" if g else
                                                                    "refused, as expected", r_up))
        c_ok &= not g
        g, _, r_up, _, _ = certify(6, F(7, 4), eps_sign=1)
        print("  phase +1 at K = 6 (max = Q(1) = %.6f): %s" % (grid_max(6, eps_sign=1)[0], "certified (UNEXPECTED)" if g
                                                               else "refused, as expected"))
        c_ok &= not g
        g, _, _, _, _ = certify(6, 1 / F("0.5750"))
        print("  a claimed 0.5750 at K = 6 (needs max |Q| < 1/0.575 = %.6f): %s"
              % (1 / 0.575, "certified (UNEXPECTED)" if g else "refused, as expected"))
        c_ok &= not g
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
