#!/usr/bin/env python3
"""7a, upper bound [ZZ2020] mu(pi) <= 7.103205334137...: the bound from Zeilberger and Zudilin's final formula,
program A (exact rational intervals). Source: arXiv 1912.06345, read through its HTML rendering (Proposition 2,
eq. (18), Lemma 3, eq. (10), and the section "World record"); their programs were not read or run.

  N_1, N_2, N_3: the roots of 108 N^3 - 2359989 N^2 + 138304 N - 2048 = 0 (eq. (18)), N_3 real (21851.69...), N_1, N_2
  complex conjugates; c = pi/(2 sqrt 3) - log(3 sqrt 3/4) (= psi(2/3) - psi(1/2), eq. (10));
  -tau = log|N_1| - (5/2) log 2 + 4 - c,  sigma = log N_3 - (5/2) log 2 + 4 - c,  mu(pi) <= 1 + sigma/tau.
Here: the cubic has one real root (its discriminant is negative), N_3 by bisection; |N_1|^2 = N_1 N_2 = 2048/(108 N_3)
(the product of the roots), so log|N_1| = (log(2048/108) - log N_3)/2 with no complex arithmetic; logs by atanh series,
pi by Machin, sqrt 3 by integer square roots; every quantity an exact rational interval. Also the identity
psi(2/3) - psi(1/2) = c from Gauss's digamma values (psi(2/3) = -gamma - (3/2) log 3 + pi/(2 sqrt 3),
psi(1/2) = -gamma - 2 log 2), which is exact algebra: (3/2) log 3 - 2 log 2 = log(3 sqrt 3/4).
Not checked: the arithmetic lemmas (denominators, the prime set P_n) and the saddle-point asymptotics, which give the
formula. Red controls (--control): the claim 7.103205334136 (below the value); the constant term 2048 replaced by 2047;
the (5/2) log 2 replaced by 2 log 2. Standard library only.
"""
import sys
from fractions import Fraction as F
from math import isqrt

CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
EPS = F(1, 10 ** 50)


def ln_iv(x):
    k = 0
    while x >= 2: x /= 2; k += 1
    while x < 1: x *= 2; k -= 1
    def at2(u):
        s = F(0); p = u; j = 0
        while p > EPS:
            s += p / (2 * j + 1); p *= u * u; j += 1
        return 2 * s, 2 * (s + p / (1 - u * u))
    lo, hi = at2((x - 1) / (x + 1))
    l2 = at2(F(1, 3))
    return (lo + k * (l2[0] if k >= 0 else l2[1]), hi + k * (l2[1] if k >= 0 else l2[0]))


def atan_series(x, terms):
    s = F(0); p = x
    for k in range(terms):
        s += (-1) ** k * p / (2 * k + 1); p *= x * x
    nxt = p / (2 * terms + 1)
    return (s, s + nxt) if terms % 2 == 0 else (s - nxt, s)


a5 = atan_series(F(1, 5), 50); a239 = atan_series(F(1, 239), 24)
PI = (16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0])
S3 = (F(isqrt(3 * 10 ** 120), 10 ** 60), F(isqrt(3 * 10 ** 120) + 1, 10 ** 60))


def dec(q, k):
    neg = q < 0; q = abs(q)
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return ("-" if neg else "") + s[:-k] + "." + s[-k:]


def main():
    d0 = 2047 if CONTROL == "const" else 2048
    a, b, c_, d = 108, -2359989, 138304, -d0
    p = lambda x: ((a * x + b) * x + c_) * x + d
    disc = 18 * a * b * c_ * d - 4 * b ** 3 * d + b * b * c_ * c_ - 4 * a * c_ ** 3 - 27 * a * a * d * d
    print("cubic 108 N^3 - 2359989 N^2 + 138304 N - %d: discriminant %s 0, so one real root and a complex pair"
          % (d0, "<" if disc < 0 else ">="))
    lo, hi = F(21000), F(23000)
    assert p(lo) < 0 < p(hi)
    for _ in range(200):
        m = (lo + hi) / 2
        if p(m) < 0: lo = m
        else: hi = m
    print("N_3 in [%s, %s]" % (dec(lo, 25), dec(hi, 25)))
    prod = F(d0, a)                      # N_1 N_2 N_3 = -d/a
    l3 = (ln_iv(lo)[0], ln_iv(hi)[1])
    lp = ln_iv(prod)
    l1 = ((lp[0] - l3[1]) / 2, (lp[1] - l3[0]) / 2)
    n1 = (prod / hi, prod / lo)
    print("|N_1|^2 = %d/(108 N_3) in [%s, %s], |N_1| ~ %.12f" % (d0, dec(n1[0], 20), dec(n1[1], 20), float(n1[0]) ** 0.5))
    k = F(2) if CONTROL == "log2" else F(5, 2)
    L2 = ln_iv(F(2))
    lc = ln_iv(3 * S3[0] / 4)[0], ln_iv(3 * S3[1] / 4)[1]
    cc = (PI[0] / (2 * S3[1]) - lc[1], PI[1] / (2 * S3[0]) - lc[0])
    print("c = pi/(2 sqrt 3) - log(3 sqrt 3/4) in [%s, %s] (paper: 0.64527561...)" % (dec(cc[0], 25), dec(cc[1], 25)))
    # psi(2/3) - psi(1/2) = pi/(2 sqrt 3) - (3/2) log 3 + 2 log 2: exact equality with c, checked in intervals
    L3 = ln_iv(F(3))
    alt = (PI[0] / (2 * S3[1]) - F(3, 2) * L3[1] + 2 * L2[0], PI[1] / (2 * S3[0]) - F(3, 2) * L3[0] + 2 * L2[1])
    same = abs(alt[0] - cc[0]) < F(1, 10 ** 40) and abs(alt[1] - cc[1]) < F(1, 10 ** 40)
    print("psi(2/3) - psi(1/2) from Gauss's values equals c: %s" % same)
    common = (-k * L2[1] + 4 - cc[1], -k * L2[0] + 4 - cc[0])
    tau = (-(l1[1] + common[1]), -(l1[0] + common[0]))
    sig = (l3[0] + common[0], l3[1] + common[1])
    print("lim sup log|I'_n|/n = log|N_1| - (5/2) log 2 + 4 - c in [%s, %s]" % (dec(-tau[1], 20), dec(-tau[0], 20)))
    print("lim log b'_n/n = log N_3 - (5/2) log 2 + 4 - c in [%s, %s]" % (dec(sig[0], 20), dec(sig[1], 20)))
    mu = (1 + sig[0] / tau[1], 1 + sig[1] / tau[0])
    print("1 + sigma/tau in [%s, %s]" % (dec(mu[0], 30), dec(mu[1], 30)))
    printed = F("7.103205334136") if CONTROL == "claim" else F("7.103205334137")
    holds = mu[1] <= printed + F(1, 10 ** 12) and mu[0] >= printed     # the row prints the leading digits
    print("row %s: the leading digits of the value: %s (the value is %s above it, so as an upper bound the safe"
          " rounding is %s)" % (printed, holds, dec(mu[0] - printed, 18), dec(printed + F(1, 10 ** 12), 12)))
    ok = disc < 0 and same and holds
    print("key %s" % dec(mu[0], 25))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
