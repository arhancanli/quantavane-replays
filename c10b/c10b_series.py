#!/usr/bin/env python3
"""Program 1 for 10b (the complex Grothendieck constant): Haagerup's upper bound [H1987], from its definition, in exact
rationals (hypergeometric series with a rigorous tail; pi by Machin's formula with a rigorous tail).

Definition (Haagerup 1987, Theorem 3.1, as restated in Pisier's survey "Grothendieck's theorem, past and present",
Section 4, and in Heilman-Jones-Malavolta, arXiv 2609.25462v1, (2.5)-(2.6)):
    phi(k) = k int_0^{pi/2} cos^2 t (1 - k^2 sin^2 t)^(-1/2) dt,   K_0 in (0, 1) the unique solution of
    phi(K_0) = (pi/8)(K_0 + 1),   and   K_G^C <= 8 / (pi (K_0 + 1)).
Expanding (1 - k^2 sin^2 t)^(-1/2) and integrating term by term (Wallis), phi(k) = (pi/4) g(k) with
    g(k) = sum_{j >= 0} a_j k^(2j+1),   a_j = C(2j, j)^2 / (16^j (j + 1)),
so the equation is 2 g(k) = k + 1 (pi cancels). The a_j decrease, so the tail after J terms lies in
[0, a_J k^(2J+1) / (1 - k^2)]. f(k) = 2 g(k) - k - 1 is increasing on (0, 1) (a_0 = 1 and a_j > 0); K_0 is
enclosed by bisection on rationals, taking f's lower enclosure at the right end and its upper one at the left.

Also (information): Haagerup's conjectured value (int_0^{pi/2} cos^2 t (1 + sin^2 t)^(-1/2) dt)^(-1), here from
g at k^2 = -1 transformed to k^2 = 1/2 (see conj() below), and Pisier's e^(1 - gamma) [P1978], gamma by
Euler-Maclaurin.

  c10b_series.py [--claim-h DEC] [--claim-p DEC] [--digits D]
Standard library only; run with python3 -I.
"""
import sys
from fractions import Fraction as Fr
from math import comb


def atan_inv(x, eps):
    """atan(1/x) for an integer x >= 2, as an interval [lo, hi] of width <= eps (alternating series)."""
    s, k, term = Fr(0), 0, Fr(1, x)
    while True:
        t = Fr(1, (2 * k + 1) * x ** (2 * k + 1))
        if t < eps:
            # alternating with decreasing terms: the sum lies between s and s + (-1)^k t
            return (s, s + t) if k % 2 == 0 else (s - t, s)
        s += t if k % 2 == 0 else -t
        k += 1


def pi_interval(eps):
    a = atan_inv(5, eps / 64)
    b = atan_inv(239, eps / 64)
    return 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]


def g_interval(k, eps):
    """g(k) = sum a_j k^(2j+1) for 0 < k < 1 (or k^2 given via k2), as [lo, hi]."""
    k2 = k * k
    s, j, a, p = Fr(0), 0, Fr(1), k
    while True:
        tail = a * p / (1 - k2)
        if tail < eps:
            return s, s + tail
        s += a * p
        a = a * Fr((2 * j + 1) ** 2, 4 * (j + 1) * (j + 2))   # a_{j+1} / a_j
        p *= k2
        j += 1


def check_coefficients(J=60):
    """a_j from the recurrence equals C(2j,j)^2/(16^j (j+1)) and the Wallis integrals give the same series."""
    a = Fr(1)
    for j in range(J):
        assert a == Fr(comb(2 * j, j) ** 2, 16 ** j * (j + 1)), j
        # int_0^{pi/2} cos^2 t sin^{2j} t dt = (pi/2) C(2j,j)/4^j * 1/(2j+2); times the binomial coefficient
        # C(-1/2, j)(-1)^j = C(2j,j)/4^j of the expansion, over pi/4:
        w = Fr(comb(2 * j, j), 4 ** j) * Fr(comb(2 * j, j), 4 ** j) * Fr(1, 2 * j + 2) * 2
        assert w == a, j
        a = a * Fr((2 * j + 1) ** 2, 4 * (j + 1) * (j + 2))
    return J


def ln_interval(eps):
    """ln 2 and ln 5 as intervals: ln 2 = 2 atanh(1/3), ln(5/4) = 2 atanh(1/9) (positive series, geometric tails)."""
    def atanh_inv(x):
        s, k = Fr(0), 0
        while True:
            t = Fr(1, (2 * k + 1) * x ** (2 * k + 1))
            tail = t * Fr(x * x, x * x - 1)
            if tail < eps:
                return s, s + tail
            s += t
            k += 1
    l2 = tuple(2 * v for v in atanh_inv(3))
    l54 = tuple(2 * v for v in atanh_inv(9))
    l5 = (2 * l2[0] + l54[0], 2 * l2[1] + l54[1])
    return l2, l5


def gamma_interval():
    """Euler's gamma by Euler-Maclaurin at n = 2000 = 2^4 5^3:
    gamma = H_n - ln n - 1/(2n) + 1/(12 n^2) - 1/(120 n^4) + 1/(252 n^6) - theta/(240 n^8), 0 < theta < 1."""
    n = 2000
    H = sum(Fr(1, i) for i in range(1, n + 1))
    l2, l5 = ln_interval(Fr(1, 10 ** 45))
    lnn = (4 * l2[0] + 3 * l5[0], 4 * l2[1] + 3 * l5[1])
    base = H - Fr(1, 2 * n) + Fr(1, 12 * n ** 2) - Fr(1, 120 * n ** 4) + Fr(1, 252 * n ** 6)
    return base - lnn[1] - Fr(1, 240 * n ** 8), base - lnn[0]


def exp_interval(x, eps):
    """e^x for a rational 0 < x < 1: partial sums below, plus a geometric tail bound above."""
    s, t, k = Fr(0), Fr(1), 0
    while True:
        if t * 2 < eps:          # tail after this point <= t (1 + x/(k+1) + ...) <= 2 t
            return s, s + 2 * t
        s += t
        k += 1
        t = t * x / k


def dec(fr, d, up=False):
    """fr >= 0 to d decimals, rounded down (or up, for the upper end of an enclosure)."""
    q = fr.numerator * 10 ** d // fr.denominator
    if up and q * fr.denominator != fr.numerator * 10 ** d:
        q += 1
    s = str(q).rjust(d + 1, "0")
    return s[:-d] + "." + s[-d:]


def main():
    claim_h = claim_p = None
    D = 30
    a = sys.argv[1:]
    while a:
        if a[0] == "--claim-h":
            claim_h = a[1]; a = a[2:]
        elif a[0] == "--claim-p":
            claim_p = a[1]; a = a[2:]
        elif a[0] == "--digits":
            D = int(a[1]); a = a[2:]
        else:
            raise SystemExit("unknown option " + a[0])
    print("coefficients a_j: recurrence, closed form and Wallis integrals agree for j < %d" % check_coefficients())
    eps = Fr(1, 10 ** (D + 12))
    f = lambda k, side: 2 * g_interval(k, eps)[side] - k - 1
    # a guess by bisection in scaled integers (not rigorous, only a starting point), then an exact bracket
    SC = 10 ** (D + 10)
    def f_approx(kk):                            # kk = k * SC
        k2 = kk * kk // SC
        s, a_num, a_den, p, j = 0, 1, 1, kk, 0
        while p * a_num // a_den > 0:
            s += p * a_num // a_den
            a_num, a_den = a_num * (2 * j + 1) ** 2, a_den * 4 * (j + 1) * (j + 2)
            p = p * k2 // SC
            j += 1
        return 2 * s - kk - SC
    a_, b_ = SC // 2, SC * 99 // 100
    while b_ - a_ > 1:
        m_ = (a_ + b_) // 2
        if f_approx(m_) > 0:
            b_ = m_
        else:
            a_ = m_
    step = Fr(1, 10 ** (D + 4))
    lo, hi = Fr(a_, SC) - step, Fr(a_, SC) + step
    assert f(lo, 1) < 0 < f(hi, 0), "the exact bracket failed"
    print("K_0 in [%s, %s]" % (dec(lo, D), dec(hi, D, True)))
    pl, ph = pi_interval(Fr(1, 10 ** (D + 12)))
    print("pi in [%s, %s] (Machin)" % (dec(pl, D + 5), dec(ph, D + 5, True)))
    Hlo, Hhi = Fr(8) / (ph * (1 + hi)), Fr(8) / (pl * (1 + lo))
    print("Haagerup's bound 8/(pi(K_0+1)) in [%s, %s]" % (dec(Hlo, D), dec(Hhi, D, True)))
    gl = 2 * g_interval(lo, eps)[0]
    print("check: phi(K_0) / (pi/4) = 2 g(K_0)... equals K_0 + 1 within the enclosure: %s"
          % (abs(gl - (1 + lo)) < Fr(1, 10 ** D)))
    ok = True
    if claim_h is not None:
        c = Fr(claim_h)
        holds = Hhi <= c
        print("printed upper bound %s [H1987]: %s (printed minus value >= %.3e)"
              % (claim_h, "holds" if holds else "FAILS", float(c - Hhi)))
        ok &= holds
    # Haagerup's conjecture (information): with m = -1, the transformation K(m) = K(m/(m-1))/sqrt(1-m),
    # E(m) = sqrt(1-m) E(m/(m-1)) gives int cos^2 t (1 + sin^2 t)^(-1/2) dt = 2K(-1) - E(-1) = sqrt(2)(K(1/2) - E(1/2)),
    # and K(1/2) - E(1/2) = (pi/2) sum_j c_j (1/2)^j with c_j = C(2j,j)^2/16^j * (1 - 1/(1-2j)) ... computed here
    # directly as (pi/2) sum_j C(2j,j)^2/16^j (1/2)^j * 2j/(2j-1) (j >= 1), positive terms, geometric tail.
    s, j, t = Fr(0), 1, Fr(1)
    tail = None
    while True:
        t = Fr(comb(2 * j, j) ** 2, 16 ** j) / 2 ** j * Fr(2 * j, 2 * j - 1)
        if t * 2 < eps:          # ratio of consecutive terms < 1/2 for j >= 1
            tail = 2 * t
            break
        s += t
        j += 1
    # sqrt(2) by integer square root
    r = 10 ** (D + 15)
    q = __import__("math").isqrt(2 * r * r)
    s2 = (Fr(q, r), Fr(q + 1, r))
    glo = s2[0] * pl / 2 * s
    ghi = s2[1] * ph / 2 * (s + tail)
    print("information: Haagerup's conjectured value (int cos^2 t (1 + sin^2 t)^(-1/2) dt)^(-1) in [%s, %s]"
          % (dec(1 / ghi, D), dec(1 / glo, D, True)))
    gm = gamma_interval()
    print("information: Euler's gamma in [%s, %s]" % (dec(gm[0], D), dec(gm[1], D, True)))
    e_lo = exp_interval(1 - gm[1], eps)[0]
    e_hi = exp_interval(1 - gm[0], eps)[1]
    print("information: e^(1 - gamma) [P1978] in [%s, %s]" % (dec(e_lo, D), dec(e_hi, D, True)))
    if claim_p is not None:
        c = Fr(claim_p)
        holds = e_hi <= c
        print("printed upper bound %s [P1978]: %s" % (claim_p, "holds" if holds else "FAILS"))
        ok &= holds
    print("VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
