#!/usr/bin/env python3
"""Exact replay of the finite part of [V2026] C_22a >= 2.124 (teorth/optimizationproblems #356, Varlamov), written
from the paper alone (paper.pdf @ 6d11254, Sections 4 to 6); its verify.py was neither read nor run.

What the paper reduces the theorem to (thickness 1, length l):
  Lemma 4.2   int_K g(x, y) ds_y <= H(l) := int_2^oo (2/t^3) min{l, (4/3)(t+1)^3} dt,
  Prop. 5.2   H(l) <= A l^(1/3) for every l > 0, with A = k (3/4)^(1/3) and k = 50657/10000,
  Thm 1.1     L >= B C^(3/4) with B = (4 pi / A)^(3/4), i.e. B^4 = 256 pi^3 / (3 k^3), and B > 531/250 = 2.124.

Our route, sharing nothing with the paper's completed-square argument except the statement:
  1. H in closed form from our own primitive of (8/3)(t+1)^3/t^3, checked as an exact identity of Laurent
     polynomials (with the log term kept apart): H(l) = l/4 for l <= 36, and for l = (4/3)(R+1)^3, R >= 2,
     H = h(R) = 4R + 8 log(R/2) + 3 - 4/R; a float quadrature of the defining integral agrees (sanity only);
  2. l <= 36: H = l/4 <= A l^(1/3) is l^2 <= 48 k^3, and l^2 <= 1296 < 48 k^3 exactly;
  3. l >= 36: A l^(1/3) = k (R+1), so we need k (R+1) - h(R) > 0 for R >= 2. On [2, 100] h and k(R+1) both increase,
     so on a dyadic interval [a, b] it suffices that hi(h(b)) < k(a+1), with hi(log) a rational upper bound
     (log 2 from the artanh(1/3) series with a geometric tail, then log x = m log 2 + log(x/2^m)); intervals are
     bisected until this holds. On [100, oo) log(R/2) <= log 50 + (R - 100)/100 (concavity), which leaves a
     linear lower bound with positive slope k - 4 - 8/100 and a positive value at R = 100;
  4. pi enclosed from Gauss's formula pi = 48 arctan(1/18) + 32 arctan(1/57) - 20 arctan(1/239) (the formula itself
     confirmed by Gaussian-integer multiplication: (18+i)^12 (57+i)^8 (239-i)^5 is a positive real multiple of 1+i),
     alternating-series bounds;
  5. B^4 = 256 pi^3/(3 k^3) compared with fourth powers of rationals.
Red controls (each must be refused): k = 50655/10000 (step 3 must find R where k(R+1) < h(R)); the printed bound
moved to 2.1241; a wrong primitive (8 log t replaced by 7 log t) in step 1.

Usage: python3 c22a_check.py [--controls]     (standard library only; a few seconds)
"""
import sys
import time
from fractions import Fraction as F

K = F(50657, 10000)


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ---------- rational bounds for log ----------
def artanh_bounds(y, terms):
    """lo <= artanh(y) <= hi for rational 0 <= y < 1: partial sum, plus the geometric tail bound."""
    s, p = F(0), y
    for j in range(terms):
        s += p / (2 * j + 1)
        p *= y * y
    tail = p / ((2 * terms + 1) * (1 - y * y))   # sum_{j>=terms} y^(2j+1)/(2j+1) <= y^(2terms+1)/((2terms+1)(1-y^2))
    return s, s + tail


LOG2_LO, LOG2_HI = (2 * b for b in artanh_bounds(F(1, 3), 40))   # log 2 = 2 artanh(1/3)


def log_bounds(x):
    """Rational lo <= log x <= hi for rational x >= 1."""
    need(x >= 1, "log_bounds needs x >= 1")
    m = 0
    while x >= 2:
        x /= 2
        m += 1
    lo, hi = artanh_bounds((x - 1) / (x + 1), 40)          # log x = 2 artanh((x-1)/(x+1)), here x in [1, 2)
    return m * LOG2_LO + 2 * lo, m * LOG2_HI + 2 * hi


# ---------- step 1: the closed form of H ----------
def primitive_identity(log_coeff=F(8)):
    """Our primitive P(t) = (8/3) t + c log t - 8/t - 4/(3 t^2); P' must equal (8/3)(t+1)^3/t^3 = (8/3)(1 + 3/t + 3/t^2 + 1/t^3),
    as Laurent polynomials {power: coefficient}. Then H = P(R) - P(2) + l/R^2, l/R^2 = (4/3)(R+1)^3/R^2, must have
    rational part 4R + 3 - 4/R (constants: -(8/3)*2 + 8/2 + 4/12 + 4) and log part 8 log(R/2)."""
    P = {1: F(8, 3), -1: F(-8), -2: F(-4, 3)}          # without the log term
    dP = {k - 1: k * c for k, c in P.items()}
    dP[-1] = dP.get(-1, 0) + log_coeff                   # d/dt (c log t) = c t^-1
    target = {0: F(8, 3), -1: F(8), -2: F(8), -3: F(8, 3)}
    ok = all(dP.get(p, 0) == target.get(p, 0) for p in set(dP) | set(target))
    # rational part of H: P(R) - P(2) + (4/3)(R^3 + 3R^2 + 3R + 1)/R^2
    P2 = F(8, 3) * 2 - F(8, 2) - F(4, 3) / 4             # P(2) without its log term
    rat = dict(P)
    for p, c in {1: F(4, 3), 0: F(4), -1: F(4), -2: F(4, 3)}.items():
        rat[p] = rat.get(p, 0) + c
    rat[0] = rat.get(0, 0) - P2
    want = {1: F(4), 0: F(3), -1: F(-4)}
    ok2 = all(rat.get(p, 0) == want.get(p, 0) for p in set(rat) | set(want))
    return ok and ok2 and log_coeff == 8   # the log part is c(log R - log 2) = 8 log(R/2) only for c = 8


def h_float(R):
    import math
    return 4 * R + 8 * math.log(R / 2) + 3 - 4 / R


def H_quadrature(l, steps=200000):
    """Float sanity check of the defining integral (not part of the proof): substitution t = 2/u, u in (0, 1]."""
    tot = 0.0
    for i in range(steps):
        u = (i + 0.5) / steps
        t = 2 / u
        tot += (2 / t ** 3) * min(l, 4 / 3 * (t + 1) ** 3) * (2 / u ** 2) / steps
    return tot


# ---------- step 3: k(R+1) - h(R) > 0 for R >= 2 ----------
def h_hi(b):
    return 4 * b + 8 * log_bounds(b / 2)[1] + 3 - F(4) / b


def h_lo(a):
    return 4 * a + 8 * log_bounds(a / 2)[0] + 3 - F(4) / a


def prop52(k, R0=F(100), max_intervals=200000):
    """Returns (True, intervals) if k(R+1) > h(R) on [2, oo), else (False, R) with a point where it fails.
    Breadth first: every interval of one size is tried before any smaller one, and each split point is first tested
    for an outright failure, so a failing k is caught at a grid point instead of being chased toward a crossing."""
    from collections import deque
    queue, done = deque([(F(2), R0)]), 0
    while queue:
        a, b = queue.popleft()
        if h_hi(b) < k * (a + 1):
            done += 1
            continue
        mid = (a + b) / 2
        for r in (a, mid):
            if h_lo(r) >= k * (r + 1):        # a point where the inequality fails outright
                return False, r
        if done + len(queue) > max_intervals:
            return False, None
        queue.append((a, mid))
        queue.append((mid, b))
    # tail: gap(R) >= (k - 4 - 8/R0) R + k - 3 - 8 hi(log(R0/2)) + 8 for R >= R0 (tangent line of log at R0, -4/R... dropped)
    slope = k - 4 - F(8) / R0
    at_R0 = (k - 4) * R0 + k - 3 - 8 * log_bounds(R0 / 2)[1]
    need(slope > 0 and at_R0 > 0, "tail bound on [R0, oo) fails")
    return True, done


# ---------- step 4: pi ----------
def arctan_bounds(x, terms):
    s, p, sign = F(0), x, 1
    sums = []
    for j in range(terms + 1):
        s += sign * p / (2 * j + 1)
        sums.append(s)
        p *= x * x
        sign = -sign
    return min(sums[-1], sums[-2]), max(sums[-1], sums[-2])   # alternating, decreasing terms: the value is between


def gauss_formula_holds():
    def mul(z, w):
        return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

    def pw(z, n):
        r = (1, 0)
        for _ in range(n):
            r = mul(r, z)
        return r
    z = mul(mul(pw((18, 1), 12), pw((57, 1), 8)), pw((239, -1), 5))
    return z[0] == z[1] and z[0] > 0       # a positive real multiple of 1 + i: argument exactly pi/4 (mod 2 pi)


def pi_bounds():
    need(gauss_formula_holds(), "Gauss's arctan formula failed the Gaussian-integer check")
    a = arctan_bounds(F(1, 18), 30)
    b = arctan_bounds(F(1, 57), 20)
    c = arctan_bounds(F(1, 239), 15)
    lo = 48 * a[0] + 32 * b[0] - 20 * c[1]
    hi = 48 * a[1] + 32 * b[1] - 20 * c[0]
    need(F(3) < lo < hi < F(4), "pi enclosure out of range")
    return lo, hi


def dec(x, d=14):
    q = x.numerator * 10 ** d // x.denominator
    s = str(q)
    return s[:-d] + "." + s[-d:]


def main():
    t0 = time.time()
    controls = "--controls" in sys.argv
    failed = False
    print("[V2026] C_22a >= 2.124 (#356): exact replay of the finite part, k = %s" % K)
    ok1 = primitive_identity()
    print("1. closed form of H: our primitive and the rational/log parts of h(R) as identities:", ok1)
    for R in (2.0, 3.0, 8.0, 20.0):
        l = 4 / 3 * (R + 1) ** 3
        print(f"   sanity, float quadrature at R = {R:g}: {H_quadrature(l):.9f} vs h(R) = {h_float(R):.9f}")
    print(f"   sanity, l = 20 (below 36): {H_quadrature(20.0):.9f} vs l/4 = 5")
    ok2 = 36 ** 2 < 48 * K ** 3
    print(f"2. l <= 36: l^2 <= 1296 < 48 k^3 = {dec(48 * K ** 3, 6)}: {ok2}")
    ok3, n_int = prop52(K)
    print(f"3. k (R+1) - h(R) > 0 for every R >= 2: {ok3} ({n_int} dyadic intervals on [2, 100], tangent-line tail beyond)")
    plo, phi = pi_bounds()
    print(f"4. pi in [{dec(plo, 30)}, {dec(phi, 30)}] (Gauss's formula confirmed by Gaussian integers)")
    B4lo = 256 * plo ** 3 / (3 * K ** 3)
    B4hi = 256 * phi ** 3 / (3 * K ** 3)
    ok5 = F(531, 250) ** 4 < B4lo
    enc = F(212403897858, 10 ** 11) ** 4 < B4lo and B4hi < F(212403897859, 10 ** 11) ** 4
    print(f"5. B^4 = 256 pi^3/(3 k^3): (531/250)^4 < B^4: {ok5}; 2.12403897858 < B < 2.12403897859: {enc}")
    # Remark 5.3 of the paper (not used in the proof): h(8)/9 > 5.0655, so this method cannot give k <= 5.0655
    h8_lo = F(69, 2) + 16 * LOG2_LO
    rem = h8_lo / 9 > F(50655, 10000)
    print(f"   Remark 5.3: h(8) = 69/2 + 16 log 2 > {dec(h8_lo, 9)}, h(8)/9 > 5.0655: {rem}")
    failed |= not (ok1 and ok2 and ok3 and ok5 and enc)
    if controls:
        c1 = prop52(F(50655, 10000))
        r1 = not c1[0]
        print(f"   control: k = 5.0655 refused: {r1}" + (f" (fails at R = {float(c1[1]):.6f})" if c1[1] is not None else ""))
        r2 = not (F(21241, 10000) ** 4 < B4lo)
        print(f"   control: the printed bound moved to 2.1241 refused: {r2}")
        r3 = not primitive_identity(F(7))
        print(f"   control: a wrong primitive (7 log t) refused: {r3}")
        failed |= not (r1 and r2 and r3)
    print(f"run time {time.time() - t0:.1f} s")
    print("RESULT:", "FAIL" if failed else "PASS: C_22a >= 531/250 = 2.124, given the paper's Sections 2-4 (not machine-checked here)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
