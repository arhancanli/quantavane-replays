#!/usr/bin/env python3
"""Independent exact check of the 59a row [BPWW2026], K_2 < 0.3177 (Baran, Pikul, Woerdeman, Wojtylak,
arXiv:2504.03236, Theorem 6.4), from the construction in the paper alone.

The paper: f = p~/p with p(z) = det(I_2 - D diag(z1, z2)) = 1 + a01 z1 + a10 z2 + a00 z1 z2 (a01 = -D11, a10 = -D22,
a00 = det D) and p~(z) = z1 z2 p(1/z1, 1/z2) = z1 z2 + a01 z2 + a10 z1 + a00, for the matrix D with the paper's
decimals (read as exact rationals); then sum_{a1, a2 <= 12} |c_a| rho0^(a1 + a2) > 1 at rho0 = 0.3177.
We check, exactly:
  1. ||D|| < 1 (I - D^T D positive definite: its (1,1) entry and determinant are positive), so p has no zero on the
     closed bidisc (|det(I - D Z)| > 0 for ||Z|| <= 1) and f is a rational inner function: |f| = 1 on the torus,
     so |f| < 1 on the open bidisc (f is not constant), i.e. f is in the Schur class of the page;
  2. the coefficients c_jk, 0 <= j, k <= 12, from p f = p~ (exact rationals; the coefficients are real), and the
     identity p * (series) = p~ up to degree (12, 12);
  3. sum |c_jk| rho0^(j+k) > 1 exactly, at rho0 = 3177/10000. All terms of the full Bohr sum are nonnegative, so
     ||f_rho0||_1 > 1 and K_2 < rho0 (by continuity also slightly below it).
  python3 -I bpww_check.py [--controls]
Python 3.9+, standard library only."""
import sys
import time
from fractions import Fraction as F

D = [[F("0.854373111798292"), F("-0.518782521594128")], [F("0.518794363700548"), F("0.848187547437653")]]
RHO = F(3177, 10000)
N = 12


def coeffs(Dm, n):
    a01, a10 = -Dm[0][0], -Dm[1][1]
    a00 = Dm[0][0] * Dm[1][1] - Dm[0][1] * Dm[1][0]
    p = {(0, 0): F(1), (1, 0): a01, (0, 1): a10, (1, 1): a00}
    pt = {(1, 1): F(1), (0, 1): a01, (1, 0): a10, (0, 0): a00}
    c = {}
    for j in range(n + 1):
        for k in range(n + 1):
            v = pt.get((j, k), F(0))
            if j:
                v -= p[(1, 0)] * c[j - 1, k]
            if k:
                v -= p[(0, 1)] * c[j, k - 1]
            if j and k:
                v -= p[(1, 1)] * c[j - 1, k - 1]
            c[j, k] = v                                  # p(0, 0) = 1
    return c, p, pt


def bohr(c, rho, n):
    return sum(abs(v) * rho ** (j + k) for (j, k), v in c.items() if j <= n and k <= n)


def norm_below_one(Dm):
    # M = I - D^T D
    m11 = 1 - (Dm[0][0] ** 2 + Dm[1][0] ** 2)
    m22 = 1 - (Dm[0][1] ** 2 + Dm[1][1] ** 2)
    m12 = -(Dm[0][0] * Dm[0][1] + Dm[1][0] * Dm[1][1])
    return m11 > 0 and m11 * m22 - m12 * m12 > 0, m11, m11 * m22 - m12 * m12


def main():
    t0 = time.time()
    ok = True
    print("59a row [BPWW2026]: K_2 < 0.3177 (Theorem 6.4; f = p~/p, p = det(I - D diag(z)), coefficients up to 12)")
    good, m11, det = norm_below_one(D)
    print(f"  ||D|| < 1: I - D^T D has (1,1) entry {float(m11):.6f} > 0 and determinant {float(det):.3e} > 0: {good}")
    ok &= good
    c, p, pt = coeffs(D, N)
    # the identity p * series = p~ for all monomials of bidegree <= (N, N)
    resid = 0
    for j in range(N + 1):
        for k in range(N + 1):
            s = sum(pv * c.get((j - a, k - b), 0) for (a, b), pv in p.items() if j - a >= 0 and k - b >= 0)
            resid += s != pt.get((j, k), 0)
    print(f"  {len(c)} coefficients c_jk (0 <= j, k <= 12), exact rationals; p * series = p~ through bidegree (12, 12): "
          f"{resid == 0}")
    ok &= resid == 0
    S = bohr(c, RHO, N)
    good = S > 1
    print(f"  sum |c_jk| 0.3177^(j+k) = 1 + {float(S - 1):.10e} > 1 (exact): {good}")
    ok &= good
    if "--controls" in sys.argv:
        print("Red controls (each must leave the finite sum at or below 1):")
        for what, (Dm, rho) in (("rho0 = 0.3170", (D, F(3170, 10000))),
                                ("rho0 = 0.3006, Knese's lower bound", (D, F(3006, 10000))),
                                ("D with its off-diagonal entries set to 0",
                                 ([[D[0][0], F(0)], [F(0), D[1][1]]], RHO))):
            cc, _, _ = coeffs(Dm, N)
            Sc = bohr(cc, rho, N)
            good = Sc <= 1
            print(f"  {what}: the finite sum is 1 {'+' if Sc > 1 else '-'} {abs(float(Sc - 1)):.3e}: "
                  f"{'not above 1, as expected' if good else 'ABOVE 1 (control failed)'}")
            ok &= good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
