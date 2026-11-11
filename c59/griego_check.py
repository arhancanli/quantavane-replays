#!/usr/bin/env python3
"""Independent exact check of the 59a row [G2026], K_2 < 0.3174541 (Griego, PR #75), from the construction stated in
the pull request alone (its integers L, A, C, D, N and r; the script in the pull request was neither read nor run).

The construction: Q(z, w) = L - A z - C w + D z w, G(z, w) = D - C z - A w + L z w, F = G/Q, and p = the double Fejer
mean of F of order M = N + 1, i.e. p(z, w) = sum_{0<=j,k<=N} (1 - j/M)(1 - k/M) a_jk z^j w^k with F = sum a_jk z^j w^k.
We check, exactly (integers and rationals):
  1. Q has no zero on the closed bidisc: for |w| <= 1 the root z = (L - C w)/(A - D w) of Q(., w) has modulus > 1,
     because |A - D w| < |L - C w| on |w| = 1 (|L - Cw|^2 - |A - Dw|^2 = L^2 + C^2 - A^2 - D^2 - 2 (LC - AD) Re w,
     smallest at Re w = 1 since LC - AD > 0, where it is (L - C)^2 - (D - A)^2 > 0) and by the maximum principle for
     the analytic (A - D w)/(L - C w) (L > C). So F is analytic on a neighbourhood of the closed bidisc;
     G = z w conj(Q) on the torus (coefficient identity), so |F| = 1 there and |F| <= 1 on the bidisc; the Fejer
     kernel is a probability kernel, so |p| <= 1 on the torus and on the bidisc, and p is not constant, so p is in the
     Schur class S_2 of the page.
  2. The 251 x 251 Taylor coefficients a_jk, as integers nu_jk = a_jk L^(j+k+1), by two routes that share no code:
     (i) the recurrence from Q F = G; (ii) the coefficients of 1/Q from Q (1/Q) = 1, then the product with G. They must
     agree at every (j, k).
  3. ||p_r||_1 = sum_{j,k<=N} (M-j)(M-k)/M^2 |a_jk| r^(j+k) > 1 at r = 3174541/10^7, as an integer comparison over the
     common denominator M^2 L^(2N+1) T^(2N), T = 10^7; and the stated margin 1/7307638490. Then K_2 <= r, and by
     continuity of s -> ||p_s||_1, K_2 < r.
Controls: r = 0.31745 (below the row) gives ||p_r||_1 < 1; D + 1 breaks the strict inequality of step 1 (Q then has a
zero on the torus); a G with A and C swapped is not z w conj(Q) on the torus; one coefficient of route (ii) changed is
caught by the comparison.
  python3 -I griego_check.py [--controls]
Python 3.9+, standard library only."""
import sys
import time
from fractions import Fraction

L, A, C, D = 10 ** 6, 641550, 660190, 981359
N = 250
M = N + 1
R, T = 3174541, 10 ** 7
MARGIN_DEN = 7307638490


def no_zero_on_closed_bidisc(L, A, C, D):
    """step 1: True iff |A - D w| < |L - C w| on |w| = 1 (with LC - AD > 0) and L > C"""
    if not (L > C and L * C - A * D > 0):
        return False, None
    gap = (L - C) ** 2 - (D - A) ** 2
    return gap > 0, gap


def g_is_reflection(L, A, C, D, G):
    """z w conj(Q(z, w)) on the torus has coefficients {zw: L, w: -A, z: -C, 1: D}, all real"""
    want = {(1, 1): L, (0, 1): -A, (1, 0): -C, (0, 0): D}
    return G == want


def nu_route_i(L, A, C, D, G, n):
    """nu_jk = a_jk L^(j+k+1) from Q F = G:  L a_jk = A a_{j-1,k} + C a_{j,k-1} - D a_{j-1,k-1} + g_jk"""
    nu = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for k in range(n + 1):
            v = G.get((j, k), 0) * L ** (j + k)
            if j:
                v += A * nu[j - 1][k]
            if k:
                v += C * nu[j][k - 1]
            if j and k:
                v -= D * L * nu[j - 1][k - 1]
            nu[j][k] = v
    return nu


def nu_route_ii(L, A, C, D, G, n):
    """beta_jk = b_jk L^(j+k+1) for 1/Q = sum b_jk z^j w^k, then a = G * b"""
    beta = [[0] * (n + 1) for _ in range(n + 1)]
    for s in range(2 * n + 1):                     # by anti-diagonals, a different order from route (i)
        for j in range(max(0, s - n), min(n, s) + 1):
            k = s - j
            v = L ** (j + k) if (j, k) == (0, 0) else 0
            if j:
                v += A * beta[j - 1][k]
            if k:
                v += C * beta[j][k - 1]
            if j and k:
                v -= D * L * beta[j - 1][k - 1]
            beta[j][k] = v
    nu = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for k in range(n + 1):
            v = G[(0, 0)] * beta[j][k]
            if j:
                v += G[(1, 0)] * L * beta[j - 1][k]
            if k:
                v += G[(0, 1)] * L * beta[j][k - 1]
            if j and k:
                v += G[(1, 1)] * L * L * beta[j - 1][k - 1]
            nu[j][k] = v
    return nu


def bohr_sum_vs_one(nu, r_num, r_den, n=N, m=M):
    """returns (S_num, S_den) with ||p_r||_1 = S_num / S_den exactly (not reduced)"""
    TL = r_den * L
    top = 2 * n
    s = 0
    pw_r = [r_num ** e for e in range(top + 1)]
    pw_tl = [TL ** e for e in range(top + 1)]
    for j in range(n + 1):
        wj = m - j
        row = nu[j]
        for k in range(n + 1):
            e = j + k
            s += wj * (m - k) * abs(row[k]) * pw_r[e] * pw_tl[top - e]
    den = m * m * L * TL ** top
    return s, den


def dec(num, den, digits=25):
    q = num * 10 ** digits // den
    s = str(q).rjust(digits + 1, "0")
    return s[:-digits] + "." + s[-digits:]


def main():
    controls = "--controls" in sys.argv[1:]
    t0 = time.time()
    ok = True
    G = {(0, 0): D, (1, 0): -C, (0, 1): -A, (1, 1): L}
    nz, gap = no_zero_on_closed_bidisc(L, A, C, D)
    print(f"[G2026] (PR #75): L = {L}, A = {A}, C = {C}, D = {D}, N = {N}, M = {M}, r = {R}/{T}")
    print(f"step 1: LC - AD = {L * C - A * D} > 0, (L - C)^2 - (D - A)^2 = {gap} > 0, L > C: {nz}; so Q has no zero on "
          f"the closed bidisc")
    refl = g_is_reflection(L, A, C, D, G)
    print(f"        G = z w conj(Q) on the torus (coefficients): {refl}; so |F| = 1 there, |F| <= 1 on the bidisc, and "
          f"the Fejer mean p has |p| <= 1 (p is not constant)")
    ok &= nz and refl
    nu1 = nu_route_i(L, A, C, D, G, N)
    nu2 = nu_route_ii(L, A, C, D, G, N)
    same = nu1 == nu2
    big = max(abs(x) for row in nu1 for x in row)
    print(f"step 2: 251 x 251 coefficients nu_jk = a_jk L^(j+k+1) by two routes (Q F = G; 1/Q then G): identical at "
          f"every (j, k): {same}; a_00 = D/L = {Fraction(nu1[0][0], L)}, a_10 = {Fraction(nu1[1][0], L ** 2)}, "
          f"a_01 = {Fraction(nu1[0][1], L ** 2)}; largest |nu| has {len(str(big))} digits")
    ok &= same
    sn, sd = bohr_sum_vs_one(nu1, R, T)
    gt1 = sn > sd
    margin_ok = sn * MARGIN_DEN > sd * (MARGIN_DEN + 1)
    print(f"step 3: ||p_r||_1 at r = 0.3174541 = {dec(sn, sd)}...; > 1: {gt1}; > 1 + 1/{MARGIN_DEN} (the stated "
          f"margin): {margin_ok}; ||p_r||_1 - 1 = {float(Fraction(sn - sd, sd)):.6e}")
    ok &= gt1 and margin_ok
    print(f"so K_2 <= 0.3174541, and K_2 < 0.3174541 by continuity in r ({time.time() - t0:.0f} s)")
    cok = True
    if controls:
        print("\ncontrols (each must come out as stated):")
        sn2, sd2 = bohr_sum_vs_one(nu1, 31745, 10 ** 5)
        c = sn2 < sd2
        print(f"[ctl] r = 0.31745 (below the row): ||p_r||_1 = {dec(sn2, sd2, 15)}... < 1: {c}")
        cok &= c
        nzx, gapx = no_zero_on_closed_bidisc(L, A, C, D + 1)
        c = not nzx
        qz = L + A - C - (D + 1)                    # Q(-1, 1) with D + 1
        c = c and qz == 0
        print(f"[ctl] D + 1: (L - C)^2 - (D + 1 - A)^2 = {gapx}, not > 0, so step 1 is refused; indeed Q(-1, 1) = "
              f"{qz}, a zero on the torus: {c}")
        cok &= c
        Gs = {(0, 0): D, (1, 0): -A, (0, 1): -C, (1, 1): L}
        c = not g_is_reflection(L, A, C, D, Gs)
        print(f"[ctl] G with A and C swapped is not z w conj(Q): {c}")
        cok &= c
        nu3 = [row[:] for row in nu2]
        nu3[137][41] += 1
        c = nu3 != nu1
        print(f"[ctl] one coefficient of route (ii) changed by 1 (j = 137, k = 41) is caught by the comparison: {c}")
        cok &= c
        print(f"controls all as expected: {cok}")
    res = ok and cok
    print(f"\nRESULT: {'PASS' if res else 'FAIL'} (|p| <= 1 on the bidisc and ||p_r||_1 > 1 + 1/{MARGIN_DEN} at "
          f"r = 0.3174541, so K_2 < 0.3174541)")
    return 0 if res else 1


if __name__ == "__main__":
    sys.exit(main())
