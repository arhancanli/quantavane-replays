"""Independent exact check of K_2 < 0.302825279490988 (C59, the Bohr radius of the bidisc; teorth/optimizationproblems
PR #297), from the certificate's data and the statement in its README alone: t7428_c59.py was neither read nor run.

  ./fetch.sh && python3 -I c59_check.py [CERT] [--controls]

The certificate gives a 3 x 3 matrix U = Ut/den of Gaussian rationals, the truncation N and r = R/E. With
U = [[A, B], [C, D]] (blocks of sizes 1 and 2) and Delta = diag(z, w), the transfer function is
f(z, w) = A + B Delta (I - D Delta)^-1 C. We check, in integers only:
  1. U* U = I and U U* = I exactly. Then f maps the bidisc into the unit disc: with g = (I - D Delta)^-1 C, unitarity
     gives 1 - |f|^2 = g* (I - Delta* Delta) g >= 0 (I - D Delta is invertible since ||D|| <= 1 and ||Delta|| < 1), and f
     is not constant, so |f| < 1 by the maximum principle: f is in the Schur class S_2 of the page.
  2. The coefficients f_jk, 0 <= j, k <= N, from G = C + D Delta G and f = A + B Delta G: the vector coefficient g_jk
     is D (e1 g_{j-1,k} + e2 g_{j,k-1}) and f_jk = B (e1 g_{j-1,k} + e2 g_{j,k-1}) for (j, k) != (0, 0), all exact
     Gaussian integers over den^(j+k+1). They are checked against the closed form of f at two points.
  3. sum_{j,k <= N} |f_jk| r^(j+k) > 1, from below: each |f_jk| is rounded down by an integer square root at 10^-80
     relative to the term. All terms of the full Bohr sum are nonnegative, so ||f_r||_1 > 1, and since r -> ||f_r||_1
     is increasing and continuous, K_2 < r (page definition: K_2 = sup{r : ||f_r||_1 <= 1 for all f in S_2}).
Python 3.9+, standard library only."""
import json
import math
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PREV = Fraction(302825279492, 10 ** 12)      # P2026 (Patel), the page's best upper bound
P = 10 ** 80                                  # fixed-point scale for the square roots


# Gaussian integers as pairs (re, im)
def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def conj(a):
    return (a[0], -a[1])


def load(path):
    d = json.load(open(path))
    den = int(d["den"])
    Ut = [[(int(e[0]), int(e[1])) for e in row] for row in d["U"]]
    R, E = (int(v) for v in d["r"])
    return d, den, Ut, int(d["N"]), Fraction(R, E)


def unitary(Ut, den):
    """U* U = I and U U* = I, exactly: sum_j conj(Ut_ji) Ut_jk = delta_ik den^2, and the same with rows."""
    n = len(Ut)
    for i in range(n):
        for k in range(n):
            s1 = (0, 0)
            s2 = (0, 0)
            for j in range(n):
                s1 = gadd(s1, gmul(conj(Ut[j][i]), Ut[j][k]))
                s2 = gadd(s2, gmul(Ut[i][j], conj(Ut[k][j])))
            want = (den * den if i == k else 0, 0)
            if s1 != want or s2 != want:
                return False
    return True


def coefficients(Ut, N):
    """F[j][k]: the Gaussian-integer numerator of f_jk over den^(j+k+1)."""
    A, B, C, D = Ut[0][0], (Ut[0][1], Ut[0][2]), (Ut[1][0], Ut[2][0]), ((Ut[1][1], Ut[1][2]), (Ut[2][1], Ut[2][2]))
    G = [[None] * (N + 1) for _ in range(N + 1)]
    F = [[None] * (N + 1) for _ in range(N + 1)]
    G[0][0], F[0][0] = C, A
    for L in range(1, 2 * N + 1):
        for j in range(max(0, L - N), min(N, L) + 1):
            k = L - j
            v = (G[j - 1][k][0] if j >= 1 else (0, 0), G[j][k - 1][1] if k >= 1 else (0, 0))
            G[j][k] = (gadd(gmul(D[0][0], v[0]), gmul(D[0][1], v[1])), gadd(gmul(D[1][0], v[0]), gmul(D[1][1], v[1])))
            F[j][k] = gadd(gmul(B[0], v[0]), gmul(B[1], v[1]))
    return F


def closed_form(Ut, den, z, w):
    """f(z, w) = A + B Delta (I - D Delta)^-1 C for Gaussian-rational z, w given as (Fraction, Fraction) pairs."""
    q = lambda e: (Fraction(e[0], den), Fraction(e[1], den))
    fm = lambda a, b: (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    fa = lambda a, b: (a[0] + b[0], a[1] + b[1])
    fs = lambda a, b: (a[0] - b[0], a[1] - b[1])

    def fdiv(a, b):
        n = b[0] * b[0] + b[1] * b[1]
        return ((a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n)
    A, B1, B2 = q(Ut[0][0]), q(Ut[0][1]), q(Ut[0][2])
    C1, C2 = q(Ut[1][0]), q(Ut[2][0])
    D11, D12, D21, D22 = q(Ut[1][1]), q(Ut[1][2]), q(Ut[2][1]), q(Ut[2][2])
    one, zero = (Fraction(1), Fraction(0)), (Fraction(0), Fraction(0))
    # M = I - D Delta = [[1 - D11 z, -D12 w], [-D21 z, 1 - D22 w]]; solve M g = C by Cramer's rule
    m11, m12 = fs(one, fm(D11, z)), fs(zero, fm(D12, w))
    m21, m22 = fs(zero, fm(D21, z)), fs(one, fm(D22, w))
    det = fs(fm(m11, m22), fm(m12, m21))
    g1 = fdiv(fs(fm(C1, m22), fm(m12, C2)), det)
    g2 = fdiv(fs(fm(m11, C2), fm(m21, C1)), det)
    return fa(A, fa(fm(B1, fm(z, g1)), fm(B2, fm(w, g2))))


def series(F, den, N, z, w):
    """sum_{j,k <= N} f_jk z^j w^k in exact Gaussian rationals."""
    fm = lambda a, b: (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    zp, wp = [(Fraction(1), Fraction(0))], [(Fraction(1), Fraction(0))]
    for _ in range(N):
        zp.append(fm(zp[-1], z))
        wp.append(fm(wp[-1], w))
    s = (Fraction(0), Fraction(0))
    for j in range(N + 1):
        for k in range(N + 1):
            c = (Fraction(F[j][k][0], den ** (j + k + 1)), Fraction(F[j][k][1], den ** (j + k + 1)))
            t = fm(c, fm(zp[j], wp[k]))
            s = (s[0] + t[0], s[1] + t[1])
    return s


def bohr_sum(F, den, N, r):
    """(lower, upper) integer bounds on P * sum_{j,k <= N} |f_jk| r^(j+k): each term rounded down (up) by isqrt."""
    R, E = r.numerator, r.denominator
    lo = hi = 0
    for j in range(N + 1):
        for k in range(N + 1):
            L = j + k
            a, b = F[j][k]
            m2 = a * a + b * b                              # |numerator|^2
            num, den_t = R ** L, den ** (L + 1) * E ** L    # term = sqrt(m2) * num / den_t
            # floor and ceiling of P * sqrt(m2) * num / den_t, via isqrt of (P * num / den_t)^2 * m2 scaled exactly
            x2_num, x2_den = m2 * (P * num) ** 2, den_t ** 2
            fl = math.isqrt(x2_num // x2_den)               # <= P * term
            lo += fl
            hi += fl + 1                                     # floor(sqrt(floor(x))) + 1 > sqrt(x)
    return lo, hi


def run(path, controls):
    t0 = time.time()
    d, den, Ut, N, r = load(path)
    ok = True

    def check(what, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  check {what}: {'PASS' if cond else 'FAIL'}{' (' + detail + ')' if detail else ''}")

    print(f"C59 certificate: U = Ut/den, den has {len(str(den))} digits; N = {N}; r = {r.numerator}/{r.denominator}")
    check("block sizes 1, 1, 1 as stated (n = m = 1), a 3 x 3 matrix", (d["n"], d["m"]) == (1, 1) and len(Ut) == 3
          and all(len(row) == 3 for row in Ut))
    check("U* U = I and U U* = I exactly", unitary(Ut, den))
    F = coefficients(Ut, N)
    check("f is not constant", any(F[j][k] != (0, 0) for j in range(N + 1) for k in range(N + 1) if j + k))
    for z, w in (((Fraction(1, 10), Fraction(0)), (Fraction(0), Fraction(1, 5))),
                 ((Fraction(-3, 50), Fraction(1, 25)), (Fraction(1, 20), Fraction(-1, 10)))):
        exact = closed_form(Ut, den, z, w)
        trunc = series(F, den, N, z, w)
        err = abs(complex(float(exact[0] - trunc[0]), float(exact[1] - trunc[1])))
        bound = 4 * max(abs(complex(float(z[0]), float(z[1]))), abs(complex(float(w[0]), float(w[1])))) ** (N + 1)
        check(f"the coefficients sum to the closed form at z = {float(z[0]):+.2f}{float(z[1]):+.2f}i, "
              f"w = {float(w[0]):+.2f}{float(w[1]):+.2f}i", err <= bound, f"difference {err:.1e}, tail bound {bound:.1e}")
    lo, hi = bohr_sum(F, den, N, r)
    print(f"  sum_(j,k <= {N}) |f_jk| r^(j+k) - 1 in [{(lo - P) / P:.4e}, {(hi - P) / P:.4e}]"
          f"  (README: the lower bound minus 1 is 2.2e-42)")
    check(f"the finite Bohr sum at r = {float(r):.15f} exceeds 1, so K_2 < r", lo > P)
    check(f"and r is below the page's {float(PREV):.12f} (P2026)", r < PREV, f"by {float(PREV - r):.3e}")
    if controls:
        print("Red controls (the same check functions, on altered data):")
        r2 = Fraction(30282527949097, 10 ** 14)
        lo2, hi2 = bohr_sum(F, den, N, r2)
        refused = not lo2 > P
        print(f"  R1 r lowered to {float(r2):.14f} (below the search's limit 0.3028252794909793): the sum is "
              f"{'below 1, so nothing is certified, as expected' if refused else 'still above 1 (control failed)'} "
              f"(minus 1 in [{(lo2 - P) / P:.3e}, {(hi2 - P) / P:.3e}])")
        ok &= refused
        bad = [list(row) for row in Ut]
        bad[1][1] = (bad[1][1][0] + 1, bad[1][1][1])
        refused = not unitary(bad, den)
        print(f"  R2 one entry of Ut changed by one unit: {'refused by the unitarity check, as expected' if refused else 'ACCEPTED (control failed)'}")
        ok &= refused
        lo3, _ = bohr_sum(F, den, 20, r)
        refused = not lo3 > P
        print(f"  R3 the sum cut at N = 20: {'below 1, so the certificate needs the higher terms, as expected' if refused else 'already above 1'}"
              f" (minus 1 = {(lo3 - P) / P:.3e})")
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.1f} s)")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = args[0] if args else os.path.join(HERE, "certs", "gomyway", "c59_cert.json")
    sys.exit(0 if run(path, "--controls" in sys.argv) else 1)


if __name__ == "__main__":
    main()
