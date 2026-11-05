"""Independent check of the two bounds that pin C16 (Brezis-Gallouet-Wainger remainder constant on T^2).

Written from the definition on the constant page (teorth/optimizationproblems, constants/16a.md): with T^2 = [0, 2pi]^2,
u = sum_k c_k e^{ik.x} zero-mean, |grad u|^2 = 4 pi^2 sum |k|^2 |c_k|^2, |Lap u|^2 = 4 pi^2 sum |k|^4 |c_k|^2,
delta = |Lap u|^2 / |grad u|^2 >= 1, and L = max_{delta >= 1} [4 pi Theta(delta) - ln delta - ln(1 + ln delta)].

Upper bound (certificate: intervals [d0, d1] with one mu each). Cauchy-Schwarz with weights |k|^2 (1 + mu |k|^2) gives
4 pi |u(0)|^2 / |grad u|^2 <= G(mu) (1 + mu delta) / pi, G(mu) = sum_{k != 0} 1 / (|k|^2 (1 + mu |k|^2)), and by translation
|u(0)| may stand for the sup norm. U(delta) = G (1 + mu delta)/pi - ln delta - ln(1 + ln delta) is convex in delta (a line
plus two convex terms), so on each interval it is at most its larger endpoint value. For delta >= 2^18 we use our own
lattice-count bound (see tail_bound) with mu = 1/(delta ln delta).

Lower bound: the untruncated trial function c_k = 1/(|k|^2 (1 + mu |k|^2)), for which u(0) = S0 = G(mu) is the sup norm
(all c_k > 0), S1 = sum |k|^2 c_k^2 = G - H/mu, S2 = sum |k|^4 c_k^2 = H/mu^2 with H = sum_{k != 0} (|k|^2 + 1/mu)^-2, so
L >= S0^2 / (pi S1) - ln(S2/S1) - ln(1 + ln(S2/S1)).

Lattice sums by rows: sum_y 1/(y^2 + s^2) = (pi/s) coth(pi s), sum_y 1/(y^2 + s^2)^2 = (pi/(2 s^3)) coth(pi s)
+ (pi^2/(2 s^2)) csch(pi s)^2. Rows x >= N are bracketed with convexity (trapezoid below, midpoint above) using closed-form
integrals. Arithmetic: decimal at 50 digits; every enclosure is widened by SLACK = 1e-30, far above the accumulated
rounding (each operation is correctly rounded, at most 1e6 operations per sum, all values below 1e3: under 1e-40).
Standard library only.  python3 -I c16_check.py certificate.json [--lower MU_HEX_OR_DECIMAL] [--controls]
"""
import json, math, sys, time
from decimal import Decimal as D, getcontext

getcontext().prec = 50
SLACK = D("1e-30")
X_EXACT = 40          # rows x < X_EXACT use coth exactly; beyond, 0 < coth(pi s) - 1 < 3 e^{-2 pi s} <= 3 e^{-80 pi} < 1e-108
EXP_TAIL = D("1e-100")  # bounds the dropped coth - 1 and csch^2 corrections of all rows x >= X_EXACT together (each sum < 1e-107)


def pi_dec():
    """pi by Machin's formula at 60 digits (error far below 1e-55)."""
    getcontext().prec = 70
    def arctan_inv(n):
        x = D(1) / n; x2 = x * x; term = x; s = x; k = 1
        while True:
            term *= -x2; k += 2
            t = term / k
            if abs(t) < D("1e-68"):
                break
            s += t
        return s
    p = 16 * arctan_inv(5) - 4 * arctan_inv(239)
    getcontext().prec = 50
    return +p


PI = pi_dec()
assert abs(float(PI) - math.pi) < 1e-15


def coth(z):
    e = (-2 * z).exp()
    return (1 + e) / (1 - e)


def csch2(z):
    e = (-2 * z).exp()            # csch^2 z = 4 e^{-2z} / (1 - e^{-2z})^2
    return 4 * e / (1 - e) ** 2


def row_sum1(s):
    """sum over all integers y of 1/(y^2 + s^2), s > 0."""
    return PI / s * coth(PI * s)


def row_sum2(s):
    """sum over all integers y of 1/(y^2 + s^2)^2, s > 0."""
    return PI / (2 * s ** 3) * coth(PI * s) + PI ** 2 / (2 * s ** 2) * csch2(PI * s)


def n_rows(a, n_min):
    """number of explicit rows: at least n_min and 8 sqrt(a), as a power of two."""
    n = n_min
    while n < 8 * math.sqrt(float(a)) + 64:
        n *= 2
    return n


def G_bracket(mu, n):
    """[lo, hi] enclosing G(mu) = sum_{k != 0} [1/r - 1/(r + a)], a = 1/mu, with n explicit rows (x = 1..n-1)."""
    a = 1 / mu
    sa = a.sqrt()
    total = PI ** 2 / 3 - (row_sum1(sa) - 1 / a)          # row x = 0: sum_{y != 0} [1/y^2 - 1/(y^2 + a)]
    rows = D(0)
    for x in range(1, min(X_EXACT, n)):
        X = D(x)
        rows += row_sum1(X) - row_sum1((X * X + a).sqrt())
    for x in range(X_EXACT, n):
        X = D(x)
        rows += PI * (1 / X - 1 / (X * X + a).sqrt())   # coth corrections bounded below by EXP_TAIL
    N = D(n)
    def I(M):                                           # integral_M^inf (1/x - 1/sqrt(x^2 + a)) dx
        return ((M + (M * M + a).sqrt()) / (2 * M)).ln()
    phiN = 1 / N - 1 / (N * N + a).sqrt()
    tail_lo = PI * (I(N) + phiN / 2)                    # trapezoid: phi convex, so sum_{x >= N} phi >= I(N) + phi(N)/2
    tail_hi = PI * I(N - D("0.5"))                      # midpoint: sum_{x >= N} phi <= I(N - 1/2)
    lo = total + 2 * (rows + tail_lo) - 2 * EXP_TAIL - SLACK
    hi = total + 2 * (rows + tail_hi) + 2 * EXP_TAIL + SLACK
    return lo, hi


def H_bracket(mu, n):
    """[lo, hi] enclosing H = sum_{k != 0} 1/(r + a)^2, a = 1/mu."""
    a = 1 / mu
    total = row_sum2(a.sqrt()) - 1 / (a * a)            # row x = 0 without y = 0
    rows = D(0)
    for x in range(1, min(X_EXACT, n)):
        X = D(x)
        rows += row_sum2((X * X + a).sqrt())
    for x in range(X_EXACT, n):
        s = (D(x) ** 2 + a).sqrt()
        rows += PI / (2 * s ** 3)                       # coth -> 1 and csch^2 -> 0 beyond X_EXACT, bounded by EXP_TAIL
    N = D(n)
    # psi(x) = pi/(2 (x^2 + a)^{3/2}) is decreasing and convex for x > 0; integral of (x^2 + a)^{-3/2} from M to inf
    # is (1 - M/sqrt(M^2 + a)) / a.
    def I(M):
        return PI / 2 * (1 - M / (M * M + a).sqrt()) / a
    psiN = PI / (2 * (N * N + a) ** D("1.5"))
    tail_lo = I(N) + psiN / 2
    tail_hi = I(N - D("0.5"))
    return total + 2 * (rows + tail_lo) - SLACK - EXP_TAIL, total + 2 * (rows + tail_hi) + SLACK + EXP_TAIL


def hexdec(s):
    return D(float.fromhex(s)) if isinstance(s, str) and s.startswith(("0x", "-0x")) else D(s)


def U_hi(G_hi, mu, delta):
    """upper bound of G (1 + mu delta)/pi - ln delta - ln(1 + ln delta) (ln correctly rounded, then widened)."""
    ld = delta.ln()
    return G_hi * (1 + mu * delta) / PI - ld - (1 + ld).ln() + SLACK


def tail_bound(mu0_exp=-20, n=1 << 14):
    """For delta >= 2^18 with mu = 1/(delta ln delta) <= mu0 = 2^mu0_exp:
    G(mu) <= G(mu0) + pi ln(mu0/mu) + (pi/2) mu0 + (pi^2/sqrt 2) sqrt(mu0)   (lattice count #{0 < |k| <= rho} <= pi (rho + 1/sqrt 2)^2)
    so U <= A (1 + 1/t) + 1 - ln(1 + 1/t) + ln(t)/t with t = ln delta and A = G(mu0)/pi + ln mu0 + mu0/2 + (pi/sqrt 2) sqrt(mu0),
    decreasing in t once A >= 0 and ln t > 2; its value at t0 = ln 2^18 bounds the whole tail."""
    mu0 = D(2) ** mu0_exp
    t0 = (D(2) ** 18).ln()
    assert 1 / ((D(2) ** 18) * t0) <= mu0
    _, Ghi = G_bracket(mu0, n)
    A = Ghi / PI + mu0.ln() + mu0 / 2 + PI / D(2).sqrt() * mu0.sqrt() + SLACK
    assert A >= 0 and t0.ln() > 2
    return A * (1 + 1 / t0) + 1 - (1 + 1 / t0).ln() + t0.ln() / t0 + SLACK, A


def upper(cert, n_min=2048, controls=False):
    T = D(cert["T"])
    iv = [[hexdec(x) for x in row] for row in cert["intervals"]]
    assert iv[0][0] == 1 and iv[-1][1] == D(2) ** 18, "intervals must run from 1 to 2^18"
    for (a0, a1, m), (b0, b1, _) in zip(iv, iv[1:]):
        assert a1 == b0, "gap or overlap between intervals"
    worst = None
    for i, (d0, d1, mu) in enumerate(iv):
        assert d0 < d1 and mu > 0
        n = n_rows(1 / mu, n_min)
        # refine until the interval is clearly below T, or G is enclosed to 1e-19 (so the endpoint values near the
        # peak are reproduced to about 1e-19, not just shown to be below T), or 2^17 rows
        while True:
            lo, hi = G_bracket(mu, n)
            u = max(U_hi(hi, mu, d0), U_hi(hi, mu, d1))
            if u <= T - D("1e-14") or hi - lo < D("1e-19") or n >= 1 << 17:
                break
            n *= 2
        if worst is None or u > worst[0]:
            worst = (u, i, n, hi - lo)
        if u > T:
            return False, f"interval {i} [{float(d0)}, {float(d1)}] mu {float(mu)}: bound {u} > T", worst
    tb, A = tail_bound()
    if tb > T:
        return False, f"tail bound {tb} > T", worst
    return True, f"tail delta >= 2^18: U <= {tb:.12f} (A = {A:.12f})", worst


def lower(mu, n_min=2048):
    n = n_rows(1 / mu, n_min) * 16
    Glo, Ghi = G_bracket(mu, n)
    Hlo, Hhi = H_bracket(mu, n)
    # S0 = G, S1 = G - H/mu, S2 = H/mu^2; bound B = S0^2/(pi S1) - ln(S2/S1) - ln(1 + ln(S2/S1)), lower end:
    S1lo, S1hi = Glo - Hhi / mu, Ghi - Hlo / mu
    S2lo, S2hi = Hlo / mu ** 2, Hhi / mu ** 2
    dlo, dhi = S2lo / S1hi, S2hi / S1lo
    B_lo = Glo ** 2 / (PI * S1hi) - dhi.ln() - (1 + dhi.ln()).ln() - SLACK
    B_hi = Ghi ** 2 / (PI * S1lo) - dlo.ln() - (1 + dlo.ln()).ln() + SLACK
    return B_lo, B_hi, (dlo, dhi), n


def main():
    t0 = time.time()
    cert = json.load(open(sys.argv[1]))
    ok, msg, worst = upper(cert)
    u, i, n, w = worst
    print(f"RESULT upper: C16 <= {cert['T']}: {'CERTIFIED' if ok else 'FAILED'}")
    print(f"  {len(cert['intervals'])} intervals cover [1, 2^18] without gaps; largest endpoint bound {u:.20f}")
    print(f"  at interval {i} (rows {n}, G bracket width {w:.2e}); margin below T {D(cert['T']) - u:.3e}")
    print("  " + msg)
    if "--lower" in sys.argv:
        mu = hexdec(sys.argv[sys.argv.index("--lower") + 1])
        lo, hi, (dlo, dhi), n = lower(mu)
        print(f"RESULT lower: mu = {mu}: C16 >= {lo:.20f} (enclosure width {hi - lo:.2e}; delta in [{dlo:.15f}, {dhi:.15f}], rows {n})")
    if "--controls" in sys.argv:
        # red controls: a target just below our largest endpoint value must fail; one mu perturbed must fail;
        # an interval removed (a gap) must be refused.
        bad = dict(cert); bad["T"] = str(u - D("1e-17"))
        print("control: T lowered to 1e-17 below the largest endpoint value rejected:", not upper(bad)[0])
        bad = dict(cert); iv = [list(r) for r in cert["intervals"]]; m = float.fromhex(iv[i][2]); iv[i][2] = (m * 1.01).hex()
        bad["intervals"] = iv
        print(f"control: mu of interval {i} raised by 1% rejected:", not upper(bad)[0])
        bad = dict(cert); bad["intervals"] = cert["intervals"][:i] + cert["intervals"][i + 1:]
        try:
            r = upper(bad)[0]
        except AssertionError:
            r = False
        print(f"control: interval {i} removed (a gap) rejected:", not r)
    print(f"run time {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
