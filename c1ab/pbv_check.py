#!/usr/bin/env python3
"""Independent replay of the numerical part of the 1a row [PBV2026], C_1a >= 1292/1000 (Piterbarg, Bajaj, Vincent,
lower_bound_proof.tex: Theorem 1.1, the five anchors of its section 4 and Proposition 5.1), from the paper and the multiplier
coefficients alone. The coefficients a_j are read as data from the certificate JSON
(delsarte_dual/grid_bound_alt_kernel/certificates/multiscale_arcsine_1292.json); no program of the repository is read or
run.

The inequality used (the z_1-free master inequality of the paper's section 2, after Matolcsi-Vinuesa; taken as given): for every
admissible f, with M = ||f*f||_inf (int f = 1), Phi(M) := M + 1 + sqrt((M - 1)(K_2 - 1)) >= 2/u + a, where
K = sum_i lambda_i K_arc(delta_i) has Fourier transform Khat(xi) = sum_i lambda_i J_0(pi delta_i xi)^2 (Mehler-Sonine),
K_2 = int Khat^2, u = 1/2 + delta_1, Ktilde(j) = Khat(j/u), G(x) = sum_{j<=200} a_j cos(2 pi j x/u), m_G = min_[0,1/4] G,
S_1 = sum a_j^2/Ktilde(j), a = (4/u) m_G^2/S_1. Phi increases in M and in K_2, so Phi(M0) < 2/u + a at M0 = 1292/1000,
with K_2 from above and a from below, gives M > M0 for every f.

What we compute (decimal at 50 digits; every error bounded):
  J_0 by its power series for x < 20 (alternating; remainder below the first omitted term) and by Hankel's expansion
  for x >= 20, with the remainders of P and Q bounded by their first omitted terms (DLMF 10.17(iii), nu = 0, x > 0);
  K_2 <= 2 [Simpson on [0, T] from upper values, plus T h^4/180 * 16 (2 pi delta_1)^4, the bound for the fourth
  derivative of Khat^2 from |J_0^(k)| <= 1] + 8 C^2/(pi^4 T), C = sum lambda_i/delta_i, from J_0(x)^2 <= 2/(pi x)
  (x (J_0^2 + Y_0^2) increases to 2/pi; Watson 13.74); S_1 from lower bounds of Ktilde(j); m_G on adaptive cells,
  G(x) >= G(c) - |G'(c)| r - B_2 r^2/2 with B_2 = sum |a_j| (2 pi j/u)^2.
  python3 -I pbv_check.py CERT_JSON [--controls]
Python 3.9+, standard library only."""
import json
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction as F

getcontext().prec = 50
EPS = D("1e-35")


def _pi():
    def atn(n):
        x = D(1) / n
        x2, t, s, k = x * x, x, x, 1
        while True:
            t *= -x2
            k += 2
            if abs(t) < D("1e-55"):
                return s
            s += t / k
    return 16 * atn(5) - 4 * atn(239)


PI = _pi()
TWO_PI = 2 * PI


def dec(fr):
    return D(fr.numerator) / D(fr.denominator)


def cos_sin(x):
    k = (x / TWO_PI).to_integral_value()
    y = x - TWO_PI * k
    y2 = y * y
    c, s, tc, ts, n = D(1), y, D(1), y, 0
    while True:
        n += 2
        tc = -tc * y2 / (n * (n - 1))
        ts = -ts * y2 / ((n + 1) * n)
        c += tc
        s += ts
        if abs(tc) < D("1e-48") and abs(ts) < D("1e-48"):
            return c, s


# Hankel coefficients a_k(0) = (-1)^k ((2k-1)!!)^2 / (k! 8^k)
AK = [D(1)]
for _k in range(1, 40):
    AK.append(AK[-1] * (-(2 * _k - 1) ** 2) / (_k * 8))


def J0(x):
    """[lo, hi] for J_0(x), x a positive Decimal"""
    if x < 20:
        t = D(1)
        s = D(1)
        k = 0
        h2 = (x / 2) ** 2
        while True:
            k += 1
            t = -t * h2 / (k * k)
            s += t
            if k > x / 2 + 1 and abs(t) < D("1e-45"):
                nxt = abs(t) * h2 / ((k + 1) ** 2)
                return s - nxt - EPS, s + nxt + EPS
    P = D(0)
    Q = D(0)
    K = 12
    for k in range(K):
        P += (-1) ** k * AK[2 * k] / x ** (2 * k)
        Q += (-1) ** k * AK[2 * k + 1] / x ** (2 * k + 1)
    rp = abs(AK[2 * K] / x ** (2 * K))
    rq = abs(AK[2 * K + 1] / x ** (2 * K + 1))
    c, s = cos_sin(x - PI / 4)
    amp = (2 / (PI * x)).sqrt()
    v = amp * (P * c - Q * s)
    err = amp * (rp + rq) + EPS
    return v - err, v + err


def sq_hi(iv):
    lo, hi = iv
    return max(lo * lo, hi * hi)


def sq_lo(iv):
    lo, hi = iv
    if lo <= 0 <= hi:
        return D(0)
    return min(lo * lo, hi * hi)


def khat(xi, deltas, lambdas, hi=True):
    tot = D(0)
    for dl, lm in zip(deltas, lambdas):
        x = PI * dec(dl) * xi
        v = J0(x) if x > 0 else (D(1), D(1))
        tot += dec(lm) * (sq_hi(v) if hi else sq_lo(v))
    return tot


def K2_upper(deltas, lambdas, T, h):
    """2 [Simpson on [0, T] from upper values + its error bound] + the tail bound"""
    n = int(T / h)
    assert n % 2 == 0 and F(n) * h == T
    hd = dec(h)
    s = D(0)
    for k in range(n + 1):
        xi = hd * k
        g = khat(xi, deltas, lambdas, hi=True)
        w = 1 if k in (0, n) else (4 if k % 2 else 2)
        s += w * g * g
    simpson = s * hd / 3
    c1 = PI * dec(max(deltas))
    f4 = 16 * (2 * c1) ** 4
    simp_err = dec(T) * hd ** 4 / 180 * f4
    C = sum(dec(lm / dl) for dl, lm in zip(deltas, lambdas))
    tail = 8 * C * C / (PI ** 4 * dec(T))
    return 2 * (simpson + simp_err) + tail, (simpson, simp_err, tail)


def S1_upper(a, deltas, lambdas, u):
    tot = D(0)
    for j, aj in enumerate(a, start=1):
        kt = D(0)
        for dl, lm in zip(deltas, lambdas):
            kt += dec(lm) * sq_lo(J0(PI * dec(dl) * j / dec(u)))
        assert kt > 0
        tot += dec(aj) ** 2 / kt
    return tot * (1 + D("1e-40"))


def G_and_dG(a, x, u):
    """G(x) and G'(x) by the Chebyshev recurrences for cos(j t), sin(j t), t = 2 pi x/u"""
    t = TWO_PI * x / dec(u)
    c1, s1 = cos_sin(t)
    cp, sp = D(1), D(0)
    cj, sj = c1, s1
    G = D(0)
    dG = D(0)
    for j, aj in enumerate(a, start=1):
        ajd = dec(aj)
        G += ajd * cj
        dG -= ajd * sj * TWO_PI * j / dec(u)
        cj, cp = 2 * c1 * cj - cp, cj
        sj, sp = 2 * c1 * sj - sp, sj
    return G, dG


def mG_lower(a, u, target):
    """min over [0, 1/4] of G from below, on cells refined until each lower bound exceeds the target (or 2^-24 wide)"""
    B2 = sum(dec(abs(aj)) * (TWO_PI * j / dec(u)) ** 2 for j, aj in enumerate(a, start=1))
    cells = [(F(k, 4096), F(k + 1, 4096)) for k in range(1024)]
    worst = None
    count = 0
    slack = D("1e-30")
    while cells:
        lo, hi = cells.pop()
        count += 1
        c = (lo + hi) / 2
        r = dec((hi - lo) / 2)
        G, dG = G_and_dG(a, dec(c), u)
        bound = G - abs(dG) * r - B2 * r * r / 2 - slack
        if bound >= target or hi - lo <= F(1, 2 ** 24):
            worst = bound if worst is None or bound < worst else worst
        else:
            cells += [(lo, c), (c, hi)]
    return worst, count, B2


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    t0 = time.time()
    body = json.load(open(args[0], encoding="utf-8"))["body"]
    deltas = [F(x) for x in body["kernel"]["deltas_q"]]
    lambdas = [F(x) for x in body["kernel"]["lambdas_q"]]
    a = [F(x) for x in body["G"]["coeffs_q"]]
    u = F(1, 2) + deltas[0]
    ok = True
    print("1a row [PBV2026]: C_1a >= 1292/1000 (the master inequality with the three-scale arcsine kernel)")
    good = deltas == [F(138, 1000), F(55, 1000), F(25, 1000)] and lambdas == [F(85, 100), F(10, 100), F(5, 100)] and \
        sum(lambdas) == 1 and all(x > 0 for x in lambdas) and len(a) == 200 and u == F(319, 500) == F(u)
    print(f"  kernel as in the paper (deltas 138, 55, 25 /1000, weights 85, 10, 5 /100, sum 1), u = 638/1000, "
          f"200 coefficients with denominators dividing 10^8: {good and all(10 ** 8 % x.denominator == 0 for x in a)}")
    ok &= good
    t1 = time.time()
    K2, parts = K2_upper(deltas, lambdas, F(10000), F(1, 80))
    print(f"  K_2 <= {K2:.10f}: Simpson on [0, 10^4] with h = 1/80 gives {2 * parts[0]:.10f}, its error bound "
          f"{2 * parts[1]:.3e}, the tail beyond 10^4 {parts[2]:.3e} ({time.time() - t1:.0f} s)")
    good = K2 <= D(47897) / 10000
    print(f"  the paper's anchor K_2 <= 47897/10000: {good}")
    ok &= good
    S1 = S1_upper(a, deltas, lambdas, u)
    good = S1 <= D(29841) / 1000
    print(f"  S_1 <= {S1:.12f}; the paper's anchor S_1 <= 29841/1000: {good}")
    ok &= good
    t1 = time.time()
    mG, cells, B2 = mG_lower(a, u, D("0.99995"))
    good = mG >= D(998) / 1000
    print(f"  m_G >= {mG:.10f} ({cells} cells, |G''| <= {B2:.1f}); the paper's anchor m_G >= 998/1000: {good} "
          f"({time.time() - t1:.0f} s)")
    ok &= good
    ud = dec(u)
    aa = 4 / ud * mG * mG / S1
    good = aa >= D(20925) / 100000
    print(f"  a = (4/u) m_G^2 / S_1 >= {aa:.10f}; the paper's anchor a >= 20925/100000: {good}")
    ok &= good
    M0 = D(1292) / 1000
    phi = M0 + 1 + ((M0 - 1) * (K2 - 1)).sqrt() * (1 + D("1e-40"))
    tau = 2 / ud + aa
    good = phi < tau
    print(f"  Phi(1.292) <= {phi:.10f} < 2/u + a >= {tau:.10f}: {good} (margin {float(tau - phi):.3e}; "
          f"with the paper's rounded anchors the margin is 9.6e-5)")
    ok &= good
    if "--controls" in sys.argv:
        print("Red controls:")
        M1 = D(12925) / 10000
        phi1 = M1 + 1 + ((M1 - 1) * (K2 - 1)).sqrt()
        good = not (phi1 < tau)
        print(f"  M0 = 1.2925 (above the certificate's 1.29232): Phi = {phi1:.6f} >= 2/u + a: "
              f"{'not certified, as expected' if good else 'CERTIFIED (bad)'}")
        ok &= good
        b = list(a)
        b[0] = -b[0]
        low = min((G_and_dG(b, dec(F(k, 4096)), u)[0], k) for k in range(1025))
        good = low[0] + D("1e-30") < D(998) / 1000
        print(f"  the first coefficient's sign flipped: G({low[1]}/4096) = {low[0]:.6f} < 0.998, so m_G >= 0.998 fails: "
              f"{'as expected' if good else 'NOT BELOW (bad)'}")
        ok &= good
        K2t, _ = K2_upper(deltas, lambdas, F(200), F(1, 20))
        good = K2t > D(47897) / 10000
        print(f"  the same bound with T = 200 (tail 8C^2/(pi^4 T) = {float(8 * sum(dec(l / d) for d, l in zip(deltas, lambdas)) ** 2 / (PI ** 4 * 200)):.3e}): "
              f"K_2 <= {K2t:.6f} only, which does not reach 4.7897: {'as expected' if good else 'REACHED (bad)'}")
        ok &= good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
