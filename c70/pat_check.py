#!/usr/bin/env python3
"""Replay of the 70a row [Pat2026], C_70 >= 1.201927 (cube B_1 = [0, a]^n and ball B_2 = r B with n = 300000,
a/r = 0.013; the 70a page, from #226), with our own code: the full Steiner sum of the row's criterion,

    C_70 >= R := ( sum_{j=0}^{n} C(n,j) a^j kappa_{n-j} r^(n-j) )^(1/n) / ( a + r kappa_n^(1/n) ),
    kappa_m = pi^(m/2) / Gamma(1 + m/2),

enclosed rigorously (with r = 1; R depends on a/r only). It reuses the fixed-point interval primitives of our
c70_check.py (scale 2^-400, outward rounding; logs by atanh series, pi by Machin, log m! by Stirling with the remainder
between 0 and the first omitted term, exp by outward-rounded Taylor sums).

How: with k = n - j the terms are T_k = C(n,k) a^(n-k) kappa_k. Within each parity, T_{k+2}/T_k =
(n-k)(n-k-1)/((k+1)(k+2)) * a^-2 * 2 pi/(k+2), a rational times pi, decreasing in k, so each parity chain is
log-concave. Two anchors (an even k0 and k0 + 1) are evaluated directly (Gamma at a half-integer through
(2m+2)!/(4^(m+1) (m+1)!) sqrt(pi)); every other term in a window of 2W + 2 terms is reached by the exact ratios.
The window sum is a lower bound for the full sum; log-concavity bounds the rest by (n + 1) times the larger edge term.
Not checked: the criterion itself (Szasz's inequality with AM-GM, intrinsic volumes, Steiner's formula), which is
mathematics; our c70 replay of the [Pin2026] row covers the same criterion's single-term form.

Usage: python3 -I pat_check.py C70_CHECK_PY [CANDIDATE_JSON] [--controls]
(CANDIDATE_JSON: pin/certificates/c70_candidate.json as fetch.sh downloads it; its baseline fields record this row's
parameters, which are compared with ours.)
"""
import importlib.util
import sys
from fractions import Fraction

N = 300000
A = Fraction(13, 1000)
CLAIM = Fraction(1201927, 1000000)


def load(path):
    spec = importlib.util.spec_from_file_location("c70_check", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    c = load(sys.argv[1])
    controls = "--controls" in sys.argv
    P, ONE = c.P, c.ONE

    def mul(x, y):          # positive intervals at scale 2^P
        assert x[0] >= 0 and y[0] >= 0
        return (x[0] * y[0]) >> P, c.cdiv(x[1] * y[1], ONE)

    def ln_gamma_half(k):   # log Gamma(1 + k/2)
        if k % 2 == 0:
            return c.lnfact(k // 2) if k >= 2 else (0, 0)
        m = (k - 1) // 2    # Gamma(m + 3/2) = (2m+2)! sqrt(pi) / (4^(m+1) (m+1)!)
        v = c.sub(c.lnfact(2 * m + 2), c.lnfact(m + 1))
        v = c.add(v, c.qmul(c.LNPI, Fraction(1, 2)))
        return c.sub(v, c.qmul(c.LN2, 2 * (m + 1)))

    def log_term(k, a):     # log T_k = log C(n,k) + (n-k) log a + (k/2) log pi - log Gamma(1 + k/2)
        v = c.sub(c.sub(c.lnfact(N), c.lnfact(k)), c.lnfact(N - k))
        v = c.add(v, c.qmul(c.ln_q(a), N - k))
        v = c.add(v, c.qmul(c.LNPI, Fraction(k, 2)))
        return c.sub(v, ln_gamma_half(k))

    def chain(k_start, steps, a, direction):
        """sum over the parity chain from k_start (exclusive) of T_k/T_start for `steps` steps, and the last ratio."""
        tot = (0, 0)
        r = (ONE, ONE)
        k = k_start
        for _ in range(steps):
            if direction > 0:
                q = Fraction((N - k) * (N - k - 1), (k + 1) * (k + 2)) / (a * a) * Fraction(2, k + 2)
                k += 2
            else:
                k -= 2
                q = 1 / (Fraction((N - k) * (N - k - 1), (k + 1) * (k + 2)) / (a * a) * Fraction(2, k + 2))
            if direction > 0:
                r = mul(c.qmul(r, q), c.PI)
            else:
                r = mul(c.qmul(r, q), (c.cdiv(ONE * ONE, c.PI[1]) - 1, c.cdiv(ONE * ONE, c.PI[0]) + 1))
            tot = (tot[0] + r[0], tot[1] + r[1])
        return tot, r

    def R_bounds(a, k0, W):
        """rigorous [lo, hi] for R at side a (r = 1), window of W steps each way in each parity around k0 (even)."""
        lt0 = log_term(k0, a)
        lt1 = log_term(k0 + 1, a)
        up0, e_up0 = chain(k0, W, a, +1)
        dn0, e_dn0 = chain(k0, W, a, -1)
        up1, e_up1 = chain(k0 + 1, W, a, +1)
        dn1, e_dn1 = chain(k0 + 1, W, a, -1)
        even = (ONE + up0[0] + dn0[0], ONE + up0[1] + dn0[1])           # sum / T_k0
        odd = (ONE + up1[0] + dn1[0], ONE + up1[1] + dn1[1])            # sum / T_k0+1
        ratio = c.exp_iv(c.sub(lt1, lt0))                                # T_k0+1 / T_k0
        win = (even[0] + mul(ratio, odd)[0], even[1] + mul(ratio, odd)[1])
        # tails: each parity chain is log-concave, so every omitted term is at most its chain's edge term
        edge = max(e_up0[1], e_dn0[1], mul(ratio, e_up1)[1], mul(ratio, e_dn1)[1])
        tail = (0, (N + 1) * edge)
        lsum_lo = c.add(lt0, c.ln_iv((win[0], win[0])))
        lsum_hi = c.add(lt0, c.ln_iv((win[1] + tail[1], win[1] + tail[1])))
        lkn = c.sub(c.qmul(c.LNPI, Fraction(N, 2)), c.lnfact(N // 2))
        kn = c.exp_iv(c.qmul(lkn, Fraction(1, N)))                       # kappa_n^(1/n)
        den = c.add(c.qiv(a), kn)
        lden = c.ln_iv(den)
        lo = c.exp_iv((c.qmul(lsum_lo, Fraction(1, N))[0] - lden[1], 0))[0]
        hi = c.exp_iv((0, c.qmul(lsum_hi, Fraction(1, N))[1] - lden[0]))[1]
        return lo, hi, edge, win

    print("70a row [Pat2026]: C_70 >= 1.201927 (n = 300000, a/r = 13/1000, the full Steiner sum)")
    status_json = True
    if len(sys.argv) > 2 and not sys.argv[2].startswith("--"):
        import json
        with open(sys.argv[2], encoding="utf-8") as fh:
            js = json.load(fh)
        status_json = (js.get("baseline_n"), Fraction(js.get("baseline_side")), Fraction(js.get("baseline_claim"))) == \
            (N, A, CLAIM)
        print("the parameters equal the baseline fields of pin/certificates/c70_candidate.json (n = %s, side = %s, "
              "claim = %s): %s" % (js.get("baseline_n"), js.get("baseline_side"), js.get("baseline_claim"), status_json))
    k0 = 110220
    lo, hi, edge, win = R_bounds(A, k0, 3000)
    print("window: k in [%d, %d] (both parities), anchors k0 = %d and %d" % (k0 - 6000, k0 + 6001, k0, k0 + 1))
    print("window sum / T_k0 in %s; largest edge term / T_k0 <= %.3e" % (c.enc(win, 12), float(Fraction(edge, ONE))))
    print("R in [%s, %s]" % (c.dec_floor(lo, 15), c.dec_ceil(hi, 15)))
    ok = lo >= c.qiv(CLAIM)[1]
    print("the row's 1.201927 <= R: %s (margin %s)" % (ok, c.dec_floor(lo - c.qiv(CLAIM)[1], 12)))
    # the single largest term alone is not enough at this n, which is why the row keeps the whole sum
    lt0 = log_term(k0, A)
    lkn = c.sub(c.qmul(c.LNPI, Fraction(N, 2)), c.lnfact(N // 2))
    kn = c.exp_iv(c.qmul(lkn, Fraction(1, N)))
    single_hi = c.exp_iv((0, c.qmul(lt0, Fraction(1, N))[1] - c.ln_iv(c.add(c.qiv(A), kn))[0]))[1]
    print("the largest term alone gives at most %s < 1.201927: %s" % (c.dec_ceil(single_hi, 10), single_hi < c.qiv(CLAIM)[0]))
    status = ok and single_hi < c.qiv(CLAIM)[0] and status_json
    if controls:
        print("")
        print("controls (each must come out as stated):")
        ctl = True
        lo2, hi2, _, _ = R_bounds(Fraction(12, 1000), k0, 3000)
        r1 = hi2 < c.qiv(CLAIM)[0]
        print("[ctl] a/r = 0.012: R <= %s < 1.201927: %s" % (c.dec_ceil(hi2, 10), r1))
        ctl &= r1
        r2 = not (lo >= c.qiv(Fraction(1201928, 1000000))[1])
        print("[ctl] the value 1.201928 is not certified: %s" % r2)
        ctl &= r2
        lo3, hi3, _, _ = R_bounds(A, k0, 300)
        r3 = lo3 <= lo and hi3 >= hi
        print("[ctl] a narrower window (300 steps) gives a wider enclosure [%s, %s]: %s"
              % (c.dec_floor(lo3, 12), c.dec_ceil(hi3, 12), r3))
        ctl &= r3
        print("controls all as expected: %s" % ctl)
        status &= ctl
    print("RESULT: %s" % ("PASS (R = 1.2019271425... >= 1.201927 from the full Steiner sum; the single largest term "
                          "alone would not reach it)" if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
