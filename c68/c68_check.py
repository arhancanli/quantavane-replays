"""Independent check of the quantitative hypotheses behind C68 >= 54/125 (Korenblum's constant; the 68a row [Pin2026] of Pinho, added in PR #239).

What is checked, from the certificate's exact rationals (certificate-jung.json), with our own code:
Wikstrom's moment lemma (arXiv:2607.17748, Lemma 3.1, as recalled in the certificate's NOTE.md): c is admissible when a
probability measure mu on [c^2, 1] with no atom at 1 satisfies
    (I)  int t^k dmu >= 1/(k+1)            for k >= 0,
    (II) int beta(t) t^k dmu <= 1/(k+1)     for k >= 1,     beta = 1 - g,
with g(t) = 1/G the envelope of GEOMETRIC_LEMMA.md: w = -log c, A = pi^2/w, X = sinh(A/2)^2,
y = sin(pi log t/(2w))^2 / X, H = 3y(2 + y), g = 2 sqrt(H)/(sqrt(H) + sqrt(H + 4)) (increasing in y).
mu = sum_j w_j delta_{t_j} + 1_[R,1] dt, so (I) is A_k = sum w_j t_j^k - R^(k+1)/(k+1) >= 0 and (II) is
B_k = int_R^1 g t^k dt + R^(k+1)/(k+1) - sum w_j beta_j t_j^k >= 0.
On [R, 1], g(t) >= l (1 - t) for the certificate's slope l: we re-derive the admissible range of l (y >= a (1-t)^2
from -log t >= 1 - t and sin v >= v - v^3/6, then g monotone in y) and check l below it. So int_R^1 g t^k dt >=
l int_R^1 (1 - t) t^k dt = l [1/((k+1)(k+2)) - R^(k+1) (1/(k+1) - R/(k+2))].
Tails: (I) for k + 1 >= R / w_last the atom at R alone gives A_k >= R^k (w_last - R/(k+1)) >= 0; (II) for k >= K2,
int_{t0}^{t1} g t^k >= (t1 - t0) l (1 - t1) t0^k exceeds sum_j w_j t_j^k (every t_j < t0, so the ratio only falls with k).

Arithmetic: transcendental constants in decimal at 110 digits (pi by Machin, ln/exp/sqrt correctly rounded, sin by its
alternating series with remainder below 1e-100), each widened by 1e-80 (far above the rounding of a few hundred
operations); moments in 320-bit fixed point with outward rounding at every step. Standard library only.
Usage: python3 -I c68_check.py certificate-jung.json [--controls]
"""
import json, sys, time
from decimal import Decimal as D, getcontext
from fractions import Fraction as F

getcontext().prec = 110
WIDEN = D("1e-80")
P = 320
ONE = 1 << P


def pi_dec():
    def arctan_inv(n):
        x = D(1) / n; x2 = x * x; term = x; s = x; k = 1
        while True:
            term *= -x2; k += 2
            t = term / k
            if abs(t) < D("1e-108"):
                return s
            s += t
    return 16 * arctan_inv(5) - 4 * arctan_inv(239)


PI = pi_dec()


def dec(fr):
    return D(fr.numerator) / D(fr.denominator)


def sin_dec(v):
    """sin v for 0 <= v <= 4 by its alternating series; the remainder is below the last term used (< 1e-100)."""
    assert 0 <= v <= 4
    term = v; s = v; n = 1
    while abs(term) > D("1e-100") or n < 6:
        term = -term * v * v / ((n + 1) * (n + 2)); n += 2
        s += term
    return s


def g_of(t, w, X):
    """g(t) = 1/G at t (a Decimal point value; callers widen it)."""
    v = -PI * t.ln() / (2 * w)                      # in [0, pi] for t in [c^2, 1]
    s = sin_dec(v) if v <= 4 else None
    y = s * s / X
    H = 3 * y * (2 + y)
    if H <= 0:
        return D(0)
    return 2 * H.sqrt() / (H.sqrt() + (H + 4).sqrt())


def fx_floor(q):
    return (q.numerator << P) // q.denominator


def fx_ceil(q):
    return -((-q.numerator << P) // q.denominator)


def fx_dec_floor(x):
    return int((x * ONE).to_integral_value(rounding="ROUND_FLOOR"))


def fx_dec_ceil(x):
    return int((x * ONE).to_integral_value(rounding="ROUND_CEILING"))


def check(cert, note=print):
    c = F(cert["c"]); R = F(cert["R"])
    t = [F(x) for x in cert["nodes"]]; wts = [F(x) for x in cert["weights"]]
    l = F(cert["linear_slope"]); t0 = F(cert["second_tail_t0"]); t1 = F(cert["second_tail_t1"]); K2 = int(cert["second_tail_k"])
    ok = True
    def need(cond, msg):
        nonlocal ok
        if not cond:
            ok = False; note("FAIL: " + msg)
        return cond
    # the measure: support in [c^2, 1), positive weights, total mass sum w_j + (1 - R) = 1, no atom at 1
    need(len(t) == len(wts) and len(t) > 0, "nodes and weights differ in number")
    need(all(c * c <= x < 1 for x in t), "a node outside [c^2, 1)")
    need(all(x > 0 for x in wts), "a weight is not positive")
    need(sum(wts) == R, f"sum of weights {float(sum(wts))} is not R (total mass would not be 1)")
    need(0 < c < 1 and 0 < R < 1 and R <= t0 < t1 <= 1, "c, R, t0, t1 out of range")
    need(max(t) < t0, "an atom at or beyond t0 (the second tail argument needs every atom below t0)")
    # constants of the envelope and the admissible slope range
    w = -dec(c).ln()
    A = PI * PI / w
    X = ((A / 2).exp() - (-A / 2).exp()) ** 2 / 4
    U = -PI * dec(R).ln() / (2 * w)
    need(U * U < 6, "U^2 >= 6: the sine bound does not apply")
    a = (PI / (2 * w)) ** 2 * (1 - U * U / 6) ** 2 / X
    d = 1 - dec(R)
    l_max = 2 * (6 * a).sqrt() / ((3 * a * d * d * (2 + a * d * d)).sqrt() + (4 + 6 * a * d * d + 3 * a * a * d ** 4).sqrt())
    l_max_lo = l_max - WIDEN
    need(dec(l) < l_max_lo, f"slope {float(l)} not below the analytic bound {float(l_max)}")
    note(f"  envelope: w = {float(w):.15f}, X = {float(X):.10f}, a = {float(a):.12f}; slope bound {l_max:.15f}, "
         f"certificate slope {float(l):.15f} (ratio {float(dec(l) / l_max):.12f})")
    # atomic beta upper bounds: beta = 1 - g, g widened downwards (clamped at 0)
    beta_hi = []
    for x in t:
        g = g_of(dec(x), w, X) - WIDEN
        b = min(D(1), 1 - max(g, D(0)))
        beta_hi.append(fx_dec_ceil(b))
    # fixed-point powers with outward rounding: lo(a, b) = floor(ab / 2^P), hi(a, b) = ceil(ab / 2^P)
    lo = lambda x, y: x * y >> P
    hi = lambda x, y: -((-x * y) >> P)
    tl = [fx_floor(x) for x in t]; th = [fx_ceil(x) for x in t]
    wl = [fx_floor(x) for x in wts]; wh = [fx_ceil(x) for x in wts]
    Rl, Rh = fx_floor(R), fx_ceil(R)
    pl = [ONE] * len(t); ph = [ONE] * len(t)              # t_j^k at k = 0
    Rkl, Rkh = Rl, Rh                                     # R^(k+1) at k = 0
    # first-family tail: from the smallest K1 with w_last (K1 + 1) >= R, the atom at R alone gives A_k >= 0
    need(t[-1] == R, "the last node is not R (the first-family tail uses the atom at R)")
    K1 = 0
    while wts[-1] * (K1 + 1) < R:
        K1 += 1
    lfx = fx_floor(l)
    minI = minII = None
    s_tail_hi = None
    for k in range(1, max(K1, K2) + 1):
        pl = [lo(pl[j], tl[j]) for j in range(len(t))]
        ph = [hi(ph[j], th[j]) for j in range(len(t))]
        Rkl, Rkh = lo(Rkl, Rl), hi(Rkh, Rh)               # R^(k+1)
        if k < K1:                                        # (I): sum w_j t_j^k >= R^(k+1)/(k+1)
            s_lo = sum(lo(wl[j], pl[j]) for j in range(len(t)))
            target_hi = -((-Rkh) // (k + 1))
            if not need(s_lo >= target_hi, f"first family fails at k = {k}"):
                break
            m = (s_lo - target_hi) / target_hi
            if minI is None or m < minI[0]:
                minI = (m, k)
        if k < K2:                                        # (II): l I_k + R^(k+1)/(k+1) >= sum w_j beta_j t_j^k
            s_hi = sum(hi(hi(wh[j], ph[j]), beta_hi[j]) for j in range(len(t)))
            inv12_lo = ONE // ((k + 1) * (k + 2))
            br_hi = -((-Rkh) // (k + 1)) - lo(Rkl, Rl) // (k + 2)     # >= R^(k+1)/(k+1) - R^(k+2)/(k+2)
            I_lo = max(inv12_lo - br_hi, 0)               # <= int_R^1 (1 - t) t^k dt
            lhs_lo = lo(lfx, I_lo) + Rkl // (k + 1)
            if not need(lhs_lo >= s_hi, f"second family fails at k = {k}"):
                break
            m = (lhs_lo - s_hi) / s_hi
            if minII is None or m < minII[0]:
                minII = (m, k)
        if k == K2:                                       # (II) tail: sum_j w_j t_j^K2, beta <= 1
            s_tail_hi = sum(hi(wh[j], ph[j]) for j in range(len(t)))
    if ok:
        t0k = ONE; t0l = fx_floor(t0)
        for _ in range(K2):
            t0k = lo(t0k, t0l)
        tail_lo = lo(fx_floor((t1 - t0) * l * (1 - t1)), t0k)
        need(tail_lo > s_tail_hi, f"second-family tail at k = {K2} fails")
    if ok:
        note(f"  first family: k = 1..{K1 - 1} explicit, smallest relative margin A_k / (R^(k+1)/(k+1)) {minI[0]:.4e} at k = {minI[1]}; "
             f"from k = {K1} the atom at R alone (w_last (k + 1) >= R)")
        note(f"  second family: k = 1..{K2 - 1} explicit, smallest relative margin B_k / sum_j w_j beta_j t_j^k {minII[0]:.4e} at k = {minII[1]}; "
             f"from k = {K2}: (t1 - t0) l (1 - t1) t0^k = {tail_lo / ONE:.4e} > sum w_j t_j^k <= {s_tail_hi / ONE:.4e}")
    return ok


def main():
    t_start = time.time()
    cert = json.load(open(sys.argv[1]))
    ok = check(cert)
    print(f"RESULT C68 >= {cert['c']} = {float(F(cert['c']))}: {'all hypotheses hold' if ok else 'FAILED'}")
    if "--controls" in sys.argv:
        quiet = lambda m: None
        # 1. mass moved between two atoms (total mass kept): shifts moments, must fail one family
        bad = json.loads(json.dumps(cert)); wv = [F(x) for x in bad["weights"]]
        shift = wv[-1] / 2; wv[-1] -= shift; wv[-2] += shift
        bad["weights"] = [str(x) for x in wv]
        print("control: half the last atom's mass moved to the one before rejected:", not check(bad, quiet))
        # 2. slope raised 1% above the certificate's (still a rational): must exceed the analytic bound or fail (II)
        bad = json.loads(json.dumps(cert)); bad["linear_slope"] = str(F(cert["linear_slope"]) * F(101, 100))
        print("control: slope raised by 1% rejected:", not check(bad, quiet))
        # 3. a false second-tail cutoff K2 = 2000
        bad = json.loads(json.dumps(cert)); bad["second_tail_k"] = 2000
        print("control: second tail cut off at k = 2000 rejected:", not check(bad, quiet))
        # 4. c raised to 0.44 with the first node moved to the new c^2 (support kept valid): the envelope weakens
        bad = json.loads(json.dumps(cert)); bad["c"] = "11/25"; bad["nodes"][0] = str(F(11, 25) ** 2)
        print("control: c raised to 0.44 (first node at the new c^2) rejected:", not check(bad, quiet))
    print(f"run time {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
