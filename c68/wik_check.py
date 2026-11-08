#!/usr/bin/env python3
"""Independent check of the 68a row [Wik2026], C_68 >= 0.4263 (Wikstrom, arXiv:2607.17748v1, Theorem 1.1; the row was
added with PR #239), from the paper's LaTeX source alone: its Arb script was neither read nor run.

The criterion (Lemma 3.1 of the paper, re-read): c is admissible if a probability measure mu on [c^2, 1] with no atom
at 1 has  int t^k dmu >= 1/(k+1) (k >= 0)  and  int beta_c(sqrt t) t^k dmu <= 1/(k+1) (k >= 1), where, from Wang's
annular estimate (Proposition 2.1, taken as given),
  K_c(r) = 2r(1 + c/r^2) prod_{n>=1} (1 + r^2 c^(2n-1))(1 + r^-2 c^(2n+1))(1 + c^(2n))^2
                                   / ((1 + r^2 c^(2n-2))(1 + r^-2 c^(2n))(1 + c^(2n-1))^2),
  G_c = (1 + sqrt((1 + 3K^2)/(1 - K^2)))/2,  beta_c = 1 - 1/G_c,  beta_c(c) = 1.
For mu = sum_j w_j delta_(t_j) + 1_[R,1] dt (the paper's Table 1: eight terminating decimals, c = 0.4263, R = 0.99) the
two families are A_k = sum w_j t_j^k - R^(k+1)/(k+1) >= 0 and B_k = int_R^1 t^k/G dt + R^(k+1)/(k+1)
- sum w_j beta_j t_j^k >= 0 (re-derived here by integrating the Lebesgue part).

What we compute, with our own code (fixed point at 2^-320 with every operation rounded outward):
  1. the data: t_1 = c^2, t_8 = R, sum w_j = R (exact rationals), so mu is a probability measure with no atom at 1;
  2. K_c(sqrt t) at a point, from above and below: the product truncated after N = 28 factors, times exp(+-E) with
     our own tail bound E = (1 + 3c + 2c^2 + (c^3 + c^2)/L) c^(2N)/(1 - c^2) for t >= L (from |log(1 + x)| <= x on the
     omitted factors). K increases on [R, 0.99999]: we rewrite d/dt log K_N(sqrt t) = (1 - t) E(t) + b_2N - b_(2N+1)
     with E(t) = 1/(2t(1+t)) + sum_{m=1}^{2N-1} (-1)^m c^m (1+t)/((1 + t c^m) t (t + c^m)) and b_m = c^m/(t(t + c^m))
     (the large terms 1/(2t) and -1/(1+t) cancel to first order at t = 1, which the factor (1 - t) makes explicit),
     bound E from below on 1000 pieces per block (each term is monotone in t), and add the derivative tail
     (1 + c + (c^3 + c^2)/L^2) c^(2N)/(1 - c^2);
  3. A_k > 0 for 1 <= k <= 2183 (A_0 = 0), and the atom at R alone for k >= 2184 (2185 w_8 > R);
  4. beta_j from above at the nodes (beta_1 = 1), and int_R^1 t^k/G dt from below on the paper's partition
     ([R, 0.999], [0.999, 0.9999], [0.9999, 0.99999] in 30000, 5000 and 500 pieces; [0.99999, 1] discarded), with
     1/G on each piece bounded below through K at its right end, from above; B_k > 0 for 1 <= k <= 1299;
  5. k >= 1300: 1/G > 0.0000225 on [0.998, 0.999], so B_k >= 0.0000225 (0.001) 0.998^k - sum w_j t_j^k, positive at
     k = 1300 and thereafter (the positive term shrinks by 0.998 per step, every negative one by at most 0.99).
Usage: python3 -I wik_check.py [--controls]. Python 3.9+, standard library only."""
import sys
import time
from fractions import Fraction as F
from math import isqrt

P = 320
ONE = 1 << P
C = F(4263, 10000)
R = F(99, 100)
NODES = [("0.1817316900000000", "0.3924240372921149"), ("0.5786106116843371", "0.3349420206616465"),
         ("0.8163605367760075", "0.1555309603886983"), ("0.9200761835179418", "0.0649137961077427"),
         ("0.9633238309788337", "0.0271186307684793"), ("0.9813028389745606", "0.0109963471767101"),
         ("0.9881242650409469", "0.0036210554863853"), ("0.9900000000000000", "0.0004531521182229")]
NFACT = 28
BLOCKS = [(F(99, 100), F(999, 1000), 30000), (F(999, 1000), F(9999, 10000), 5000), (F(9999, 10000), F(99999, 100000), 500)]


# ------------------------------------------------------------------ fixed-point intervals (pairs of ints at 2^-P)
def q_lo(x):
    return x.numerator * ONE // x.denominator


def q_hi(x):
    return -((-x.numerator * ONE) // x.denominator)


def mul(a, b):        # nonnegative intervals
    assert a[0] >= 0 and b[0] >= 0
    return (a[0] * b[0]) >> P, -((-(a[1] * b[1])) >> P)


def div(a, b):        # nonnegative a, positive b
    assert a[0] >= 0 and b[0] > 0
    return (a[0] << P) // b[1], -((-(a[1] << P)) // b[0])


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def sqrt_iv(a):
    return isqrt(a[0] << P), isqrt(a[1] << P) + 1


def exp_small(x_hi):  # an upper bound for exp(x) with 0 <= x <= 2^-20 (fixed point): 1 + x + x^2
    return ONE + x_hi + ((x_hi * x_hi) >> P) + 1


CP = {}               # powers of c as intervals


def cpow(m):
    if m not in CP:
        CP[m] = (q_lo(C ** m), q_hi(C ** m))
    return CP[m]


def tail_bound(L):    # E(L) from above
    e = (1 + 3 * C + 2 * C * C + (C ** 3 + C * C) / L) * C ** (2 * NFACT) / (1 - C * C)
    return q_hi(e)


def K_bounds(u, v):
    """[K_lo, K_hi] for K_c(sqrt t), t in [u, v] (Fractions, c < u <= v < 1)."""
    U, V = (q_lo(u), q_hi(u)), (q_lo(v), q_hi(v))
    one = (ONE, ONE)
    su, sv = sqrt_iv(U), sqrt_iv(V)
    c = cpow(1)
    # prefactor 2 (sqrt t + c / sqrt t), increasing in t for t > c
    hi = mul((2 * ONE, 2 * ONE), add(sv, div(c, sv)))[1]
    lo = mul((2 * ONE, 2 * ONE), add(su, div(c, su)))[0]
    k_hi, k_lo = (hi, hi), (lo, lo)
    for n in range(1, NFACT + 1):
        a1, a2, a3 = cpow(2 * n - 1), cpow(2 * n + 1), cpow(2 * n)
        b1 = cpow(2 * n - 2)
        num_hi = mul(mul(add(one, mul(V, a1)), add(one, div(a2, U))), mul(add(one, a3), add(one, a3)))[1]
        den_lo = mul(mul(add(one, mul(U, b1)), add(one, div(a3, V))), mul(add(one, a1), add(one, a1)))[0]
        num_lo = mul(mul(add(one, mul(U, a1)), add(one, div(a2, V))), mul(add(one, a3), add(one, a3)))[0]
        den_hi = mul(mul(add(one, mul(V, b1)), add(one, div(a3, U))), mul(add(one, a1), add(one, a1)))[1]
        k_hi = mul(k_hi, div((num_hi, num_hi), (den_lo, den_lo)))
        k_lo = mul(k_lo, div((num_lo, num_lo), (den_hi, den_hi)))
    eh = exp_small(tail_bound(u))
    K_hi = mul(k_hi, (eh, eh))[1]
    K_lo = div(k_lo, (eh, eh))[0]
    return K_lo, K_hi


def inv_G_lo(K_hi):
    """a lower bound for 1/G from an upper bound for K (G increases with K); needs K_hi < 1."""
    assert K_hi < ONE
    K = (K_hi, K_hi)
    K2 = mul(K, K)
    num = (ONE + 3 * K2[0], ONE + 3 * K2[1])
    den = (ONE - K2[1], ONE - K2[0])
    rt = sqrt_iv(div(num, den))
    G_hi = (ONE + rt[1] + 1) // 2 + 1
    return (ONE * ONE) // G_hi


def E_lo(u, v):
    """a lower bound over t in [u, v] for E(t) = 1/(2t(1+t)) + sum_{m=1}^{2N-1} (-1)^m f_m(t),
    f_m(t) = c^m (1 + t)/((1 + t c^m) t (t + c^m)); every f_m and 1/(2t(1+t)) decrease in t."""
    U, V = (q_lo(u), q_hi(u)), (q_lo(v), q_hi(v))
    one = (ONE, ONE)
    first = div(one, mul((2 * ONE, 2 * ONE), mul(V, add(one, V))))[0]
    tot = first
    for m in range(1, 2 * NFACT):
        cm = cpow(m)
        if m % 2 == 0:        # + f_m(t) >= f_m(v)
            val = div(mul(cm, add(one, V)), mul(mul(add(one, mul(V, cm)), V), add(V, cm)))[0]
            tot += val
        else:                 # - f_m(t) >= - f_m(u)
            val = div(mul(cm, add(one, U)), mul(mul(add(one, mul(U, cm)), U), add(U, cm)))[1]
            tot -= val
    return tot


def deriv_tail(L):    # the derivative tail for the omitted factors, from above
    e = (1 + C + (C ** 3 + C * C) / (L * L)) * C ** (2 * NFACT) / (1 - C * C)
    return q_hi(e)


def monotone(blocks, n):
    """K_c(sqrt t) increases on each block: d/dt log K_N = (1 - t) E(t) + b_2N - b_(2N+1) with b_2N > b_(2N+1), and
    |d/dt log K - d/dt log K_N| <= deriv_tail; so (1 - v) E_lo > deriv_tail on every piece suffices."""
    out = []
    for a, b, _ in blocks:
        h = (b - a) / n
        worst = None
        for i in range(n):
            u, v = a + i * h, a + (i + 1) * h
            e = E_lo(u, v)
            d = ((q_lo(1 - v) * e) >> P) - deriv_tail(u)
            worst = d if worst is None or d < worst else worst
        out.append(worst)
    return out


def pieces(blocks):
    """breakpoints and, on each piece [u, v], 1/G from below through K(v) from above (K increases there)."""
    xs, gs = [], []
    for a, b, n in blocks:
        h = (b - a) / n
        for i in range(n):
            u, v = a + i * h, a + (i + 1) * h
            K_lo, K_hi = K_bounds(v, v)
            xs.append(u)
            gs.append(inv_G_lo(K_hi))
    xs.append(blocks[-1][1])
    return xs, gs


def second_family(xs, gs, beta_hi, w, t, kmax):
    """B_k from below for 1 <= k <= kmax; returns (all positive, smallest relative margin, its k, first k not certified).
    int_R^0.99999 t^k/G >= sum_p g_p (v_p^(k+1) - u_p^(k+1))/(k+1) = sum_i c_i x_i^(k+1)/(k+1), c_i = g_(i-1) - g_i."""
    coef = [(gs[i - 1] if i > 0 else 0) - (gs[i] if i < len(gs) else 0) for i in range(len(xs))]
    X = [(q_lo(x), q_hi(x)) for x in xs]
    powx = list(X)
    T = [(q_lo(x), q_hi(x)) for x in t]
    W = [(q_lo(x), q_hi(x)) for x in w]
    cur = list(T)
    RR = (q_lo(R), q_hi(R))
    Rk1 = mul(RR, RR)
    worst = kw = kf = None
    for k in range(1, kmax + 1):
        powx = [mul(p, x) for p, x in zip(powx, X)]
        integ = 0
        for cf, p in zip(coef, powx):
            integ += (cf * p[0]) >> P if cf >= 0 else (cf * p[1]) >> P     # each product rounded down
        integ_lo = integ // (k + 1)
        rk = Rk1[0] // (k + 1)
        neg = sum(mul(mul(wj, (bh, bh)), x)[1] for wj, bh, x in zip(W, beta_hi, cur))
        B_lo = integ_lo + rk - neg
        rel = F(B_lo * (k + 1), Rk1[1])
        if worst is None or rel < worst:
            worst, kw = rel, k
        if B_lo <= 0 and kf is None:
            kf = k
        cur = [mul(a, x) for a, x in zip(cur, T)]
        Rk1 = mul(Rk1, RR)
    return kf is None, worst, kw, kf


def main():
    t0 = time.time()
    controls = "--controls" in sys.argv
    ok = True
    t = [F(a) for a, _ in NODES]
    w = [F(b) for _, b in NODES]
    good = t[0] == C * C and t[-1] == R and sum(w) == R and all(x < y for x, y in zip(t, t[1:]))
    print(f"68a row [Wik2026]: C_68 >= 0.4263 (c = 4263/10000, R = 99/100, eight atoms from Table 1)")
    print(f"  t_1 = c^2 = {t[0]}, t_8 = R, sum w_j = {sum(w)} = R, nodes increasing: {good}")
    ok &= good
    # 3. the first family
    worst, kw = None, None
    pw = [(q_lo(x), q_hi(x)) for x in t]
    cur = [(ONE, ONE)] * 8
    Rk1 = (q_lo(R), q_hi(R))
    good = True
    for k in range(1, 2184):
        cur = [mul(a, b) for a, b in zip(cur, pw)]
        Rk1 = mul(Rk1, (q_lo(R), q_hi(R)))           # R^(k+1)
        s_lo = sum(mul((q_lo(wj), q_hi(wj)), x)[0] for wj, x in zip(w, cur))
        bound_hi = -((-Rk1[1]) // (k + 1))
        A_lo = s_lo - bound_hi
        rel = F(A_lo * (k + 1), Rk1[1])
        if worst is None or rel < worst:
            worst, kw = rel, k
        good &= A_lo > 0
    tail = 2185 * w[-1] > R
    print(f"  A_k > 0 for 1 <= k <= 2183: {good}; smallest A_k / (R^(k+1)/(k+1)) = {float(worst):.10e} at k = {kw} "
          f"(the paper: > 6.6452117503e-6 at k = 23); 2185 w_8 = {float(2185 * w[-1]):.16f} > R: {tail}")
    ok &= good and tail and kw == 23
    # 4. beta at the nodes, from above
    beta_hi = [ONE]
    for x in t[1:]:
        K_lo, K_hi = K_bounds(x, x)
        beta_hi.append(ONE - inv_G_lo(K_hi))
    paper = ["0.989685929536115", "0.995533783756722", "0.998127951092143", "0.999157174953505", "0.999573932287887",
             "0.999730257446456", "0.999773067105237"]
    good = all(F(b, ONE) < F(p) for b, p in zip(beta_hi[1:], paper))
    print(f"  beta_c(sqrt t_j) from above, j = 2..8: {', '.join(f'{float(F(b, ONE)):.15f}' for b in beta_hi[1:])}")
    print(f"  each below the paper's stated enclosure: {good}")
    ok &= good
    # K increases on [R, 0.99999]: the derivative of log K from below, block by block
    t1 = time.time()
    mins = monotone(BLOCKS, 1000)
    good = all(m > 0 for m in mins)
    print(f"  d/dt log K_c(sqrt t) from below on [R, 0.999], [0.999, 0.9999], [0.9999, 0.99999] (1000 pieces each): "
          f"{', '.join(f'{float(F(m, ONE)):.10e}' for m in mins)} (the paper: 3.6694461888e-7, 3.6647630565e-8, "
          f"3.6642942322e-9); all positive: {good} ({time.time() - t1:.0f} s)")
    ok &= good
    # the integral from below: breakpoints and 1/G per piece
    t1 = time.time()
    xs, gs = pieces(BLOCKS)
    print(f"  {len(gs)} pieces on [R, 0.99999]: 1/G from below, K < 1 on every piece; {time.time() - t1:.0f} s")
    t2 = time.time()
    good, worst, kw, kf = second_family(xs, gs, beta_hi, w, t, 1299)
    print(f"  B_k > 0 for 1 <= k <= 1299: {good}; smallest B_k / (R^(k+1)/(k+1)) = {float(worst):.10e} at k = {kw} "
          f"(the paper: > 6.0214605308e-6 at k = 522); {time.time() - t2:.0f} s")
    ok &= good and kw == 522
    # 5. the tail k >= 1300
    K_lo, K_hi = K_bounds(F(999, 1000), F(999, 1000))     # K increases on [R, 0.999], so K <= K(0.999) on [0.998, 0.999]
    g = inv_G_lo(K_hi)
    good1 = F(g, ONE) > F(225, 10 ** 7)
    pos = q_lo(F(225, 10 ** 7) * F(1, 1000) * F(998, 1000) ** 1300)
    negs = sum(q_hi(wj * x ** 1300) for wj, x in zip(w, t))
    good2 = pos > negs and all(x <= F(99, 100) for x in t)
    print(f"  1/G > 0.0000225 on [0.998, 0.999] (through K(0.999); our bound {float(F(g, ONE)):.7e}): {good1}; at k = 1300 "
          f"{float(F(pos, ONE)):.15e} > {float(F(negs, ONE)):.15e}, and every t_j <= 0.99 < 0.998: {good2}")
    ok &= good1 and good2
    if controls:
        print("Red controls (each must not be certified):")
        w2 = list(w)
        w2[1] = w2[1] * F(99, 100)
        cur = [(q_lo(x), q_hi(x)) for x in t]
        kf = None
        Rk1 = mul((q_lo(R), q_hi(R)), (q_lo(R), q_hi(R)))
        for k in range(1, 2184):
            s_hi = sum(mul((q_lo(wj), q_hi(wj)), x)[1] for wj, x in zip(w2, cur))
            if s_hi < Rk1[0] // (k + 1):
                kf = k
                break
            cur = [mul(a, (q_lo(x), q_hi(x))) for a, x in zip(cur, t)]
            Rk1 = mul(Rk1, (q_lo(R), q_hi(R)))
        print(f"  w_2 cut by 1%: A_k < 0 (rigorously) at k = {kf}: {'as expected' if kf else 'NOT REFUSED'}")
        ok &= kf is not None
        coarse = [(a, b, n // 100) for a, b, n in BLOCKS]
        xs2, gs2 = pieces(coarse)
        good, worst, kw, kf = second_family(xs2, gs2, beta_hi, w, t, 1299)
        print(f"  the partition 100 times coarser (300, 50, 5 pieces): B_k not certified from k = {kf}: "
              f"{'as expected' if not good else 'CERTIFIED (control failed)'}")
        ok &= not good
        b2 = list(beta_hi)
        b2[7] = b2[7] + (ONE - b2[7]) // 2          # beta at R raised halfway to 1
        good, worst, kw, kf = second_family(xs, gs, b2, w, t, 1299)
        print(f"  beta_c(sqrt t_8) raised halfway to 1: B_k not certified from k = {kf}: "
              f"{'as expected' if not good else 'CERTIFIED (control failed)'}")
        ok &= not good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
