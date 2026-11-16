#!/usr/bin/env python3
"""5a, upper bound [COT26] C_5a <= 0.94304 (Carter, O'Bryant and Tao, arXiv 2610.36207v1, Corollary 1.2 and
Lemma 2.2), program A: the WT-triple witnesses checked exactly against Definition 2.2, with the paper's index formula.

Definition 2.2 (read through the paper's HTML): (upsilon, alpha, beta) is a WT-triple if for some sigma >= 2 upsilon there
are random integrable w (support [0, upsilon]) and t (support [-upsilon, sigma]) with (i) t = 1 on (0, sigma - upsilon),
(ii) int w = 1, (iii) E ||t||_2^2 <= sigma + alpha, (iv) E ||w||_2^2 <= beta, (v) E (w autocorrelated)(x) >= 0 for all x,
(vi) E (w * t)(x) >= 1 on (0, sigma). Strong: w, t >= 0 (which gives (v)). Section 6.1: with delta = upsilon/N and
sigma = 2 upsilon, atom i (probability p_i) has w_i = sum_j g_j 1_[j delta, (j+1) delta) and
t_i(x) = sum_j h_(3N-1-j) 1_[j delta, (j+1) delta)(x + upsilon) (the stored vector h is mirrored).
Data: dcartermath/sidon-error @ ff5b032, asymptotic/ (upsilon = 10, K = 6 atoms, N = 200) and explicit/ (upsilon = 3,
K = 5, N = 32); G.txt, H.txt, p.txt as exact rationals (SHA-256 in SHA256SUMS). verify.py and the search code were
neither read nor run.
Program A: (w_i * t_i)(j delta) = delta sum_m g_m T_(j - m - 1 + N), T_k = h_(3N-1-k) (derived from the definitions: for
y in [m delta, (m+1) delta), j delta - y lies in the t-cell with k delta - upsilon = (j - m - 1) delta); w * t is linear
between knots j delta (both are step functions on the same grid), so (vi) holds iff it holds at j = 0..2N (we check the
closed interval [0, sigma], which is stronger). alpha = E ||t||^2 - sigma and beta = E ||w||^2 exactly; then
2 sqrt(alpha beta) and 4 (1 + upsilon sqrt(beta/alpha))^2 against the paper's 1.88607 and 483 (resp. 1.8898801 and 30),
and the row's sqrt(alpha beta) <= 0.94304. Theorem 2.3 (WT-triple => the diameter bound => C_5a <= sqrt(alpha beta)) is
read, not checked. Red controls (--control): the claim 0.94303; one table height in the middle third set to 0.999
(condition (i)); the mirroring ignored (t read in stored order; (vi) must fail). Standard library only.
"""
import os
import sys
from fractions import Fraction as F
from math import isqrt

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else os.path.join(os.path.dirname(os.path.abspath(__file__)), "certs")
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None


def load(d):
    rd = lambda f: open(os.path.join(ROOT, d, f)).read()
    G = [[F(x) for x in l.strip().split(",")] for l in rd("G.txt").splitlines() if l.strip()]
    H = [[F(x) for x in l.strip().split(",")] for l in rd("H.txt").splitlines() if l.strip()]
    p = [F(x) for x in rd("p.txt").strip().split(",")]
    u = F(rd("upsilon.txt").strip())
    return G, H, p, u


def check(d, claims):
    G, H, p, u = load(d)
    K = len(G); N = len(G[0]); delta = u / N; sigma = 2 * u
    ok = len(H) == K == len(p) and all(len(g) == N for g in G) and all(len(h) == 3 * N for h in H)
    ok &= all(x >= 0 for x in p) and sum(p) == 1
    if CONTROL == "middle" and d == "asymptotic":
        H[0][N + 7] = F(999, 1000)
    T = []
    for h in H:
        T.append(list(h) if CONTROL == "mirror" else [h[3 * N - 1 - k] for k in range(3 * N)])
    c1 = all(T[i][k] == 1 for i in range(K) for k in range(N, 2 * N))                   # t = 1 on [0, upsilon)
    c2 = all(delta * sum(g) == 1 for g in G)
    strong = all(x >= 0 for g in G for x in g) and all(x >= 0 for h in H for x in h)
    Et2 = sum(p[i] * delta * sum(x * x for x in H[i]) for i in range(K))
    Ew2 = sum(p[i] * delta * sum(x * x for x in G[i]) for i in range(K))
    alpha = Et2 - sigma; beta = Ew2
    # (vi) at the knots j delta, j = 0..2N
    worst = None; worst_j = None
    for j in range(2 * N + 1):
        v = F(0)
        for i in range(K):
            if p[i] == 0:
                continue
            s = F(0)
            g = G[i]; Ti = T[i]
            for m in range(N):
                k = j - m - 1 + N
                if 0 <= k < 3 * N:
                    s += g[m] * Ti[k]
            v += p[i] * delta * s
        if worst is None or v < worst:
            worst, worst_j = v, j
    c6 = worst >= 1
    print("== %s: upsilon = %s, K = %d atoms, N = %d, delta = %s, sigma = %s" % (d, u, K, N, delta, sigma))
    print("   p_i >= 0, sum p_i = 1, shapes: %s; (i) t = 1 on (0, sigma - upsilon): %s; (ii) int w_i = 1: %s; strong (w, t >= 0): %s"
          % (ok, c1, c2, strong))
    print("   (vi) min over knots of E(w * t) = %s at x = %s delta (>= 1: %s)" % (dec(worst, 15), worst_j, c6))
    print("   alpha = E||t||^2 - sigma = %s...; beta = E||w||^2 = %s..." % (dec(alpha, 15), dec(beta, 15)))
    print("   beta = %s/%s (%d and %d digits)" % (str(beta.numerator)[:20] + "...", str(beta.denominator)[:20] + "...",
          len(str(beta.numerator)), len(str(beta.denominator))))
    ab = alpha * beta
    res = ok and c1 and c2 and strong and c6
    if d == "asymptotic":
        # Appendix A.3 prints beta_0 as one fraction; compare exactly
        b0 = F(1772130384572488116897712811968136530510904755768535954960651,
               1882725390027397260273972602940280969126497650532485399838720)
        print("   beta equals the fraction printed in Appendix A.3: %s" % (beta == b0))
        res &= beta == b0
    b2, k0, row = claims
    b_ok = 4 * ab <= F(b2) ** 2                                   # 2 sqrt(alpha beta) <= b2
    # k0 = 4 (1 + upsilon sqrt(beta/alpha))^2 <= k0max  <=>  upsilon sqrt(beta/alpha) <= sqrt(k0max)/2 - 1
    r = F(isqrt(k0 * 10 ** 40), 10 ** 20) / 2 - 1                # a lower bound for sqrt(k0)/2 - 1
    k_ok = r > 0 and u * u * beta <= r * r * alpha
    s_lo = F(isqrt(int(ab * 10 ** 60)), 10 ** 30)
    print("   2 sqrt(alpha beta) >= %s... <= %s: %s; 4 (1 + upsilon sqrt(beta/alpha))^2 <= %d: %s"
          % (dec(2 * s_lo, 12), b2, b_ok, k0, k_ok))
    res &= b_ok and k_ok
    if row is not None:
        row_ok = ab <= F(row) ** 2
        print("   the row: sqrt(alpha beta) = %s... <= %s: %s" % (dec(s_lo, 12), row, row_ok))
        res &= row_ok
    print("   key %s %s %s" % (d, dec(alpha, 20), dec(beta, 20)))
    return res


def dec(q, k):
    neg = q < 0; q = abs(q)
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return ("-" if neg else "") + s[:-k] + "." + s[-k:]


def main():
    row = "0.94303" if CONTROL == "claim" else "0.94304"
    ok = check("asymptotic", ("1.88607", 483, row))
    if CONTROL is None:
        ok &= check("explicit", ("1.8898801", 30, None))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
