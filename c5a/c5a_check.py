#!/usr/bin/env python3
"""Exact replay of the Sidon set bound C_5a <= 0.9428095 (5a page, [Akw26]) through the independent certificate of
PR #335 (A. Varlamov, sidon-c5a-certificate @ 091844e), from its PROOF.md and the certificate data alone.

The certificate gives n, m = 2n, L, D, a repair r and integers V_k (0 <= k < M/2, M = Lm). The kernels are the centred
uniform ones, p^(l)_i = 1/(2l) for n - l <= i < n + l (1 <= l <= n), with weights lambda_l = 4 l^3 / n^4 (l < n) and
(2n - 1)/n^2 (l = n). The boundary weights are w^(l)_j = W^(l)_j / D with
    W^(l)_j = floor( (1/l) sum{ V_k : 0 <= k < M/2, n - l <= j - 2k < n + l } ) + r   (0 <= j < M),   w_j = 1 (j >= M).
We check, in exact integer and rational arithmetic:
  1. the kernels and weights: nonnegative, symmetric under i -> m - 1 - i, masses 1, sum of lambda_l = 1;
  2. the covering inequalities sum_l lambda_l sum_i p^(l)_i w^(l)_{q+i} >= 1 for every 0 <= q <= M, over the common
     denominator 2 n^4 D (our window sums use prefix sums over the even sites, not pulses);
  3. a = m sum_l lambda_l sum_i (p^(l)_i)^2 and b = 1 + 2((1/m) sum_l lambda_l sum_{j<M} (w^(l)_j)^2 - L), b > 0;
  4. ab <= (942809491151629 / 10^15)^2 < (9428095 / 10^7)^2, by cross-multiplication.
By the symmetric case of Hou-Zhao's Theorem 2.1 (restated in PROOF.md, section 3), these give C_5a <= sqrt(ab).

Usage: python3 -I c5a_check.py var2026 [CERTS] [--controls]   Python 3.9+, standard library only.
"""
import hashlib
import json
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
VAR2026 = dict(rel="var2026/certificate.json", sha="fece83947509adf9c342c5e0e6a5739c3915ed62a6a5d3ec8d9141ab6e4f80ab",
               claim=(942809491151629, 10 ** 15), page=(9428095, 10 ** 7),
               stated_a=Fraction(349525, 262144),
               stated_b=Fraction(1172814639933515173834491879581015223161385476153,
                                 1759218604441600000000000000000000000000000000000))


def weights(n):
    return [None] + [Fraction(4 * l ** 3, n ** 4) for l in range(1, n)] + [Fraction(2 * n - 1, n ** 2)]


def boundaries(n, m, L, V, r):
    """W^(l)_j for 0 <= j < M + m (j >= M: D is added by the caller), by prefix sums over k."""
    M, half = L * m, len(V)
    pref = [0]
    for v in V:
        pref.append(pref[-1] + v)

    def P(x):                            # sum of V_k with 0 <= k < M/2 and 2k <= x
        return 0 if x < 0 else pref[min(half, x // 2 + 1)]
    out = {}
    for l in range(1, n + 1):
        out[l] = [(P(j - n + l) - P(j - n - l)) // l + r for j in range(M)]
    return out


def check(cert, controls_label=None, verbose=True):
    """Return (ok, facts). With controls_label set, nothing is printed except the verdict line."""
    n, m, L, D, r, V = cert["n"], cert["m"], cert["L"], cert["D"], cert["repair"], cert["renewal_numerators"]
    M = L * m
    lam = weights(n)
    facts, ok = {}, True
    say = print if verbose else (lambda *a, **k: None)
    good = m == 2 * n and len(V) == M // 2 and all(isinstance(v, int) and v >= 0 for v in V) and isinstance(r, int)
    good &= sum(lam[1:]) == 1 and all(x > 0 for x in lam[1:])
    # p^(l): 1/(2l) on [n - l, n + l): symmetric under i -> m - 1 - i, mass 1 (checked on the index sets)
    good &= all({m - 1 - i for i in range(n - l, n + l)} == set(range(n - l, n + l)) for l in range(1, n + 1))
    say(f"  data: n = {n}, m = {m}, L = {L}, M = {M}, D = 10^{len(str(D)) - 1}, repair r = {r}, {len(V)} integers V_k; "
        f"kernels symmetric with mass 1, sum of lambda_l = 1: {good}")
    ok &= good
    W = boundaries(n, m, L, V, r)
    # covering sums over the common denominator 2 n^4 D: coefficient 4 l^2 (l < n), n (2n - 1) (l = n)
    coef = [None] + [4 * l * l for l in range(1, n)] + [n * (2 * n - 1)]
    target = 2 * n ** 4 * D
    tot = [0] * (M + 1)
    for l in range(1, n + 1):
        row = W[l] + [D] * m             # w_j = 1 for j >= M
        pref = [0]
        for x in row:
            pref.append(pref[-1] + x)
        c = coef[l]
        for q in range(M + 1):
            tot[q] += c * (pref[q + n + l] - pref[q + n - l])
    fails = [q for q in range(M + 1) if tot[q] < target]
    mn = min(tot)
    arg = [q for q in range(M + 1) if tot[q] == mn]
    facts.update(fails=len(fails), min_num=mn, target=target, argmin=arg)
    say(f"  covering inequalities q = 0..{M}: {M + 1 - len(fails)} of {M + 1} hold; the least sum is "
        f"{'exactly 1' if mn == target else str(Fraction(mn, target))} at q = {arg[:5]}"
        f"{'' if len(arg) <= 5 else ' ...'}; below 1 at {fails[:5]}")
    ok &= not fails
    a = m * sum(lam[l] * Fraction(1, 2 * l) for l in range(1, n + 1))          # sum_i p_i^2 = 1/(2l)
    sq = sum(lam[l] * sum(x * x for x in W[l]) for l in range(1, n + 1))       # in units of D^2
    b = 1 + 2 * (Fraction(sq, m * D * D) - L)
    ab = a * b
    facts.update(a=a, b=b, ab=ab)
    say(f"  a = {a}; b = {float(b):.15f}... (exact fraction {b.numerator.bit_length()}/{b.denominator.bit_length()} bits), b > 0: {b > 0}")
    ok &= b > 0
    cn, cd = VAR2026["claim"]
    pn, pd = VAR2026["page"]
    below = ab.numerator * cd * cd <= cn * cn * ab.denominator
    below_page = cn * cn * pd * pd < pn * pn * cd * cd
    say(f"  ab <= (0.942809491151629)^2: {below}; (0.942809491151629)^2 < (0.9428095)^2: {below_page}")
    ok &= below and below_page
    return ok, facts


def isqrt_digits(fr, digits):
    """floor(sqrt(fr) * 10^digits) for a positive Fraction, exactly."""
    from math import isqrt
    return isqrt(fr.numerator * 10 ** (2 * digits) // fr.denominator)


def check_var(root, controls):
    t0 = time.time()
    path = os.path.join(root, VAR2026["rel"])
    raw = open(path, "rb").read()
    h = hashlib.sha256(raw).hexdigest()
    print(f"Row var2026 (PR #335, ArseniyVarlamov/sidon-c5a-certificate @ 091844e, certificate.json): sha256 {h} "
          f"({'as fetched' if h == VAR2026['sha'] else 'NOT AS PINNED'})")
    assert h == VAR2026["sha"]
    cert = json.loads(raw)
    assert cert.get("format") == "sidon-centered-uniform-v1"
    ok, f = check(cert)
    same_a, same_b = f["a"] == VAR2026["stated_a"], f["b"] == VAR2026["stated_b"]
    print(f"  a and b equal PROOF.md's exact fractions: {same_a and same_b}")
    ok &= same_a and same_b
    s = isqrt_digits(f["ab"], 20)
    print(f"  sqrt(ab) = 0.{str(s)[-20:]}... (floor to 20 digits): C_5a <= 0.942809491151629 < 0.9428095 "
          f"{'(by Hou-Zhao, Theorem 2.1, as restated in PROOF.md)' if ok else 'NOT shown'}")
    if controls:
        print("Red controls (the same checks on altered data):")
        bad = dict(cert)
        bad["repair"] = cert["repair"] - 1
        okc, fc = check(bad, verbose=False)
        c1 = not okc and fc["fails"] > 0
        print(f"  C1 repair r - 1 = {bad['repair']}: {fc['fails']} covering inequalities fail: "
              f"{'refused, as expected' if c1 else 'ACCEPTED (control failed)'}")
        bad = dict(cert)
        V = list(cert["renewal_numerators"])
        V[100] = V[100] // 2
        bad["renewal_numerators"] = V
        okc, fc = check(bad, verbose=False)
        c2 = not okc and fc["fails"] > 0
        print(f"  C2 V_100 halved: {fc['fails']} covering inequalities fail: {'refused, as expected' if c2 else 'ACCEPTED (control failed)'}")
        ab = f["ab"]
        lower = (942809491151628, 10 ** 15)
        c3 = not (ab.numerator * lower[1] ** 2 <= lower[0] ** 2 * ab.denominator)
        print(f"  C3 the claim lowered by one unit in its last digit (0.942809491151628): "
              f"{'not implied, as expected' if c3 else 'IMPLIED (control failed)'}")
        ok &= c1 and c2 and c3
    print(f"RESULT: var2026 {'PASS' if ok else 'FAIL'}")
    print(f"run time {time.time() - t0:.0f} s")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args or args[0] != "var2026":
        print(__doc__)
        sys.exit(2)
    root = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    sys.exit(0 if check_var(root, "--controls" in sys.argv) else 1)


if __name__ == "__main__":
    main()
