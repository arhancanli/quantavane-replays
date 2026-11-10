#!/usr/bin/env python3
"""Replay of #284 (DennisKoshta, 26a, comment only): explicit f_d on {+-1}^{3d} of degree d, values in {0, +-1}, with
(2^d - 1) 2^d Fourier coefficients +-2^-d, so BH^{<=d} >= R_d = 2 (1 - 2^-d)^{(d+1)/(2d)} > 2^{(d-1)/d} for d >= 2.
Written from WRITEUP.md (Sections 2 and 3) and problems/bh_cube/FORMAT.md of DennisKoshta/c26a-bohnenblust-hille-cube-
family @ 2557cdb and the six instance files; core.py, cycle.py and verify/ were neither read nor run.

Two routes:
  A. the six certificates as FORMAT.md defines them: integer coefficients a_S on subsets S of {0..n-1}; F(x) =
     sum a_S x^S at all 2^n points by our own exact fast Walsh-Hadamard transform, M = max |F|, and the claim
     R = (sum |a_S|^p)^{1/p} / M >= c with p = 2d/(d+1). Every a_S here is +-1, so R >= c is the exact integer
     comparison N^(d+1) >= (c M)^(2d) (N = number of terms, c the claimed decimal as a rational);
  B. our own f_d from Section 3 (z_i = i, a_i = d + 2i, b_i = d + 2i + 1; epsilon_i = e_{z_i} + e_{z_{i+1}},
     pi_i = e_{a_i} + e_{b_i}, c = 0 for odd d and e_0 for even d, for every proper J of Z_d the coset beta_J + V_J with
     coefficients 2^-d (-1)^{sum lambda_i w_J(i)} by Lemma 1, alpha = sum lambda_i v_i, v_i = epsilon_i (i in J) or
     pi_i (i not in J)), for d = 2..7: the (2^d - 1) 2^d frequencies are distinct (Claim 1) and of weight <= d, F = 2^d f
     takes only the values 0, +-2^d (Claim 2), and F != 0 on a fraction 1 - 2^-d of the cube; for d = 5, 6, 7 it is
     compared term by term with the instance file.
R_d is compared with the closed form and with the address value 2^{(d-1)/d} exactly (R_d^(2d) = 4^d (1-2^-d)^(d+1),
rational). Red controls (each must be refused): one coefficient's sign flipped; the claim raised by 10^-12; a term of
degree d + 1 added; the parity choice c = 0 for even d (then two live sets meet and F takes other values).

Usage: python3 bh_check.py INSTANCE_DIR [--controls]     (standard library only; about two minutes)
"""
import hashlib
import json
import os
import re
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations

FILES = {  # file: (sha256 the writeup states, n, d, terms, sup)
    "d2_n5.json": ("13691d36694b3bc1fcc2ed8ba62d785f97c3bb0900c0e9f5f2d91f0a7444385d", 5, 2, 12, 4),
    "d3_n8.json": ("af475381723e51c82340c9d618f99cf9769af2d0c0062e1ae4cb3eaa91fb2e53", 8, 3, 56, 8),
    "d4_n11.json": ("4d303d60ce61ddea5435648bec1cf7410c8bea7992127c906727c7f9122277ad", 11, 4, 240, 16),
    "cycle_d5_n15.json": ("bbbd7e89949ffa937ad05f2dd57842b6b5c608d82573121e5a6024d6c2dcb8c5", 15, 5, 992, 32),
    "cycle_d6_n18.json": ("f9b01bcf183adb0aabcead0de136f0512d29f470da06f6dc58eda369c9637459", 18, 6, 4032, 64),
    "cycle_d7_n21.json": ("5b22ba5c5b60dd5f6074ee6547f26f3727b2cf5064dd40998d578eb201a929d3", 21, 7, 16256, 128),
}


def need(c, msg):
    if not c:
        raise AssertionError(msg)


def wht_values(n, coeffs):
    """F at every point: coeffs maps a bitmask S to an integer; returns the list F[m], x_i = -1 iff bit i of m."""
    v = [0] * (1 << n)
    for S, a in coeffs.items():
        v[S] += a
    h = 1
    while h < len(v):
        for start in range(0, len(v), 2 * h):
            for j in range(start, start + h):
                x, y = v[j], v[j + h]
                v[j], v[j + h] = x + y, x - y
        h *= 2
    return v


def claim_holds(N, M, d, c):
    """R = N^((d+1)/(2d)) / M >= c (every |a_S| = 1)  <=>  N^(d+1) >= (c M)^(2d), exactly"""
    return Fr(N) ** (d + 1) >= (c * M) ** (2 * d)


def check_file(path, mutate=None, log=print):
    raw = open(path, "rb").read()
    name = os.path.basename(path)
    sha = hashlib.sha256(raw).hexdigest()
    d_ = json.loads(raw.decode("utf-8"))
    n, d = d_["n"], d_["d"]
    ok = d_["format"] == "bh_cube/v1" and d_["constant"] == "C26a" and isinstance(n, int) and isinstance(d, int)
    coeffs = {}
    for S, a in d_["terms"]:
        ok &= isinstance(a, str) and re.fullmatch(r"-?[1-9][0-9]*", a) is not None
        ok &= all(isinstance(i, int) for i in S) and list(S) == sorted(set(S)) and all(0 <= i < n for i in S)
        ok &= len(S) <= d
        mask = sum(1 << i for i in S)
        ok &= mask not in coeffs
        coeffs[mask] = int(a)
    if mutate:
        coeffs, d = mutate(coeffs, d)
        ok &= all(bin(m).count("1") <= d for m in coeffs)
    F = wht_values(n, coeffs)
    M = max(abs(x) for x in F)
    N = len(coeffs)
    ok &= M > 0 and all(abs(a) == 1 for a in coeffs.values())
    c = Fr(d_["claimed_lower_bound"])
    sup_ok = str(M) == d_.get("claimed_sup", str(M))
    holds = claim_holds(N, M, d, c)
    vals = sorted(set(F))
    frac_nonzero = Fr(sum(1 for x in F if x != 0), len(F))
    want = FILES.get(name)
    log("  %s: sha256 %s; n = %d, d = %d, %d terms, M = %d (claimed_sup %s), values of F %s, F != 0 on %s of the "
        "cube; claim %s: %s" % (name, "as stated" if want and sha == want[0] else "DIFFERENT", n, d, N, M,
                                d_.get("claimed_sup"), vals if len(vals) <= 5 else "%d distinct" % len(vals),
                                frac_nonzero, d_["claimed_lower_bound"], "holds" if holds else "FAILS"))
    return ok and sup_ok and holds and (want is None or (sha == want[0] and (n, d, N, M) == want[1:])), (coeffs, n, d, c)


def build_fd(d, parity_ok=True):
    """Section 3: the coefficients of F = 2^d f_d as a dict bitmask -> +-1 (n = 3d)."""
    z = list(range(d))
    a = [d + 2 * i for i in range(d)]
    b = [d + 2 * i + 1 for i in range(d)]
    eps = [(1 << z[i]) ^ (1 << z[(i + 1) % d]) for i in range(d)]
    pi = [(1 << a[i]) | (1 << b[i]) for i in range(d)]
    cvec = [0] * d
    if d % 2 == 0 and parity_ok:
        cvec[0] = 1
    coeffs = {}
    clash = 0
    for r in range(0, d):
        for J in combinations(range(d), r):          # proper subsets of Z_d (size < d)
            Js = set(J)
            w = [((1 if (i - 1) % d in Js else 0) + cvec[i]) % 2 for i in range(d)]
            gens = [eps[i] if i in Js else pi[i] for i in range(d)]
            # beta_J: e_{a_i} for i not in J; for each path of J (maximal runs of consecutive i in J), e_u if its
            # number of edges m is odd (C_P = e_u + E_P), nothing if m is even (C_P = E_P)
            beta = 0
            for i in range(d):
                if i not in Js:
                    beta |= 1 << a[i]
            for i in Js:
                if (i - 1) % d in Js:
                    continue                          # not the first edge of its run
                m, k = 0, i
                while k in Js and m < d:
                    m += 1
                    k = (k + 1) % d
                if m % 2 == 1:
                    beta ^= 1 << z[i]                 # u = the run's first vertex z_i
            for lam in range(1 << d):
                alpha, sgn = 0, 0
                for i in range(d):
                    if lam >> i & 1:
                        alpha ^= gens[i]
                        sgn ^= w[i]
                key = beta ^ alpha
                if key in coeffs:
                    clash += 1
                coeffs[key] = coeffs.get(key, 0) + (-1 if sgn else 1)
    return coeffs, clash


def main():
    t0 = time.time()
    root = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    controls = "--controls" in sys.argv
    print("26a (#284): the six certificates (route A) and our own f_d from Section 3 (route B)")
    ok = True
    files = {}
    print("Route A, the certificates as FORMAT.md defines them:")
    for name in FILES:
        good, data = check_file(os.path.join(root, name))
        files[name] = data
        ok &= good
    print("Route B, f_d built from Section 3:")
    for d in range(2, 8):
        coeffs, clash = build_fd(d)
        n = 3 * d
        N = len(coeffs)
        weights_ok = all(bin(m).count("1") <= d for m in coeffs)
        F = wht_values(n, coeffs)
        vals = sorted(set(F))
        frac = Fr(sum(1 for x in F if x != 0), len(F))
        good = clash == 0 and N == (2 ** d - 1) * 2 ** d and weights_ok and vals == [-2 ** d, 0, 2 ** d] and \
            frac == 1 - Fr(1, 2 ** d) and all(abs(v) == 1 for v in coeffs.values())
        # R_d^(2d) = N^(d+1) / 2^(2d^2) = 4^d (1 - 2^-d)^(d+1); compare with the address value 2^((d-1)/d)
        R2d = Fr(N) ** (d + 1) / Fr(2) ** (2 * d * d)
        closed = Fr(4) ** d * (1 - Fr(1, 2 ** d)) ** (d + 1)
        addr = Fr(2) ** (2 * (d - 1))                       # (2^((d-1)/d))^(2d)
        same = ""
        name = "cycle_d%d_n%d.json" % (d, 3 * d)
        if name in files:
            fc = files[name][0]
            same = "; equal to %s term by term: %s" % (name, fc == coeffs)
            good &= fc == coeffs
        good &= R2d == closed and R2d > addr
        print("  d = %d: %d frequencies, distinct: %s, weights <= d: %s; F takes the values %s, non-zero on %s of the cube; "
              "R_d^(2d) = 4^d (1-2^-d)^(d+1): %s; R_d > 2^((d-1)/d): %s (R_d = %.6f)%s"
              % (d, N, clash == 0, weights_ok, vals, frac, R2d == closed, R2d > addr,
                 float(R2d) ** (1 / (2 * d)), same))
        ok &= good
    if controls:
        print("Red controls (each must be refused):")
        name = "cycle_d5_n15.json"

        def flip(cf, d):
            cf = dict(cf)
            k = sorted(cf)[17]
            cf[k] = -cf[k]
            return cf, d
        g1, _ = check_file(os.path.join(root, name), mutate=flip, log=lambda s: None)
        coeffs, n, d, c = files[name]
        F = wht_values(n, flip(coeffs, d)[0])
        print("  1. one coefficient's sign flipped in %s: M = %d, the claim %s: %s"
              % (name, max(abs(x) for x in F), "fails" if not g1 else "HOLDS (bad)", "refused" if not g1 else "ACCEPTED"))
        ok &= not g1
        up_ok = all(not claim_holds(len(files[nm][0]), FILES[nm][4], files[nm][2], files[nm][3] + Fr(1, 10 ** 12))
                    for nm in FILES)
        print("  2. every claim raised by 10^-12 (the claims are 12-digit truncations): refused for all six: %s" % up_ok)
        ok &= up_ok

        def add_deg(cf, d):
            cf = dict(cf)
            cf[(1 << (d + 1)) - 1] = 1                      # a term of degree d + 1
            return cf, d
        g3, _ = check_file(os.path.join(root, name), mutate=add_deg, log=lambda s: None)
        print("  3. a term of degree d + 1 added: %s" % ("refused" if not g3 else "ACCEPTED (bad)"))
        ok &= not g3
        coeffs, clash = build_fd(6, parity_ok=False)
        F = wht_values(18, coeffs)
        vals = sorted(set(F))
        bad = vals != [-64, 0, 64]
        print("  4. d = 6 with c = 0 (d + |c| even): %d frequencies collide and F takes %d values (not only 0, +-64): %s"
              % (clash, len(vals), "refused" if bad else "ACCEPTED (bad)"))
        ok &= bad
    print("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("run time %.0f s" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
