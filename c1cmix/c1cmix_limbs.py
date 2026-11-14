#!/usr/bin/env python3
"""1c, the remark of PR #354 on test-measure upper bounds: program B (16-bit limbs, exact 64-bit correlations).

Written separately from program A. Each weight is split into 16-bit limbs, w = sum_p L_p 2^(16 p). Every limb
correlation c_pq[k] = sum_j L_p[j] L_q[j + k] is at most n (2^16 - 1)^2 < 2^45 for n < 8192 cells, so numpy's
int64 arithmetic is exact, and a_k = sum_pq c_pq[k] 2^(16 (p + q)) is recombined in Python integers. The value
V = min_k sum_i lam_i K a_k(i) / S_i^2 is then formed with integers only (a common denominator), and the fraction is
reduced at the end with math.gcd, so no Fraction arithmetic is shared with program A.

Needs numpy.   usage: c1cmix_limbs.py CERT.json [--bump COMPONENT INDEX] [--swap-lam]
"""
import json
import sys
from math import gcd

import numpy as np


def frac(text):
    num, _, den = text.partition("/")
    return int(num), int(den or "1")


def autocorr(row, K):
    n = len(row)
    assert n < 8192
    limbs = []
    rest = list(row)
    while any(rest):
        limbs.append(np.array([x & 0xFFFF for x in rest], dtype=np.int64))
        rest = [x >> 16 for x in rest]
    out = [0] * (K + 1)
    for p, Lp in enumerate(limbs):
        for q, Lq in enumerate(limbs):
            # full[n - 1 + k] = sum_j Lp[j] Lq[j + k]
            full = np.correlate(Lq, Lp, mode="full")
            seg = full[n - 1 : n + K].tolist()
            sh = 16 * (p + q)
            for k in range(K + 1):
                out[k] += seg[k] << sh
    return out


def main(argv):
    doc = json.loads(open(argv[1], "rb").read().decode("utf-8"))
    K = doc["K"]
    lam = [frac(x) for x in doc["lam"]]
    w = [[int(x) for x in row] for row in doc["w"]]
    for row_in, row in zip(doc["w"], w):
        if any(type(x) is not int for x in row_in):
            raise SystemExit("B: a weight is not an integer")
    if "--bump" in argv:
        j = argv.index("--bump")
        w[int(argv[j + 1])][int(argv[j + 2])] += 1
    if "--swap-lam" in argv:
        lam = lam[::-1]
    (n1, d1), (n2, d2) = lam
    good = n1 > 0 and n2 > 0 and d1 > 0 and d2 > 0 and n1 * d2 + n2 * d1 == d1 * d2
    good = good and min(min(r) for r in w) >= 0 and all(len(r) >= K + 1 for r in w)
    print("B: K = %d, lam sums to 1 and is positive, weights non-negative, at least K + 1 cells: %s" % (K, good))
    S = [sum(r) for r in w]
    a = [autocorr(r, K) for r in w]
    spot = all(a[i][k] == sum(w[i][j] * w[i][j + k] for j in range(len(w[i]) - k)) for i in (0, 1) for k in (0, 1, K))
    print("B: S = %s; lags 0, 1 and K equal the direct sums: %s" % (S, spot))
    # V_k = (n1 K a1 d2 S2^2 + n2 K a2 d1 S1^2) / (d1 d2 S1^2 S2^2)
    den = d1 * d2 * S[0] ** 2 * S[1] ** 2
    c1 = n1 * K * d2 * S[1] ** 2
    c2 = n2 * K * d1 * S[0] ** 2
    nums = [c1 * a[0][k] + c2 * a[1][k] for k in range(K + 1)]
    low = min(nums)
    at = [k for k, x in enumerate(nums) if x == low]
    g = gcd(low, den)
    print("B: V = %d / %d" % (low // g, den // g))
    print("B: V = 0.%030d..., minimum at k = %s" % (low * 10 ** 30 // den, at))
    print("B: VERDICT %s" % ("PASS" if good and spot else "REFUSED"))


if __name__ == "__main__":
    main(sys.argv)
