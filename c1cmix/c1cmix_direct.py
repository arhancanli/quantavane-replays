#!/usr/bin/env python3
"""1c, the remark of PR #354 on test-measure upper bounds: program A (direct sums).

The certificate gives two step functions f_1, f_2 on K equal steps of width 1/K (non-negative integer weights,
K + 1 cells) and rationals lam_1 + lam_2 = 1. For a step function the autocorrelation t -> int f(x) f(x + t) dx is
linear between multiples of 1/K, and at t = k/K the normalised value is K a_k / S^2, with a_k = sum_j w_j w_(j+k)
and S = sum_j w_j. So

    V = min over 0 <= t <= 1 of sum_i lam_i (f_i * f_i)(t) / (int f_i)^2 = min over 0 <= k <= K of
        sum_i lam_i K a_k(i) / S_i^2,

and for every probability measure nu on [0, 1], sup_f int (f * f) d nu >= V (f >= 0, int f = 1).

This program computes every a_k by the direct double sum in Python integers and V in exact rationals.
Standard library only (python3 -I).

usage: c1cmix_direct.py CERT.json [--bump COMPONENT INDEX] [--swap-lam]
"""
import json
import sys
from fractions import Fraction


def main(argv):
    with open(argv[1], encoding="utf-8") as fh:
        doc = json.load(fh)
    K = doc["K"]
    lam = [Fraction(x) for x in doc["lam"]]
    w = [list(row) for row in doc["w"]]
    if "--bump" in argv:
        j = argv.index("--bump")
        w[int(argv[j + 1])][int(argv[j + 2])] += 1
    if "--swap-lam" in argv:
        lam.reverse()
    ok = type(K) is int and K > 0 and len(w) == len(lam) == 2
    ok = ok and all(type(x) is int and x >= 0 for row in w for x in row) and all(len(row) > K for row in w)
    ok = ok and all(x > 0 for x in lam) and sum(lam) == 1
    print("A: K = %d, cells %s, weights non-negative integers, lam = %s, sum 1 and positive: %s" % (
        K, [len(row) for row in w], [str(x) for x in lam], ok))
    mix = [Fraction(0)] * (K + 1)
    for i, row in enumerate(w):
        n = len(row)
        S = sum(row)
        scale = lam[i] * K / (S * S)
        print("A: component %d: S = %d" % (i + 1, S))
        for k in range(K + 1):
            a = 0
            for j in range(n - k):
                a += row[j] * row[j + k]
            mix[k] += scale * a
    V = min(mix)
    where = [k for k in range(K + 1) if mix[k] == V]
    print("A: V = %d / %d" % (V.numerator, V.denominator))
    print("A: V = 0.%030d..., minimum at k = %s" % (V.numerator * 10 ** 30 // V.denominator, where))
    print("A: VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main(sys.argv)
