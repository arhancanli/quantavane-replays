#!/usr/bin/env python3
"""Exact replay of the 26a lower-bound row [ADGP2025] (teorth/optimizationproblems, main 87b4ebb): "Degree d address
function achieves the bound 2^((d-1)/d)", so C_26a = sup_d BH^{<=d} >= sup_d 2^((d-1)/d) = 2.

The degree-d address function on n = k + 2^k bits (k = d - 1): f(a, y) = y_a, the address a in {+-1}^k selecting one of
the 2^k data bits.  For d = 1, ..., 5 (n up to 20) the full truth table is transformed exactly (fast Walsh-Hadamard,
integers): every nonzero coefficient must have |fhat(S)| = 2^-k, there must be 4^k of them, the degree must be d and
||f||_inf = 1.  Then (sum |fhat(S)|^p)^(1/p) with p = 2d/(d+1) is N^(1/p) 2^-k = 2^(2k/p - k) = 2^((d-1)/d) exactly
(the exponent computed as a rational).  Written from the page's row and the standard definition; no download.
usage: python3 -I address_check.py [--controls]
"""
import sys
from fractions import Fraction

CONTROLS = "--controls" in sys.argv
DMAX = 5


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def table(k, wiring=None):
    """Values (+1/-1) of f(a, y) = y_{wiring(a)} on {0,1}^n (bit i set means x_i = -1); bits 0..k-1 are the address."""
    n = k + (1 << k)
    out = []
    for x in range(1 << n):
        a = x & ((1 << k) - 1)
        b = wiring(a) if wiring else a
        out.append(-1 if (x >> (k + b)) & 1 else 1)
    return n, out


def wht(v):
    v = list(v)
    h = 1
    while h < len(v):
        for i in range(0, len(v), 2 * h):
            for j in range(i, i + h):
                a, b = v[j], v[j + h]
                v[j], v[j + h] = a + b, a - b
        h *= 2
    return v


def analyse(k, wiring=None):
    n, f = table(k, wiring)
    F = wht(f)
    need(sum(c * c for c in F) == 4 ** n, "Parseval fails")          # sum fhat^2 = 1 for a +-1 function
    nz = [(S, c) for S, c in enumerate(F) if c]
    mags = {Fraction(abs(c), 2 ** n) for _, c in nz}
    deg = max(bin(S).count("1") for S, _ in nz)
    return n, len(nz), mags, deg, max(abs(x) for x in f)


def main():
    ok = True
    print("d  k  n   nonzero coefficients  |fhat|     degree  ||f||_inf  (sum |fhat|^p)^(1/p), p = 2d/(d+1)")
    for d in range(1, DMAX + 1):
        k = d - 1
        n, N, mags, deg, sup = analyse(k)
        good = N == 4 ** k and mags == {Fraction(1, 2 ** k)} and deg == d and sup == 1
        p = Fraction(2 * d, d + 1)
        expo = Fraction(2 * k) / p - k          # log2 of N^(1/p) * 2^-k with N = 4^k
        good &= expo == Fraction(d - 1, d)
        ok &= good
        print(f"{d}  {k}  {n:<3} {N:<21} {'2^-' + str(k):<10} {deg:<7} {sup:<10} 2^({expo}) = {2 ** float(expo):.12f}"
              f"{'' if good else '   FAILS'}")
    print(f"RESULT [ADGP2025]: BH^(<=d) >= 2^((d-1)/d) for d = 1..{DMAX} exactly, increasing to 2, so C_26a >= 2 (not "
          f"attained): {'HOLDS' if ok else 'FAILS'}")
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        # two addresses wired to the same data bit: no longer the address function
        n, N, mags, deg, sup = analyse(2, wiring=lambda a: 0 if a == 1 else a)
        res.append(("k = 2 with addresses 0 and 1 both reading y_0", not (N == 16 and mags == {Fraction(1, 4)})))
        print(f"  k = 2, addresses 0 and 1 both read y_0: {N} nonzero coefficients, magnitudes {sorted(mags)}")
        res.append(("claimed exponent d/(d+1) at d = 3", Fraction(2 * 2) / Fraction(6, 4) - 2 != Fraction(3, 4)))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
