#!/usr/bin/env python3
"""Second exact check of a C_1c step-function certificate, by a different route from check.py.

Definition (teorth/optimizationproblems, constants/1c.md): C_1c is the smallest constant with
min_{0<=t<=1} int f(x) f(x+t) dx <= C_1c (int f)^2 for every non-negative f in L^1(R).  For
f = sum_i w_i 1_[i/K,(i+1)/K) with integers w_i >= 0 the autocorrelation is linear between multiples of 1/K, so
C_1c >= K * min_{0<=k<=K} A_k / S^2 with A_k = sum_i w_i w_{i+k} and S = sum_i w_i.

Route (shares no code with check.py, which uses one big-integer product and a random-point identity):
  1. every A_k modulo six NTT primes below 2^31 by number-theoretic transforms in numpy int64 (all products of two
     residues stay below 2^62);
  2. Garner's algorithm turns the six residues of each A_k into mixed-radix digits; since every A_k is below
     n * max(w)^2, which is checked to be below the product of the primes, the digits determine A_k exactly, and the
     least A_k over 0 <= k <= K is found by comparing digit vectors lexicographically (np.lexsort);
  3. that A_k is rebuilt as a Python integer and summed again directly; the bound is the reduced fraction.
A self-test against the definition (brute force) on small random inputs runs first.

Usage: python3 ntt_check.py CERT.json [--claim DECIMAL] [--fraction N/D] [--controls]
Exit status 0 only if the self-test passes, the fraction equals --fraction (when given), --claim is implied, and
every red control is refused as it should be.
"""
import json
import sys
import time
from fractions import Fraction

import numpy as np

# (prime, primitive root); each prime is c * 2^k + 1 with k >= 23, so transforms of length up to 2^23 exist
PRIMES = [(2013265921, 31), (1811939329, 13), (2113929217, 5), (1711276033, 29), (998244353, 3), (754974721, 11)]


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def ntt(a, p, g, invert=False):
    """In-place style iterative radix-2 transform of the int64 array a (length a power of two) modulo p."""
    n = len(a)
    bits = n.bit_length() - 1
    idx = np.arange(n, dtype=np.int64)
    rev = np.zeros(n, dtype=np.int64)
    for b in range(bits):
        rev |= ((idx >> b) & 1) << (bits - 1 - b)
    a = a[rev] % p
    root = pow(g, (p - 1) // n, p)
    if invert:
        root = pow(root, p - 2, p)
    tw = np.empty(n // 2, dtype=np.int64)  # root^j for j < n/2
    tw[0] = 1
    step = 1
    while step < n // 2:  # doubling: tw[step:2 step] = tw[:step] * root^step
        rs = pow(root, step, p)
        tw[step:2 * step] = tw[:step] * rs % p
        step *= 2
    length = 2
    while length <= n:
        half = length // 2
        w = tw[::n // length][:half]
        blk = a.reshape(-1, length)
        u = blk[:, :half].copy()
        v = blk[:, half:] * w % p
        blk[:, :half] = (u + v) % p
        blk[:, half:] = (u - v) % p
        a = blk.reshape(-1)
        length *= 2
    if invert:
        a = a * pow(n, p - 2, p) % p
    return a


def autocorr_mod(w_res, n, p, g):
    """A_k mod p for k = 0..n-1, from the residues w_res (int64) of the weights."""
    L = 1 << (2 * n - 1).bit_length()
    x = np.zeros(L, dtype=np.int64)
    y = np.zeros(L, dtype=np.int64)
    x[:n] = w_res[::-1]   # u_j = w_{n-1-j}
    y[:n] = w_res
    c = ntt(ntt(x, p, g) * ntt(y, p, g) % p, p, g, invert=True)
    return c[n - 1:2 * n - 1]  # (u * w)_{n-1+k} = sum_i w_i w_{i+k}


def all_lags(w, upto):
    """Mixed-radix digits (Garner) of A_k for k < upto, the moduli, and their product."""
    n = len(w)
    res = []
    for p, g in PRIMES:
        wr = np.array([x % p for x in w], dtype=np.int64)
        res.append(autocorr_mod(wr, n, p, g)[:upto])
    ps = [p for p, _ in PRIMES]
    digits = []
    for j, p in enumerate(ps):
        t = res[j].copy()
        for i in range(j):
            t = (t - digits[i]) % p * pow(ps[i], p - 2, p) % p
        digits.append(t)
    prod = 1
    for p in ps:
        prod *= p
    return digits, ps, prod


def from_digits(digits, ps, k):
    val, rad = 0, 1
    for d, p in zip(digits, ps):
        val += int(d[k]) * rad
        rad *= p
    return val


def bound(K, w):
    n = len(w)
    need(n > K, "need more than K weights, so that every A_k with k <= K is a full sum")
    need(all(isinstance(x, int) and x >= 0 for x in w) and any(w), "weights must be non-negative integers, not all 0")
    upto = K + 1
    digits, ps, prod = all_lags(w, upto)
    need(n * max(w) ** 2 < prod, "the moduli are too small for these weights")
    order = np.lexsort(tuple(digits))  # last key (the top digit) is the primary key
    kmin = int(order[0])
    amin = from_digits(digits, ps, kmin)
    ties = 1
    while ties < upto and all(int(d[order[ties]]) == int(d[kmin]) for d in digits):
        ties += 1
    direct = sum(w[i] * w[i + kmin] for i in range(n - kmin))
    need(direct == amin, "the direct sum at the least lag disagrees with the transforms")
    nxt = from_digits(digits, ps, int(order[ties])) if ties < upto else None
    S = sum(w)
    return Fraction(K * amin, S * S), kmin, amin, S, ties, nxt


def selftest():
    rng = np.random.default_rng(20261109)
    for trial in range(12):
        n = int(rng.integers(3, 40))
        K = int(rng.integers(1, n))
        top = int(rng.choice([3, 1000, 2 ** 40, 2 ** 60]))
        w = [int(x) for x in rng.integers(0, top, size=n)]
        w[0] = max(w[0], 1)
        brute = [sum(w[i] * w[i + k] for i in range(n - k)) for k in range(K + 1)]
        b, kmin, amin, S, ties, _ = bound(K, w)
        need(amin == min(brute) and kmin == brute.index(min(brute)), "self-test: least lag differs from brute force")
        need(b == Fraction(K * min(brute), sum(w) ** 2), "self-test: bound differs from brute force")
    return True


def dec(fr, digits):
    q = fr.numerator * 10 ** digits // fr.denominator
    s = str(q).rjust(digits + 1, "0")
    return s[:-digits] + "." + s[-digits:]


def main():
    path = sys.argv[1]
    opts = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
    t0 = time.time()
    print("self-test against brute force (12 random inputs):", selftest())
    d = json.load(open(path))
    K = int(d["K"])
    w = [int(x) for x in d["weights_int"]]
    b, kmin, amin, S, ties, nxt = bound(K, w)
    failed = False
    print(f"{path}: K = {K}, n = {len(w)}, S = {S}, max weight = {max(w)}")
    print(f"  least A_k at k = {kmin} (lags with that value: {ties}), A = {amin}")
    print(f"  bound = {b.numerator}/{b.denominator}")
    print(f"  decimal (floor): {dec(b, 22)}")
    if nxt is not None:
        print(f"  next lag value K A / S^2 = {dec(Fraction(K * nxt, S * S), 22)}")
    if "--fraction" in opts:
        ok = Fraction(opts["--fraction"]) == b
        print("  stated fraction:", "EQUAL" if ok else "DIFFERS")
        failed |= not ok
    if "--claim" in opts:
        ok = Fraction(opts["--claim"]) <= b
        print(f"  claim {opts['--claim']}:", "implied" if ok else "NOT implied")
        failed |= not ok
    if "--controls" in sys.argv:
        i = max(range(len(w)), key=lambda j: (w[j] > 0, j % 7919))  # some non-zero weight
        bad = list(w)
        bad[i] += 1
        b2 = bound(K, bad)[0]
        r1 = b2 != b
        print(f"  control: weight {i} raised by 1 changes the bound: {r1}")
        if "--claim" in opts:
            c = Fraction(opts["--claim"])
            step = Fraction(1, 10 ** len(opts["--claim"].split(".")[1]))
            r2 = not (c + step <= b)
            print(f"  control: the claim moved one unit up ({c + step}) is refused: {r2}")
            failed |= not r2
        if "--fraction" in opts:
            r3 = Fraction(opts["--fraction"]) != b2
            print(f"  control: the changed certificate fails the stated fraction: {r3}")
            failed |= not r3
        failed |= not r1
    print(f"  run time {time.time() - t0:.0f} s")
    print("RESULT:", "FAIL" if failed else "PASS")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
