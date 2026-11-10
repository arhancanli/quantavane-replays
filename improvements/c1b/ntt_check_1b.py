#!/usr/bin/env python3
"""Second exact check of a C_1b (Erdos minimum overlap) upper-bound certificate, by a different route from
check_1b.py (which uses one big-integer product and a random-point identity).

Definition (teorth/optimizationproblems, constants/1b.md): for f, g: [-1,1] -> [0,1], f + g = 1 on [-1,1],
int f = 1, extended by zero, C_1b <= sup_x int f(t) g(x+t) dt.  Certificate: n equal cells, f = a_i / D,
g = (D - a_i) / D, integers 0 <= a_i <= D, sum a_i = n D / 2.  The overlap is linear in x between grid shifts, so
    C_1b <= 2 max_k X_k / (n D^2),   X_k = sum_i a_i (D - a_{i+k})  (zero outside the cells), |k| < n.

Route:
  1. every X_k modulo six NTT primes below 2^31 by number-theoretic transforms in numpy int64;
  2. Garner's algorithm gives the mixed-radix digits of each X_k (exact, since X_k <= n D^2 is checked to be below
     the product of the primes), and the largest X_k is found by comparing digit vectors (np.lexsort);
  3. that X_k is rebuilt as a Python integer and summed again directly; the bound is the reduced fraction.
A self-test against brute force on small random inputs runs first.

Usage: python3 ntt_check_1b.py CERT.json [--claim=DECIMAL] [--fraction=N/D] [--controls]
"""
import json
import sys
import time
from fractions import Fraction

import numpy as np

PRIMES = [(2013265921, 31), (1811939329, 13), (2113929217, 5), (1711276033, 29), (998244353, 3), (754974721, 11)]


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def bitrev(n):
    bits = n.bit_length() - 1
    idx = np.arange(n, dtype=np.int64)
    rev = np.zeros(n, dtype=np.int64)
    for b in range(bits):
        rev |= ((idx >> b) & 1) << (bits - 1 - b)
    return rev


def ntt(a, p, g, rev, invert=False):
    n = len(a)
    a = a[rev] % p
    root = pow(g, (p - 1) // n, p)
    if invert:
        root = pow(root, p - 2, p)
    tw = np.empty(max(n // 2, 1), dtype=np.int64)
    tw[0] = 1
    step = 1
    while step < n // 2:
        tw[step:2 * step] = tw[:step] * pow(root, step, p) % p
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


def cross_digits(a, b):
    """Garner digits of X_k = sum_i a_i b_{i+k}, k = -(n-1)..(n-1) (array index k + n - 1)."""
    n = len(a)
    L = 1 << (2 * n - 1).bit_length()
    rev = bitrev(L)
    res = []
    for p, g in PRIMES:
        x = np.zeros(L, dtype=np.int64)
        y = np.zeros(L, dtype=np.int64)
        x[:n] = np.array([v % p for v in a[::-1]], dtype=np.int64)   # u_j = a_{n-1-j}
        y[:n] = np.array([v % p for v in b], dtype=np.int64)
        c = ntt(ntt(x, p, g, rev) * ntt(y, p, g, rev) % p, p, g, rev, invert=True)
        res.append(c[:2 * n - 1])                                     # (u*b)_{n-1+k} = X_k
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


def value_at(digits, ps, idx):
    v, rad = 0, 1
    for d, p in zip(digits, ps):
        v += int(d[idx]) * rad
        rad *= p
    return v


def bound(a, D):
    n = len(a)
    need(all(isinstance(x, int) and 0 <= x <= D for x in a), 'heights must be integers in [0, D]')
    need(2 * sum(a) == n * D, 'mass: sum a must equal n D / 2')
    b = [D - x for x in a]
    digits, ps, prod = cross_digits(a, b)
    need(n * D * D < prod, 'moduli too small')
    order = np.lexsort(tuple(digits))          # ascending; the last is the largest X_k
    top = int(order[-1])
    xmax = value_at(digits, ps, top)
    ties = 1
    while ties < len(order) and all(int(d[order[-1 - ties]]) == int(d[top]) for d in digits):
        ties += 1
    k = top - (n - 1)
    direct = sum(a[i] * b[i + k] for i in range(max(0, -k), min(n, n - k)))
    need(direct == xmax, 'direct sum at the largest shift disagrees with the transforms')
    nxt = value_at(digits, ps, int(order[-1 - ties])) if ties < len(order) else None
    tied = sorted(int(order[-1 - t]) - (n - 1) for t in range(ties))
    return Fraction(2 * xmax, n * D * D), k, xmax, ties, nxt, tied


def selftest():
    rng = np.random.default_rng(20261110)
    for trial in range(16):
        n = int(rng.integers(2, 40)) * 2
        D = int(rng.choice([2, 10, 1000, 2 ** 40]))
        a = [int(x) for x in rng.integers(0, D + 1, size=n)]
        # repair the mass to n D / 2
        diff = n * D // 2 - sum(a)
        i = 0
        while diff != 0:
            j = i % n
            if diff > 0 and a[j] < D:
                s = min(D - a[j], diff); a[j] += s; diff -= s
            elif diff < 0 and a[j] > 0:
                s = min(a[j], -diff); a[j] -= s; diff += s
            i += 1
        b = [D - x for x in a]
        brute = {k: sum(a[i] * b[i + k] for i in range(n) if 0 <= i + k < n) for k in range(-(n - 1), n)}
        mx = max(brute.values())
        bd, k, xmax, ties, _, _t = bound(a, D)
        need(xmax == mx and brute[k] == mx, 'self-test: largest shift differs from brute force')
        need(bd == Fraction(2 * mx, n * D * D), 'self-test: bound differs')
    return True


def dec(fr, digits):
    q = fr.numerator * 10 ** digits // fr.denominator
    s = str(q).rjust(digits + 1, '0')
    return s[:-digits] + '.' + s[-digits:]


def main():
    path = sys.argv[1]
    opts = dict(x.split('=', 1) for x in sys.argv[2:] if '=' in x)
    t0 = time.time()
    print('self-test against brute force (16 random inputs):', selftest())
    d = json.load(open(path))
    n, D = int(d['n']), int(d['D'])
    a = [int(x) for x in d['a']]
    need(len(a) == n, 'length differs from n')
    bd, k, xmax, ties, nxt, tied = bound(a, D)
    failed = False
    print(f'{path}: n = {n}, D = {D}, mass n D / 2 checked')
    print(f'  largest X_k at k = {k} (shifts with that value: {ties}), X = {xmax}')
    print(f'  bound = {bd.numerator}/{bd.denominator}')
    print(f'  decimal (floor to 22 places): {dec(bd, 22)}')
    if nxt is not None:
        print(f'  next shift value 2 X / (n D^2) = {dec(Fraction(2 * nxt, n * D * D), 22)}')
    if '--fraction' in opts:
        ok = Fraction(opts['--fraction']) == bd
        print('  stated fraction:', 'EQUAL' if ok else 'DIFFERS'); failed |= not ok
    if '--claim' in opts:
        ok = bd <= Fraction(opts['--claim'])
        print(f'  claim C_1b <= {opts["--claim"]}:', 'implied' if ok else 'NOT implied'); failed |= not ok
    if '--controls' in sys.argv:
        # one unit moved from j0 to i0 changes every deciding X_k by g_k(i0) - g_k(j0), g_k(c) = (D - a[c+k]) - a[c-k];
        # i0 maximises min_k g_k and j0 minimises max_k g_k over the tied deciding shifts, so the bound rises
        def gk(c, kk):
            return (D - (a[c + kk] if 0 <= c + kk < n else 0)) - (a[c - kk] if 0 <= c - kk < n else 0)
        i0 = max((c for c in range(n) if a[c] < D), key=lambda c: min(gk(c, kk) for kk in tied))
        j0 = min((c for c in range(n) if a[c] > 0 and c != i0), key=lambda c: max(gk(c, kk) for kk in tied))
        bad = list(a); bad[i0] += 1; bad[j0] -= 1
        b2 = bound(bad, D)[0]
        r1 = b2 != bd
        print(f'  control: height {i0} +1 and height {j0} -1 changes the bound: {r1}'); failed |= not r1
        try:
            bad2 = list(a); bad2[i0] += 1; bound(bad2, D); r2 = False
        except AssertionError:
            r2 = True
        print(f'  control: broken mass refused: {r2}'); failed |= not r2
        if '--claim' in opts:
            c = Fraction(opts['--claim'])
            step = Fraction(1, 10 ** len(opts['--claim'].split('.')[1]))
            r3 = not (bd <= c - step)
            print(f'  control: the claim moved one unit down ({dec(c - step, len(opts["--claim"].split(".")[1]))}) is refused: {r3}')
            failed |= not r3
    print(f'  run time {time.time() - t0:.0f} s')
    print('RESULT:', 'FAIL' if failed else 'PASS')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
