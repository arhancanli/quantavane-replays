#!/usr/bin/env python3
"""Second exact check of a C_1a or C_1d upper-bound certificate, by a different route from check.py (which uses
big-integer products split into positive and negative parts and a random-point identity).

Definitions (teorth/optimizationproblems, constants/1a.md and 1d.md): for f = sum_i v_i 1_{cell i} on n equal cells of
[-1/4,1/4] (width h = 1/(2n)) with integers v_i (v_i >= 0 for C_1a), f*f is linear between multiples of h with knot
values h c_m, c = v*v, so
    C_1a, C_1d <= 2n max_m |c_m| / S^2,   S = sum_i v_i.
Route: every c_m modulo six NTT primes below 2^31 (numpy int64 transforms), Garner's algorithm for the mixed-radix
digits, the signed value recovered from the symmetric residue range (2 n max|v|^2 < product of the primes is checked),
max |c_m| taken over exact integers, and the deciding knot summed again directly.  A self-test against brute force
on small random inputs (signed and non-negative) runs first.
Usage: python3 ntt_check_conv.py CERT.json [--fraction=N/D] [--claim=DECIMAL] [--controls]
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


def conv_exact(v):
    """c_m = sum_i v_i v_{m-i}, m = 0..2n-2, as exact Python integers (numpy object array)."""
    n = len(v)
    L = 1 << (2 * n - 1).bit_length()
    rev = bitrev(L)
    res = []
    for p, g in PRIMES:
        x = np.zeros(L, dtype=np.int64)
        x[:n] = np.array([t % p for t in v], dtype=np.int64)
        X = ntt(x, p, g, rev)
        res.append(ntt(X * X % p, p, g, rev, invert=True)[:2 * n - 1])
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
    need(2 * n * max(abs(t) for t in v) ** 2 < prod, 'moduli too small for these values')
    c = np.zeros(2 * n - 1, dtype=object)
    rad = 1
    for d, p in zip(digits, ps):
        c = c + d.astype(object) * rad
        rad *= p
    half = prod // 2
    c = np.where(c > half, c - prod, c)
    return c


def bound(v, nonneg):
    n = len(v)
    need(all(isinstance(t, int) for t in v), 'values must be integers')
    if nonneg:
        need(all(t >= 0 for t in v), 'C_1a needs non-negative values')
    S = sum(v)
    need(S > 0, 'the integral must be positive')
    c = conv_exact(v)
    a = np.abs(c)
    cmax = max(a)
    knots = [int(m) for m in np.nonzero(a == cmax)[0]]
    m = knots[0]
    direct = sum(v[i] * v[m - i] for i in range(max(0, m - n + 1), min(m, n - 1) + 1))
    need(abs(direct) == cmax, 'direct sum at the deciding knot disagrees with the transforms')
    return Fraction(2 * n * cmax, S * S), knots, int(c[m]), S


def selftest():
    rng = np.random.default_rng(20261110)
    for trial in range(16):
        n = int(rng.integers(2, 40))
        nonneg = bool(trial % 2)
        top = int(rng.choice([3, 1000, 2 ** 40, 2 ** 56]))
        lo = 0 if nonneg else -top
        v = [int(t) for t in rng.integers(lo, top, size=n)]
        if sum(v) <= 0:
            v[0] = abs(v[0]) + 1 + abs(sum(v))
        brute = [sum(v[i] * v[m - i] for i in range(n) if 0 <= m - i < n) for m in range(2 * n - 1)]
        b, knots, cm, S = bound(v, nonneg)
        need(b == Fraction(2 * n * max(abs(t) for t in brute), sum(v) ** 2), 'self-test: bound differs from brute force')
        need(abs(brute[knots[0]]) == max(abs(t) for t in brute), 'self-test: knot differs')
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
    v = [int(t) for t in d['v_int']]
    nonneg = d.get('constant') == 'C_1a'
    kind = '1a' if nonneg else '1d'
    bd, knots, cm, S = bound(v, nonneg)
    failed = False
    print(f'{path}: C_{kind}, n = {len(v)}, S = {S}, max |v| = {max(abs(t) for t in v)}')
    print(f'  max |c_m| at m = {knots[:6]} (knots with that value: {len(knots)}), c_m = {cm}')
    print(f'  bound = {bd.numerator}/{bd.denominator}')
    print(f'  decimal (floor to 22 places): {dec(bd, 22)}')
    if '--fraction' in opts:
        ok = Fraction(opts['--fraction']) == bd
        print('  stated fraction:', 'EQUAL' if ok else 'DIFFERS'); failed |= not ok
    if '--claim' in opts:
        ok = bd <= Fraction(opts['--claim'])
        print(f'  claim C_{kind} <= {opts["--claim"]}:', 'implied' if ok else 'NOT implied'); failed |= not ok
    if '--controls' in sys.argv:
        n = len(v)
        m = knots[0]
        sgn = 1 if cm > 0 else -1
        # move v_j by one unit, in the direction sgn(c_m) sgn(v_{m-j}), where |v_{m-j}| is largest: c_m then moves by
        # 2 |v_{m-j}| away from zero, so |c_m| and hence the maximum grow
        js = [j for j in range(max(0, m - n + 1), min(m, n - 1) + 1) if 2 * j != m and v[m - j] != 0]
        j = max(js, key=lambda t: abs(v[m - t]))
        step = sgn * (1 if v[m - j] > 0 else -1)
        bad = list(v); bad[j] += step
        b2 = bound(bad, nonneg)[0]
        r1 = b2 > bd
        print(f'  control: value {j} moved by {step} (|c_m| at m = {m} grows) raises the bound: {r1}'); failed |= not r1
        if '--claim' in opts:
            c = Fraction(opts['--claim'])
            dd = len(opts['--claim'].split('.')[1])
            r3 = not (bd <= c - Fraction(1, 10 ** dd))
            print(f'  control: the claim moved one unit down ({dec(c - Fraction(1, 10 ** dd), dd)}) is refused: {r3}')
            failed |= not r3
    print(f'  run time {time.time() - t0:.0f} s')
    print('RESULT:', 'FAIL' if failed else 'PASS')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
