#!/usr/bin/env python3
"""Exact check of a step-function certificate for an upper bound on C_1b (Erdos minimum overlap constant).

Definition (teorth/optimizationproblems, constants/1b.md): C_1b is the largest constant with
    sup_{x in [-2,2]} int_{-1}^{1} f(t) g(x+t) dt >= C_1b
for all f, g: [-1,1] -> [0,1] with f + g = 1 on [-1,1] and int f = 1 (f, g extended by zero).  Any admissible f
therefore gives C_1b <= sup_x int f(t) g(x+t) dt.

Certificate: n equal cells of [-1,1] (width 2/n), f = a_i / D on cell i with integers 0 <= a_i <= D and
sum_i a_i = n D / 2 (so int f = 1 exactly); g = (D - a_i) / D = b_i / D.  For two step functions on the same grid
the overlap int f(t) g(x+t) dt is linear in x between multiples of 2/n, so its supremum is a maximum over the grid
shifts x = 2k/n, k = -(n-1)..(n-1) (k = +-n give 0):
    C_1b <= (2/n) max_k X_k / D^2,   X_k = sum_i a_i b_{i+k}   (b = 0 outside 0..n-1).

Method (Python standard library only; every number an integer or a Fraction):
  1. every X_k at once from one big-integer product (Kronecker substitution of A(x) and the reversed B(x));
  2. the whole sequence confirmed by its generating-polynomial identity at random points modulo 2^61 - 1;
  3. the deciding shift summed again directly;
  4. the bound is the reduced fraction; it must equal the fraction stored in the certificate, imply the stated
     decimal claim (bound <= claim) and beat every --record strictly (bound < record).
Red controls (--controls), each must be refused: one unit of mass moved between the two cells where the move raises
every deciding X_k most (the old fraction stated); 1% of a cell moved between the same cells (its own exact
fraction stated, so the claim must fail); the mass broken (one height raised by one unit with no compensation); a
height above D; the claim moved one unit of its last digit downwards.

Usage: python3 -I check_1b.py CERT.json [--record DECIMAL_OR_FRACTION ...] [--controls]
"""
import json
import random
import sys
import time
from fractions import Fraction

P61 = (1 << 61) - 1


def slot_bits(bound):
    return (bound.bit_length() + 2 + 7) // 8 * 8


def pack(seq, slot):
    nb = slot // 8
    return int.from_bytes(b"".join(int(x).to_bytes(nb, "little") for x in seq), "little")


def unpack(big, slot, count):
    nb = slot // 8
    raw = big.to_bytes(nb * (count + 1), "little")
    return [int.from_bytes(raw[i * nb:(i + 1) * nb], "little") for i in range(count)]


def cross_all(a, b):
    """X_k = sum_i a_i b_{i+k} for k = -(n-1)..(n-1); returned list index k + n - 1."""
    n = len(a)
    slot = slot_bits(n * max(max(a), 1) * max(max(b), 1))
    prod = pack(a, slot) * pack(b[::-1], slot)          # coefficient e = i + (n-1-j): X_{j-i} with j-i = n-1-e
    coef = unpack(prod, slot, 2 * n - 1)
    return coef[::-1]                                      # index k + n - 1 = 2n-2-e  <=>  k = n-1-e


def horner(seq, r):
    acc = 0
    for x in reversed(seq):
        acc = (acc * r + x) % P61
    return acc


def identity_ok(a, b, X, rounds=3):
    """sum_k X_k r^(k+n-1) == A(r) * Brev(r) mod P61 at random r."""
    n = len(a)
    rng = random.SystemRandom()
    brev = b[::-1]
    Xrev = X[::-1]
    for _ in range(rounds):
        r = rng.randrange(2, P61 - 1)
        if horner(Xrev, r) != horner(a, r) * horner(brev, r) % P61:
            return False
    return True


def direct(a, b, k):
    n = len(a)
    return sum(a[i] * b[i + k] for i in range(max(0, -k), min(n, n - k)))


def evaluate(a, D):
    n = len(a)
    if n < 2 or any((not isinstance(x, int)) or x < 0 or x > D for x in a):
        return None, 'heights must be integers in [0, D]'
    if 2 * sum(a) != n * D:
        return None, f'mass: 2 sum(a) = {2 * sum(a)} differs from n D = {n * D}'
    b = [D - x for x in a]
    X = cross_all(a, b)
    if not identity_ok(a, b, X):
        return None, 'generating-polynomial identity failed'
    mx = max(X)
    ks = [i - (n - 1) for i, x in enumerate(X) if x == mx]
    for k in ks[:3]:
        if direct(a, b, k) != mx:
            return None, f'direct sum at k = {k} disagrees'
    bound = Fraction(2 * mx, n * D * D)
    return {'bound': bound, 'max': mx, 'ks': ks, 'X': X}, None


def parse_num(s):
    s = s.strip()
    return Fraction(s) if '/' in s else Fraction(s)


def judge(cert, records, quiet=False):
    a, D = cert['a'], int(cert['D'])
    n = int(cert['n'])
    if len(a) != n:
        return False, f'length {len(a)} != n {n}', None
    res, err = evaluate(a, D)
    if err:
        return False, err, None
    bound = res['bound']
    stated = Fraction(int(cert['score_num']), int(cert['score_den']))
    if bound != stated:
        return False, f'fraction {bound} differs from the stated {stated}', res
    claim = parse_num(cert['claim'])
    if not bound <= claim:
        return False, f'bound {float(bound)!r} does not imply the claim {cert["claim"]}', res
    for r in records:
        if not bound < r:
            return False, f'bound does not beat the record {r} ({float(r)!r})', res
    return True, 'PASS', res


def decimal(fr, digits=22):
    q, r = divmod(fr.numerator * 10 ** digits, fr.denominator)
    s = str(q).rjust(digits + 1, '0')
    return s[:-digits] + '.' + s[-digits:]


def main():
    args = sys.argv[1:]
    path = args[0]
    records = []
    controls = False
    i = 1
    while i < len(args):
        if args[i] == '--record':
            records.append(parse_num(args[i + 1])); i += 2
        elif args[i] == '--controls':
            controls = True; i += 1
        else:
            raise SystemExit(f'unknown argument {args[i]}')
    cert = json.load(open(path))
    t0 = time.time()
    ok, msg, res = judge(cert, records)
    n = int(cert['n'])
    print(f'certificate: {path}')
    print(f'n = {n} cells, D = {cert["D"]}, sum(a) = n D / 2 checked')
    if res:
        print(f'max_k X_k = {res["max"]} at k = {res["ks"][:10]} ({len(res["ks"])} shift(s))')
        print(f'bound = {res["bound"].numerator}/{res["bound"].denominator}')
        print(f'      = {decimal(res["bound"])}...')
        X = res['X']; mx = res['max']
        near = sorted(X, reverse=True)
        second = next((x for x in near if x < mx), None)
        if second is not None:
            print(f'next largest X_k is lower by {Fraction(2 * (mx - second), n * int(cert["D"]) ** 2)} in the bound')
    print(f'claim: C_1b <= {cert["claim"]}; records beaten strictly: {[str(r) for r in records]}')
    print(f'RESULT: {msg} ({time.time() - t0:.1f} s)')
    if not ok:
        sys.exit(1)
    if controls:
        D = int(cert['D']); a = list(cert['a'])
        k0 = res['ks'][0]
        B0 = res['bound']
        # Moving u units of mass from cell j to cell i changes X_k, for every deciding shift k, by
        # u (G_k(i) - G_k(j)) + O(u^2) with G_k(c) = (D - a_{c+k}) - a_{c-k} (a = 0 outside the cells).  For the
        # controls we take i where min_k G_k is largest and j where max_k G_k is smallest, so the move raises every
        # deciding X_k and hence the exact bound.
        ks = res['ks']
        def G(c, k):
            up = a[c + k] if 0 <= c + k < n else 0
            dn = a[c - k] if 0 <= c - k < n else 0
            return (D - up) - dn
        cand_i = [c for c in range(n) if a[c] + D // 100 <= D]
        cand_j = [c for c in range(n) if a[c] >= D // 100]
        i0 = max(cand_i, key=lambda c: min(G(c, k) for k in ks))
        j0 = min((c for c in cand_j if c != i0), key=lambda c: max(G(c, k) for k in ks))
        # control 1: one unit of mass moved from j0 to i0 (mass kept), the old fraction stated
        t = a[:]; t[i0] += 1; t[j0] -= 1
        c1 = dict(cert); c1['a'] = t
        # control 1b: 1% of a cell moved from j0 to i0, its own exact fraction stated: the claim must then fail
        t = a[:]; t[i0] += D // 100; t[j0] -= D // 100
        r2, _e = evaluate(t, D)
        c1b = dict(cert); c1b['a'] = t
        c1b['score_num'] = str(r2['bound'].numerator); c1b['score_den'] = str(r2['bound'].denominator)
        print(f'controls move mass from cell {j0} to cell {i0}; the 1% move gives the bound {float(r2["bound"])!r}')
        # control 2: mass broken
        c2 = dict(cert); c2['a'] = a[:]; c2['a'][i0] += 1
        # control 3: a height above D
        c3 = dict(cert); c3['a'] = a[:]; jj = next(j for j in range(n) if a[j] == D) if D in a else 0
        c3['a'][jj] = D + 1; c3['a'][next(j for j in range(n) if a[j] > 0 and j != jj)] -= 1
        # control 4: the claim one unit of its last digit lower
        cl = cert['claim']; dec = len(cl.split('.')[1])
        lower = Fraction(cl) - Fraction(1, 10 ** dec)
        c4 = dict(cert); c4['claim'] = decimal(lower, dec)
        for name, cc in [('one unit of mass moved, old fraction stated', c1),
                         ('1% of a cell moved, its own fraction stated', c1b), ('mass broken', c2),
                         ('height above D', c3), (f'claim moved to {c4["claim"]}', c4)]:
            ok2, msg2, _ = judge(cc, records)
            print(f'control [{name}]: {"PASSED (BAD)" if ok2 else "refused"}: {msg2[:160]}')
            if ok2:
                print('RED CONTROL FAILED'); sys.exit(2)
        print('CONTROLS: all refused as expected')


if __name__ == '__main__':
    main()
