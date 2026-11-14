#!/usr/bin/env python3
"""88a, program B: a second, separately written check of the admissible tuples, by a different route.

Program A lists the occupied residue classes. This program instead produces a witness and checks it with gcds:

  * for each prime p <= k the occupied classes are collected in one integer bitmask; if the mask is full the tuple is
    refused; otherwise the lowest unset bit a_p is a free class;
  * the Chinese remainder theorem gives n with n = -a_p (mod p) for every prime p <= k;
  * then gcd(n + h, P) = 1 is verified for every entry h, with P the product of the primes p <= k. That equality says
    directly that no n + h is divisible by a prime p <= k, i.e. that the class -n mod p is free for each such p.
    (For p > k nothing is needed: k entries cannot occupy p classes.)

For the [Zha14] row (the k0 primes after k0) the route is also different from program A: the primes come from an
odd-only numpy sieve, and the freedom of class 0 mod p is verified directly, for EVERY prime p <= k0, by striking the
multiples of each such p across the interval [first entry, last entry] and confirming that no entry is struck.

Tuple rows: standard library only (python3 -I). The zhang mode needs numpy.

usage:  c88a_crt.py list|lean|json FILE K D [--swap OLD NEW]
        c88a_crt.py zhang K0 LIMIT [--from-first]
"""
import sys
from math import gcd


def load(kind, path):
    raw = open(path, "rb").read().decode("utf-8")
    if kind == "list":
        return [int(x) for x in raw.split()]
    if kind == "lean":
        key = "def admissibleTuple"
        if raw.count(key) != 1:
            raise SystemExit("B: admissibleTuple is not defined exactly once")
        tail = raw[raw.index(key) :]
        body = tail[tail.index("{") + 1 : tail.index("}")]
        return [int(x) for x in body.split(",")]
    if kind == "json":
        import json

        tup = json.loads(raw)["admissible_tuple"]["tuple"]
        for v in tup:
            if not isinstance(v, int) or isinstance(v, bool):
                raise SystemExit("B: non-integer entry")
        return tup[:]
    raise SystemExit("B: unknown kind")


def is_prime(m):
    if m < 2:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            return False
        d += 1
    return True


def tuple_mode(kind, path, k, diam, swap):
    H = load(kind, path)
    if swap:
        old, new = swap
        pos = [i for i, h in enumerate(H) if h == old]
        if len(pos) != 1 or new in H:
            raise SystemExit("B: bad --swap")
        H[pos[0]] = new
    good = True
    H2 = sorted(H)
    strictly = all(H2[i] < H2[i + 1] for i in range(len(H2) - 1))
    print("B: size %d (k = %d), strictly increasing after sorting: %s" % (len(H2), k, strictly))
    good = good and strictly and len(H2) == k
    print("B: diameter %d (claimed %d)" % (H2[-1] - H2[0], diam))
    good = good and (H2[-1] - H2[0] == diam)
    n, modulus, full = 0, 1, []
    for p in range(2, k + 1):
        if not is_prime(p):
            continue
        mask = 0
        for h in H2:
            mask |= 1 << (h % p)
        if mask == (1 << p) - 1:
            full.append(p)
            continue
        a = 0
        while (mask >> a) & 1:
            a += 1
        # solve n = -a (mod p) together with the congruences so far
        target = (-a) % p
        t = ((target - n) * pow(modulus, -1, p)) % p
        n += modulus * t
        modulus *= p
    if full:
        print("B: no free class modulo %s" % full)
        print("B: admissible: False")
        good = False
    else:
        worst = max(gcd(n + h, modulus) for h in H2)
        print("B: witness n = %d (mod the product %d of the primes <= %d)" % (n, modulus, k))
        print("B: largest gcd(n + h, product) over the entries: %d" % worst)
        print("B: admissible: %s" % (worst == 1))
        good = good and worst == 1
    print("B: VERDICT %s" % ("PASS" if good else "REFUSED"))


def zhang_mode(k0, limit, from_first):
    import numpy as np

    # odd-only sieve: index i stands for 2 i + 1
    half = (limit + 1) // 2
    odd = np.ones(half, dtype=np.bool_)
    odd[0] = False
    i = 1
    while (2 * i + 1) ** 2 <= limit:
        if odd[i]:
            q = 2 * i + 1
            odd[(q * q) // 2 :: q] = False
        i += 1
    primes = np.concatenate((np.array([2], dtype=np.int64), 2 * np.nonzero(odd)[0].astype(np.int64) + 1))
    primes = primes[primes <= limit]
    pi_limit = int(primes.size)
    pi_k0 = int(np.searchsorted(primes, k0, side="right"))
    print("B: pi(%d) = %d, pi(%d) = %d" % (k0, pi_k0, limit, pi_limit))
    first = 0 if from_first else pi_k0
    H = primes[first : first + k0]
    good = H.size == k0
    lo, hi = int(H[0]), int(H[-1])
    print("B: %d entries, from %d to %d, diameter %d" % (H.size, lo, hi, hi - lo))
    good = good and (hi - lo) < limit
    if not from_first:
        ineq = (pi_limit - pi_k0) > k0
        print("B: pi(%d) - pi(%d) = %d > %d: %s" % (limit, k0, pi_limit - pi_k0, k0, ineq))
        good = good and ineq
    # strike the multiples of every prime p <= k0 on [lo, hi]; an entry that is struck occupies class 0 mod p
    struck = np.zeros(hi - lo + 1, dtype=np.bool_)
    for p in primes[:pi_k0].tolist():
        start = ((lo + p - 1) // p) * p
        struck[start - lo :: p] = True
    hit = int(np.count_nonzero(struck[H - lo]))
    print("B: entries divisible by some prime p <= k0: %d" % hit)
    if hit == 0:
        print("B: class 0 is free modulo every prime p <= k0; k0 entries cannot fill p > k0 classes")
        print("B: admissible: True")
    else:
        full = []
        for p in primes[:15].tolist():
            if int(np.unique(H % p).size) == p:
                full.append(p)
        print("B: every class occupied modulo (primes <= 47): %s" % (full if full else "none"))
        print("B: admissible: %s" % ("False" if full else "not decided by this route"))
        good = False
    print("B: VERDICT %s" % ("PASS" if good else "REFUSED"))


def main(argv):
    if len(argv) >= 4 and argv[1] == "zhang":
        zhang_mode(int(argv[2]), int(argv[3]), "--from-first" in argv)
        return 0
    if len(argv) >= 5 and argv[1] in ("list", "lean", "json"):
        swap = None
        if "--swap" in argv:
            j = argv.index("--swap")
            swap = (int(argv[j + 1]), int(argv[j + 2]))
        tuple_mode(argv[1], argv[2], int(argv[3]), int(argv[4]), swap)
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
