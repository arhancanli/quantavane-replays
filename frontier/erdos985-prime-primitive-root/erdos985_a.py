#!/usr/bin/env python3
"""Program A (Python, exact integers) for the finite part of Erdos problem 985 as claimed by
Howe, Huffman, Kohlmeyer, Kruer (github.com/H0W4RD-dev/erdos-985-prime-primitive-root @ 99fd0a7).

Claim checked: for every odd prime p in the stated finite ranges there is a prime q < p of
multiplicative order p - 1 modulo p.

Method A (written without reading the claimants' verifier):
  * primes by a bytearray sieve; p - 1 factored with a smallest-prime-factor table;
  * q is a primitive root mod p  iff  q^((p-1)/r) != 1 (mod p) for every prime r | p - 1;
  * (a) every odd prime p <= 2791: least prime primitive root g(p), g(p) < p; the table's q is
    prime, q < p and a primitive root (also its factorisation column and residues column);
  * (b) companion certificate: the candidate moduli (odd n, 2792 <= n < u_k^2, omega(n-1) = k,
    k = 3..9) are enumerated by a depth-first search over factorisations of n - 1; n is tested
    for primality by trial division; for every prime n the least prime primitive root is
    computed and must be < n; the certificate's word for n must be type 2 with a prime
    primitive root b < n; for composite n the word must not be type 2;
  * (c) extension: every odd prime p <= BOUND has a prime primitive root < p.

Usage: erdos985_a.py SMALL_TSV STRIP_DIR [--bound N] [--control witness|strip]
Exit 0 and 'A: VERDICT PASS' when every check holds; 'A: VERDICT REFUSED' otherwise.
"""
import sys
import os
from array import array

UPPER = {3: 143 ** 2, 4: 363 ** 2, 5: 871 ** 2, 6: 2019 ** 2, 7: 4321 ** 2, 8: 9693 ** 2, 9: 20443 ** 2}
LOWER = 2792


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    i = 2
    while i * i <= n:
        if s[i]:
            s[i * i::i] = bytes(len(range(i * i, n + 1, i)))
        i += 1
    return s


def spf_table(n, primes_small):
    spf = array('I', [0]) * (n + 1)
    # largest primes first, so the smallest prime factor is written last
    for p in reversed(primes_small):
        cnt = len(range(p, n + 1, p))
        spf[p::p] = array('I', [p]) * cnt
    return spf


def distinct_primes_spf(m, spf):
    out = []
    while m > 1:
        r = spf[m]
        if r == 0:
            r = m
        out.append(r)
        while m % r == 0:
            m //= r
    return out


def is_prim_root(q, p, rs):
    if q % p == 0:
        return False
    for r in rs:
        if pow(q, (p - 1) // r, p) == 1:
            return False
    return True


def least_prime_root(p, rs, plist):
    for q in plist:
        if q >= p:
            return None
        if is_prim_root(q, p, rs):
            return q
    return None


def is_prime_td(n, plist):
    if n < 2:
        return False
    for r in plist:
        if r * r > n:
            return True
        if n % r == 0:
            return n == r
    raise ValueError("prime list too short for %d" % n)


def enumerate_candidates(k, lo, hi, P):
    """all n with lo <= n < hi, n - 1 even with exactly k distinct prime factors; returns sorted [(n, primes)]"""
    out = []
    lim = hi - 2  # n - 1 <= hi - 2

    def rec(N, t, ps):
        need = k - len(ps)
        if need == 0:
            if N + 1 >= lo:
                out.append((N + 1, tuple(ps)))
            return
        while t + need <= len(P):
            m = 1
            for x in P[t:t + need]:
                m *= x
            if N * m > lim:
                return
            p = P[t]
            pe = p
            while N * pe <= lim:
                rec(N * pe, t + 1, ps + [p])
                pe *= p
            t += 1
        raise ValueError("prime list too short")

    e2 = 2
    while e2 <= lim:
        rec(e2, 1, [2])
        e2 *= 2
    out.sort()
    return out


def main():
    args = sys.argv[1:]
    bound = 10 ** 7
    control = None
    if "--bound" in args:
        i = args.index("--bound"); bound = int(args[i + 1]); del args[i:i + 2]
    if "--control" in args:
        i = args.index("--control"); control = args[i + 1]; del args[i:i + 2]
    tsv, stripdir = args[0], args[1]
    fails = []

    top = max(bound, 50000)
    s = sieve(top)
    plist = [i for i in range(2, 50001) if s[i]]
    sbound = max(bound, 2791)
    spf = spf_table(sbound, [p for p in plist if p * p <= sbound])

    # (a) the small table
    rows = []
    with open(tsv) as f:
        head = f.readline().rstrip("\n").split("\t")
        assert head[:2] == ["p", "q"], head
        for line in f:
            if line.strip():
                rows.append(line.rstrip("\n").split("\t"))
    if control == "witness":
        for r in rows:
            if r[0] == "7":
                r[1] = "2"  # 2 has order 3 modulo 7
        print("A: CONTROL witness: table row p=7 given q=2 (order 3 mod 7)")
    small = [p for p in range(3, 2792) if s[p]]
    tab = {}
    for r in rows:
        p = int(r[0])
        if p in tab:
            fails.append("duplicate row p=%d" % p)
        tab[p] = r
    if sorted(tab) != small:
        fails.append("table rows are not exactly the odd primes <= 2791")
    nleast = 0
    gmax = (0, 0)
    for p in small:
        rs = distinct_primes_spf(p - 1, spf)
        g = least_prime_root(p, rs, plist)
        if g is None or g >= p:
            fails.append("(a) no prime primitive root below p=%d" % p)
            continue
        gmax = max(gmax, (g, p))
        r = tab.get(p)
        if r is None:
            continue
        q = int(r[1])
        if not (q < p and s[q] and is_prim_root(q, p, rs)):
            fails.append("(a) table row p=%d: q=%d is not a prime primitive root below p" % (p, q))
        if q == g:
            nleast += 1
        # extra: factorisation column and residue column
        fac = {}
        for term in r[2].split("*"):
            b, e = term.split("^")
            fac[int(b)] = int(e)
        prod = 1
        for b, e in fac.items():
            prod *= b ** e
        if prod != p - 1 or sorted(fac) != sorted(rs):
            fails.append("(a) factorisation column wrong at p=%d" % p)
        res = dict(tuple(map(int, t.split(":"))) for t in r[3].split(";"))
        for rr in rs:
            if res.get(rr) != pow(q, (p - 1) // rr, p):
                fails.append("(a) residue column wrong at p=%d r=%d" % (p, rr))
    print("A (a): odd primes p <= 2791: %d; table rows %d; least prime primitive root < p for all: %s; "
          "table q is the least one in %d rows; largest least prime primitive root %d (p = %d)"
          % (len(small), len(rows), "yes" if not [x for x in fails if x.startswith("(a) no")] else "NO",
             nleast, gmax[0], gmax[1]))

    # (b) the companion certificate
    P = plist
    tot = dict(cand=0, prime=0, w0=0, w1=0, w2=0)
    bmax = 0
    gmaxb = (0, 0)
    for k in range(3, 10):
        path = os.path.join(stripdir, "strip-%02d.u16" % (k - 3))
        raw = open(path, "rb").read()
        words = [raw[i] | (raw[i + 1] << 8) for i in range(0, len(raw), 2)]
        cands = enumerate_candidates(k, LOWER, UPPER[k], P)
        if control == "strip" and k == 3:
            for i, (n, ps) in enumerate(cands):
                if words[i] >> 14 == 2:
                    q = next(x for x in P if not is_prim_root(x, n, ps))
                    print("A: CONTROL strip: rank 3, n=%d: witness %d replaced by the prime %d (not a primitive root)"
                          % (n, words[i] & 0x3FFF, q))
                    words[i] = (2 << 14) | q
                    break
        if len(cands) != len(words):
            fails.append("(b) rank %d: %d candidates but %d words" % (k, len(cands), len(words)))
        np_ = 0
        cnt = [0, 0, 0, 0]
        for (n, ps), wd in zip(cands, words):
            a, b = wd >> 14, wd & 0x3FFF
            cnt[a] += 1
            if is_prime_td(n, P):
                np_ += 1
                g = least_prime_root(n, ps, P)
                if g is None or g >= n:
                    fails.append("(b) no prime primitive root below n=%d" % n)
                else:
                    gmaxb = max(gmaxb, (g, n))
                if a != 2 or not (1 < b < n and is_prime_td(b, P) and is_prim_root(b, n, ps)):
                    fails.append("(b) prime n=%d: certificate word (%d, %d) is not a prime primitive root" % (n, a, b))
                bmax = max(bmax, b)
            else:
                if a == 2:
                    fails.append("(b) composite n=%d has a type-2 word" % n)
                elif a == 0 and not (1 < b < n and n % b == 0):
                    fails.append("(b) composite n=%d: type-0 word %d is not a proper divisor" % (n, b))
                elif a == 1 and not (pow(b, n - 1, n) != 1 or 1 < __import__('math').gcd(pow(b, (n - 1) // 2, n) - 1, n) < n):
                    fails.append("(b) composite n=%d: type-1 word %d is not an Euler/Fermat witness" % (n, b))
        tot["cand"] += len(cands); tot["prime"] += np_
        tot["w0"] += cnt[0]; tot["w1"] += cnt[1]; tot["w2"] += cnt[2]
        print("A (b): rank %d, %d <= n < %d: candidates %d, words %d, prime moduli %d, word types 0/1/2 = %d/%d/%d"
              % (k, LOWER, UPPER[k], len(cands), len(words), np_, cnt[0], cnt[1], cnt[2]))
    print("A (b): total candidates %d, prime moduli %d, word types 0/1/2 = %d/%d/%d, largest certificate witness %d, "
          "largest least prime primitive root %d (n = %d)"
          % (tot["cand"], tot["prime"], tot["w0"], tot["w1"], tot["w2"], bmax, gmaxb[0], gmaxb[1]))

    # (c) extension to every odd prime p <= bound
    nprimes = 0
    gsum = 0
    gm = (0, 0)
    bad = 0
    for p in range(3, bound + 1, 2):
        if not s[p]:
            continue
        nprimes += 1
        rs = distinct_primes_spf(p - 1, spf)
        g = least_prime_root(p, rs, plist)
        if g is None or g >= p:
            bad += 1
            fails.append("(c) no prime primitive root below p=%d" % p)
            continue
        gsum += g
        if g > gm[0]:
            gm = (g, p)
    print("A (c): odd primes p <= %d: %d; each has a prime primitive root < p: %s; sum of least prime primitive "
          "roots %d; largest %d (first at p = %d)" % (bound, nprimes, "yes" if bad == 0 else "NO", gsum, gm[0], gm[1]))

    for x in fails[:20]:
        print("A: FAIL " + x)
    if fails:
        print("A: VERDICT REFUSED (%d failures)" % len(fails))
        return 1
    print("A: VERDICT PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
