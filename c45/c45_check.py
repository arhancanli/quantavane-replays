#!/usr/bin/env python3
"""Independent replay of upper bounds for Romanoff's constant C45 (the 45a page), from the page's obstruction inequality
and the prime sets alone. No program of any certificate was read or run; the base histograms of vidigoat's packet are
read as data, only to compare them with ours at the end.

The inequality (45a page). P a finite set of odd primes, M = prod P, T = lcm of ord_q(2). For a mod M,
F(a) = {r mod T : gcd(a - 2^r, M) = 1}, nu(a) = |F(a)|, delta(nu) = #{a mod M : nu(a) = nu}, and
    dbar(A) <= S := sum_nu delta(nu) min( 1/(2M), nu/(T phi(M) log 2) ).
(Re-derived: an odd n = a mod M with n = p + 2^k, p > max P, has k mod T in F(a); there are about nu(a) log2(x)/T such
k <= log2 x, each with at most (1 + o(1)) x/(phi(M) log x) primes p <= x in the class a - 2^k, and there are about
x/(2M) odd n <= x in the class.) Replacing log 2 by a rational lambda < log 2 can only increase S.

How S is computed. Per prime q of order m, a mod q is a power 2^j (exactly one j mod m) for m of the q residues,
which removes the class r = j (mod m) from F(a), and removes nothing for the other q - m residues.
  * Base primes have orders dividing 5040. Work on Z/5040 (a multiple of their lcm T0; nu and T scale together). The
    primes are grouped by the modulus their order divides, G1 | 12, G2 | 36, G3 | 60, G4 | 48, GA | 28. For t in Z/12
    the residues r = t (mod 12) correspond one to one to their components mod 36, 60, 48, 28 above t (3, 5, 4, 7 of
    them; 3*5*4*7 = 420 = 5040/12), so nu0 = sum_t prod_g u_g(t), u_g(t) = the uncovered residues of group g above t.
    Each group's distribution of (u_g(t))_t is exact (here, by a DP over its primes on the covered set); c45_comb.c
    combines them (affine symmetry of Z/12, exact 128-bit weights) into the histogram delta0 of nu0.
  * The histogram is checked: total M0, and the first two moments against closed forms (sum nu0 = 5040 prod (q - 1);
    sum nu0^2 = 5040 sum_d prod_q (q - 1 if ord_q | d else q - 2)); and on small prime sets against brute force over
    every a mod M0 from the definition (selftest).
  * Cluster primes have prime order c coprime to 5040 and to the other cluster orders; by CRT, nu = nu0 prod_c h_c with
    h_c = c minus the number of distinct classes mod c that the cluster's primes remove (exact distribution by a DP).
  * S = E[min(1/2, kappa nu0 H)], H = prod h_c, kappa = (M/phi(M)) / (T lambda), T = 5040 prod c, lambda =
    sum_{k <= 120} 1/(k 2^k) < log 2: exact rational arithmetic, the clusters split in two halves, a two-pointer sweep.

Usage: python3 -I c45_check.py ROW BINDIR [CERTS] [--controls]     rows: griego, yoo, pin21
       python3 -I c45_check.py selftest BINDIR
"""
import hashlib
import os
import subprocess
import sys
import time
from fractions import Fraction
from math import gcd, prod

ROWS = {
    # 45a page, [G2026] (superseded row): the 36-prime set printed on the page
    "griego": dict(
        primes=[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 61, 73, 89, 167, 223, 233, 263, 359, 383, 431, 439, 479,
                1103, 1433, 1913, 2089, 2351, 4513, 5737, 8191, 9719, 176383, 178481],
        printed="0.490249407811155", where="45a page, [G2026] (36 primes, as printed there)",
        hist=("base_griego13.txt", "f35f2f3d0a14c04ecc0a8f65a0b2e8b469f358a344b3e0038bee83e1e65c7d08")),
    # 45a page, [Y2026] (the current record row): the 38-prime set of Yoo's README (andrew-yoo/romanoff, tag v3)
    "yoo": dict(
        primes=[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 61, 73, 89, 167, 223, 233, 257, 263, 359, 383, 431, 439,
                479, 503, 1103, 1399, 1433, 2089, 2351, 2687, 3391, 4513, 6361, 8191, 9719],
        printed="0.490180063290061", where="45a page, [Y2026] (38 primes, Yoo's README at tag v3)",
        hist=("base_yoo14.txt", "4ac1d42106f662d810a47876ff60fbc93c9e7b2d0f96c4b7c4fa8c3351d0aeee")),
    # PR #281, row [P2026]: Yoo's 14 base primes + 43, 127, 241, 113, 109, 97, 151, and the union of the clusters
    "pin21": dict(
        primes=[3, 5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 61, 73, 257, 43, 127, 241, 113, 109, 97, 151,
                23, 89, 8191, 47, 178481, 233, 1103, 2089, 223, 431, 9719, 2351, 4513, 6361, 439, 2687, 167, 3391, 263,
                359, 1433, 383, 1399, 479, 1913, 5737, 176383, 503],
        printed="0.489596115074071", where="PR #281, row [P2026] (21 base primes and 28 cluster primes; the row says 39 primes)",
        hist=("base_final.txt", "238174d0239c9cec3718238416d62d201405471c2a63127ba2397ebba1fcc41d")),
}

SLOTS = (("G1", 12), ("G2", 36), ("G3", 60), ("G4", 48), ("GA", 28))
FRAME = 5040


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def order2(q):
    k, x = 1, 2 % q
    while x != 1:
        x = x * 2 % q
        k += 1
    return k


def split(primes):
    """base primes (order | 5040) and clusters {c: [q]} (c prime, coprime to 5040); anything else is refused."""
    if len(set(primes)) != len(primes):
        raise ValueError("a prime is listed twice")
    base, clusters = [], {}
    for q in primes:
        if not (is_prime(q) and q > 2):
            raise ValueError(f"{q} is not an odd prime")
        m = order2(q)
        if FRAME % m == 0:
            base.append(q)
        elif is_prime(m) and gcd(m, FRAME) == 1:
            clusters.setdefault(m, []).append(q)
        else:
            raise ValueError(f"ord_{q}(2) = {m} neither divides 5040 nor is a prime coprime to it")
    return base, clusters


def assign(base):
    groups = {name: [] for name, _ in SLOTS}
    for q in base:
        m = order2(q)
        for name, L in SLOTS:
            if L % m == 0:
                groups[name].append(q)
                break
        else:
            raise ValueError(f"ord_{q}(2) = {m} divides none of 12, 36, 60, 48, 28")
    return groups


def group_vectors(primes, L):
    """exact distribution of the vector u(t), t mod 12 (t mod 4 when L = 28), of uncovered residues of Z/L above t."""
    states = {0: 1}
    for q in primes:
        m = order2(q)
        assert L % m == 0
        aps = [sum(1 << r for r in range(j, L, m)) for j in range(m)]
        new = {}
        for s, w in states.items():
            new[s] = new.get(s, 0) + w * (q - m)
            for a in aps:
                new[s | a] = new.get(s | a, 0) + w
        states = new
    g = gcd(L, 12)
    n = 12 if g == 12 else 4
    assert g in (12, 4)
    out = {}
    for s, w in states.items():
        v = tuple(sum(1 for r in range(t, L, g) if not (s >> r) & 1) for t in range(n))
        out[v] = out.get(v, 0) + w
    assert sum(out.values()) == prod(primes)
    return out


def base_histogram(base, bindir, log):
    groups = assign(base)
    M0 = prod(base)
    if M0 >= 2 ** 127:
        raise ValueError("M0 does not fit the 128-bit weights")
    t0 = time.time()
    lines = []
    for name, L in SLOTS:
        d = group_vectors(groups[name], L)
        lines.append(f"{name} {len(d)}")
        lines += [" ".join(map(str, v)) + f" {w}" for v, w in sorted(d.items())]
        log(f"  {name} (orders dividing {L}): {groups[name]} -> {len(d)} distinct vectors")
    t1 = time.time()
    run = subprocess.run([os.path.join(bindir, "c45_comb")], input="\n".join(lines) + "\n", capture_output=True,
                         text=True, check=True)
    hist = {}
    for line in run.stdout.split("\n"):
        if line.strip():
            nu, c = line.split()
            hist[int(nu)] = int(c)
    for line in run.stderr.strip().split("\n"):
        if not line.startswith("  "):
            log(f"  c45_comb: {line}")
    log(f"  histogram of nu0 on Z/5040: {len(hist)} values ({t1 - t0:.0f} s for the groups, {time.time() - t1:.0f} s "
        f"for the combination)")
    return hist, M0


def moment_checks(hist, base, log):
    M0 = prod(base)
    ms = [order2(q) for q in base]
    s0 = sum(hist.values())
    s1 = sum(nu * c for nu, c in hist.items())
    s2 = sum(nu * nu * c for nu, c in hist.items())
    e1 = FRAME * prod(q - 1 for q in base)
    e2 = FRAME * sum(prod((q - 1) if d % m == 0 else (q - 2) for q, m in zip(base, ms)) for d in range(FRAME))
    ok = (s0, s1, s2) == (M0, e1, e2)
    log(f"  checks: total = M0, sum nu0 = 5040 prod (q - 1), sum nu0^2 = 5040 sum_d prod_q (q - 1 or q - 2): "
        f"{'all three exact' if ok else 'MISMATCH'}")
    return ok


def brute_histogram(base):
    """from the definition: every a mod M0, F(a) = {r mod T0 : gcd(a - 2^r, M0) = 1}; nu scaled to Z/5040."""
    M0 = prod(base)
    T0 = 1
    for q in base:
        m = order2(q)
        T0 = T0 * m // gcd(T0, m)
    pw = [pow(2, r, M0) for r in range(T0)]
    hist = {}
    for a in range(M0):
        nu = sum(1 for r in range(T0) if gcd((a - pw[r]) % M0, M0) == 1) * (FRAME // T0)
        hist[nu] = hist.get(nu, 0) + 1
    return hist


def cluster_dist(c, qs):
    d = {0: 1}
    for q in qs:
        if order2(q) != c:
            raise ValueError(f"ord_{q}(2) = {order2(q)}, not the cluster order {c}")
        new = {}
        for D, w in d.items():
            new[D] = new.get(D, 0) + w * (q - c) + w * D
            new[D + 1] = new.get(D + 1, 0) + w * (c - D)
        d = new
    out = {c - D: w for D, w in d.items()}
    assert sum(out.values()) == prod(qs)
    return out


def product_dist(dists):
    out = {1: 1}
    for d in dists:
        new = {}
        for h, w in out.items():
            for k, v in d.items():
                new[h * k] = new.get(h * k, 0) + w * v
        out = new
    return out


def lam():
    """lambda = sum_{k <= 120} 1/(k 2^k), a partial sum of log 2 = sum_k 1/(k 2^k), so lambda < log 2."""
    return sum(Fraction(1, k * 2 ** k) for k in range(1, 121))


def evaluate(hist, base, clusters, lam_value):
    """exact S = E[min(1/2, kappa nu0 H)]; returns S (a Fraction)."""
    M0 = prod(base)
    allq = base + [q for qs in clusters.values() for q in qs]
    T = FRAME * prod(clusters)
    kappa = Fraction(prod(allq), prod(q - 1 for q in allq)) / (T * lam_value)
    P, Q = kappa.numerator, kappa.denominator
    dists = sorted((cluster_dist(c, qs) for c, qs in clusters.items()), key=len)
    half1, half2 = [], []
    n1 = n2 = 1
    for d in reversed(dists):     # balance the numbers of combinations
        if n1 <= n2:
            half1.append(d)
            n1 *= len(d)
        else:
            half2.append(d)
            n2 *= len(d)
    H1 = product_dist(half1)
    H2 = sorted(product_dist(half2).items())
    if len(H1) > len(H2):
        H1, H2 = dict(H2), sorted(H1.items())
    W2, X2 = [0], [0]
    for h, w in H2:
        W2.append(W2[-1] + w)
        X2.append(X2[-1] + w * h)
    tot1, tot2 = sum(H1.values()), W2[-1]
    nus = sorted(hist.items())
    A = B = 0
    for h1, w1 in H1.items():
        i = len(H2)
        for nu, w0 in nus:
            while i > 0 and 2 * P * nu * h1 * H2[i - 1][0] > Q:   # capped: kappa nu h1 h2 > 1/2
                i -= 1
            A += w0 * w1 * nu * h1 * X2[i]
            B += w0 * w1 * (tot2 - W2[i])
    return (Fraction(P * A, Q) + Fraction(B, 2)) / (M0 * tot1 * tot2), len(H1) * len(H2)


def dec(x, n):
    i = x.numerator * 10 ** n // x.denominator
    s = str(i).rjust(n + 1, "0")
    return s[:-n] + "." + s[-n:]


def read_hist(path, sha):
    raw = open(path, "rb").read()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha:
        raise SystemExit(f"SHA-256 of {path} is {got}, expected {sha}")
    h = {}
    for line in raw.decode().split("\n"):
        if line.strip():
            nu, c = line.split()
            h[int(nu)] = int(c)
    return h


def selftest(bindir):
    ok = True
    print("Self-test: the histogram pipeline against brute force over every a mod M0 (from the definition)")
    for base in ([3, 5, 7, 31, 73], [3, 11, 17, 127], [5, 13, 19, 29], [3, 7, 241, 11], [37, 97, 3]):
        t0 = time.time()
        h1, _ = base_histogram(base, bindir, log=lambda s: None)
        h2 = brute_histogram(base)
        good = h1 == h2
        print(f"  base {base} (M0 = {prod(base)}): pipeline {len(h1)} values, brute force {len(h2)} values, "
              f"{'identical' if good else 'DIFFERENT'} ({time.time() - t0:.0f} s)")
        ok &= good
    print("  cluster distributions against direct enumeration of (a mod q) for q in the cluster:")
    for c, qs in ((11, [23, 89]), (29, [233, 1103])):
        d = cluster_dist(c, qs)
        direct = {}
        pw = [{pow(2, j, q): j for j in range(c)} for q in qs]
        for a1 in range(qs[0]):
            for a2 in range(qs[1]):
                fired = {p[a] for p, a in zip(pw, (a1, a2)) if a in p}
                direct[c - len(fired)] = direct.get(c - len(fired), 0) + 1
        good = d == direct
        print(f"    order {c}, primes {qs}: {'identical' if good else 'DIFFERENT'}")
        ok &= good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}")
    return ok


def check_row(key, bindir, certs, controls):
    row = ROWS[key]
    t0 = time.time()
    ok = True
    log = print
    base, clusters = split(row["primes"])
    log(f"{row['where']}: C45 < {row['printed']}")
    log(f"  {len(row['primes'])} primes: {len(base)} base primes (orders dividing 5040), {sum(map(len, clusters.values()))} "
        f"in {len(clusters)} clusters of prime orders {sorted(clusters)}")
    hist, M0 = base_histogram(base, bindir, log)
    ok &= moment_checks(hist, base, log)
    path = os.path.join(certs, key, row["hist"][0]) if certs else None
    if path and os.path.exists(path):
        theirs = read_hist(path, row["hist"][1])
        T0 = 1
        for q in base:
            T0 = T0 * order2(q) // gcd(T0, order2(q))
        scaled = {nu * (FRAME // T0): c for nu, c in theirs.items()}
        same = scaled == hist
        log(f"  the packet's {row['hist'][0]} (sha256 {row['hist'][1][:8]}..., {len(theirs)} entries, on Z/{T0}), "
            f"read as data: {'identical to ours' if same else 'DIFFERENT FROM OURS'}")
        ok &= same
    else:
        log("  (the packet's histogram file is not present; comparison skipped)")
    lv = lam()
    t1 = time.time()
    S, ncomb = evaluate(hist, base, clusters, lv)
    log(f"  lambda = sum_(k <= 120) 1/(k 2^k) < log 2; {ncomb:,} (h1, h2) pairs of the cluster halves; "
        f"{time.time() - t1:.0f} s")
    log(f"  S = {dec(S, 30)}... (exact rational)")
    b = Fraction(row["printed"])
    nd = len(row["printed"].split(".")[1])
    good = S <= b and S > b - Fraction(1, 10 ** nd)
    log(f"  S <= {row['printed']}, and the printed value is S rounded up at its last digit: {good}")
    ok &= good
    if controls:
        log("Red controls:")
        drop = base[-1]
        hb, _ = base_histogram([q for q in base if q != drop], bindir, log=lambda s: None)
        Sd, _ = evaluate(hb, [q for q in base if q != drop], clusters, lv)
        good = Sd > b
        log(f"  without the base prime {drop}: S = {dec(Sd, 15)}..., above {row['printed']}: {good}")
        ok &= good
        c0 = max(clusters, key=lambda c: len(clusters[c]))
        cl = {c: (qs[:-1] if c == c0 else qs) for c, qs in clusters.items()}
        Sc, _ = evaluate(hist, base, cl, lv)
        good = Sc > b
        log(f"  without the cluster prime {clusters[c0][-1]} (order {c0}): S = {dec(Sc, 15)}..., above: {good}")
        ok &= good
        try:
            split(row["primes"] + [53])
            log("  the prime 53 (order 52) added to the set: NOT REFUSED")
            ok = False
        except ValueError:
            log("  the prime 53 added to the set: refused (order 52 neither divides 5040 nor is prime), as expected")
        try:
            wrong = {c: list(qs) for c, qs in clusters.items()}
            wrong[min(wrong)] = wrong[min(wrong)] + [53]
            evaluate(hist, base, wrong, lv)
            log("  53 placed in the cluster of order %d: NOT REFUSED" % min(wrong))
            ok = False
        except ValueError:
            log(f"  53 placed in the cluster of order {min(clusters)}: refused (its order is 52), as expected")
        bad = dict(hist)
        k = sorted(bad)[len(bad) // 2]
        bad[k] -= 1
        bad[k + 1] = bad.get(k + 1, 0) + 1
        good = not moment_checks(bad, base, log=lambda s: None)
        log(f"  one a moved from nu0 = {k} to {k + 1} in our histogram: the moment checks fail, as expected: {good}")
        ok &= good
    log(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) < 2:
        print(__doc__)
        return 2
    if args[0] == "selftest":
        return 0 if selftest(args[1]) else 1
    if args[0] not in ROWS:
        print(f"unknown row {args[0]}; rows: {', '.join(ROWS)}")
        return 2
    certs = args[2] if len(args) > 2 else None
    return 0 if check_row(args[0], args[1], certs, "--controls" in sys.argv) else 1


if __name__ == "__main__":
    sys.exit(main())
