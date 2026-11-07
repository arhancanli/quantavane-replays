#!/usr/bin/env python3
"""Independent check of the finite inputs of C_3a > 1.21418 [O2026] (carson-olaf, PR #148; constants/3a.md), written
from section 6 of c3a_two_temperature_proof.md (carson-olaf/c3a-two-temperature-bound @ 033952dc, the commit Palomar
registered) alone: its verify_c3a_two_temperature_121418.py, the GMP program, test_universal_quotient_small.py and the
Lean files were neither read nor run.

The finite proof (section 6): M = <4, 5> cap [0, 80], D = M - M = [-80, 80], kappa(d) = min{a + b : a - b = d},
N = 10^13, x = 8149/10000, y = 42571/50000, Z = P_-(x) = 1 + 2 sum_{d=1..80} x^kappa(d),
n_d = n_-d = floor(N x^kappa(d) / Z), n_0 = N - 2 sum n_d, L = sum_{d>0} n_d kappa(d);
W = the type class of D^N with counts (n_d), V = {a in M^N : sum a_i <= L}. Then
  W subset V - V (symmetric minimum witnesses: both budgets equal L), |W| = N! / prod n_d!,
  |V + V| y^(2L) <= P_+(y)^N  and  |V - V| x^(2L) <= Z^N  (tilted counting),
and the note asserts the integer bounds
  log|W| >= a = 32160395000000,  log|V + V| <= b = 25272245000000,  log|V - V| <= r = 32160396000000,
and the final comparison (6.1)  alpha (r + 8) + b + 2 < a  with alpha = 10709/50000, slack 36384716.28656.
Cyclic compression to a prime p (Bertrand), the dilation V = 2A and Kleinwaks's GHR principle turn these into
C_3a > 1 + alpha = 1.21418.

We compute M, kappa, the counts, n_0 and L in exact integers; log|W| from below with Robbins' bounds on log n!, and
log|V +- V| from above, in 90-digit decimal interval arithmetic rounded outward (the code of our c3e replay, whose
B = 80 row has the same M, N, x and y); then a, b, r and (6.1) exactly. The three finite inclusions and counting bounds
are also checked by brute force on small cases. Python 3.9+, standard library only.

  python3 -I o2026_check.py [--controls]
"""
import decimal
import hashlib
import itertools
import json
import os
import sys
import time
from collections import Counter
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PREC = 90
G = (4, 5)
B, N, X, Y = 80, 10 ** 13, (8149, 10000), (42571, 50000)
STATED = dict(n0=1917647564302, twoL=76435323145630, a=32160395000000, b=25272245000000, r=32160396000000,
              alpha=Fraction(10709, 50000), slack=Fraction("36384716.28656"),
              M=75, MM=155, DD=161, kappa_gaps={1: 9, 2: 18, 3: 13, 6: 14, 7: 17, 11: 19})

LO = decimal.Context(prec=PREC, rounding=decimal.ROUND_FLOOR)
HI = decimal.Context(prec=PREC, rounding=decimal.ROUND_CEILING)
NEAR = decimal.Context(prec=PREC, rounding=decimal.ROUND_HALF_EVEN)


class I:
    """A closed interval [lo, hi] of decimals; every operation rounds outward."""

    def __init__(self, lo, hi=None):
        self.lo, self.hi = Decimal(lo), Decimal(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(LO.add(self.lo, o.lo), HI.add(self.hi, o.hi))

    def __sub__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(LO.subtract(self.lo, o.hi), HI.subtract(self.hi, o.lo))

    def __mul__(self, o):
        o = o if isinstance(o, I) else I(o)
        ps = [(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(LO.multiply(a, b) for a, b in ps), max(HI.multiply(a, b) for a, b in ps))


def iln(x):
    """ln of a positive interval: decimal's ln is correctly rounded, so one unit either way encloses the true value."""
    assert x.lo > 0
    return I(NEAR.next_minus(NEAR.ln(x.lo)), NEAR.next_plus(NEAR.ln(x.hi)))


def ratio(num, den, digits=PREC + 10):
    """The positive rational num/den as an interval, by integer division."""
    shift = max(0, digits - (len(str(num)) - len(str(den))))
    s = 10 ** shift
    f = num * s // den
    return I(LO.divide(Decimal(f), Decimal(s)), HI.divide(Decimal(f + 1), Decimal(s)))


def pi_interval():
    """pi = 16 atan(1/5) - 4 atan(1/239) in integers scaled by 10^(PREC+10); the alternating tails bound the error."""
    s = 10 ** (PREC + 10)

    def atan_inv(n):
        total, term, k, sign = 0, s // n, 1, 1
        while term:
            total += sign * (term // k)
            term //= n * n
            k += 2
            sign = -sign
        return total

    p = 16 * atan_inv(5) - 4 * atan_inv(239)
    return I(LO.divide(Decimal(p - 200), Decimal(s)), HI.divide(Decimal(p + 200), Decimal(s)))


PI = pi_interval()


def lnfact(n):
    """ln n!: exact for n <= 20, else Robbins' bounds 1/(12n+1) < ln n! - (n ln n - n + ln(2 pi n)/2) < 1/(12n)."""
    if n <= 20:
        f = 1
        for k in range(2, n + 1):
            f *= k
        return iln(I(f))
    base = I(n) * iln(I(n)) - I(n) + iln(I(2) * PI * I(n)) * I("0.5")
    return I(LO.add(base.lo, LO.divide(1, 12 * n + 1)), HI.add(base.hi, HI.divide(1, 12 * n)))


def semigroup_mask(Bm, gens=G):
    reach = bytearray(Bm + 1)
    reach[0] = 1
    for v in range(1, Bm + 1):
        reach[v] = 1 if any(v >= g and reach[v - g] for g in gens) else 0
    return [v for v in range(Bm + 1) if reach[v]]


def kappa_table(M, Bm):
    kap, wit = [None] * (Bm + 1), [None] * (Bm + 1)
    for i, a in enumerate(M):
        for b in M[:i + 1]:
            d, c = a - b, a + b
            if kap[d] is None or c < kap[d]:
                kap[d], wit[d] = c, (a, b)
    return kap, wit


def counts(Nn, kap, Bm, xn, xd):
    """n_d = floor(N x^kappa(d) / Z(x)) exactly, Z(x) = 1 + 2 sum x^kappa(d)."""
    K = max(kap[1:])
    Znum = xd ** K + 2 * sum(xn ** kap[d] * xd ** (K - kap[d]) for d in range(1, Bm + 1))
    cnt = {d: Nn * xn ** kap[d] * xd ** (K - kap[d]) // Znum for d in range(1, Bm + 1)}
    n0 = Nn - 2 * sum(cnt.values())
    L = sum(cnt[d] * kap[d] for d in cnt)
    return cnt, n0, L, Fraction(Znum, xd ** K)


def logs(Nn, Mset, kap, cnt, n0, L, x, y):
    """(log|W| interval, upper bound for log|V+V|, upper bound for log|V-V|)."""
    lnW = lnfact(Nn) - lnfact(n0)
    for n, m in Counter(cnt.values()).items():
        lnW = lnW - lnfact(n) * I(2 * m)
    S = sorted({a + b for a in Mset for b in Mset})
    smax = max(S)
    yn, yd = y
    Pplus = ratio(sum(yn ** s * yd ** (smax - s) for s in S), yd ** smax)
    lnVpV = I(Nn) * iln(Pplus) - I(2 * L) * iln(ratio(yn, yd))
    xn, xd = x
    K = max(kap[1:])
    Znum = xd ** K + 2 * sum(xn ** kap[d] * xd ** (K - kap[d]) for d in range(1, len(kap)))
    lnVmV = I(Nn) * iln(ratio(Znum, xd ** K)) - I(2 * L) * iln(ratio(xn, xd))
    return lnW, lnVpV, lnVmV


def brute(Bm=10, Nmax=4):
    """The three finite statements by enumeration on small cases, for M = <4,5> cap [0, Bm] = {0,4,5,8,9,10}: for every
    N <= Nmax and every symmetric count vector (n_-d = n_d in {0, 1} on at most two positive d, n_0 = N - 2 sum n_d),
    with L = sum_{d>0} n_d kappa(d) and V = {v in M^N : sum v <= L}: W subset V - V, and at x in {1/2, 4/5}
    |V + V| x^(2L) <= P_+(x)^N and |V - V| x^(2L) <= Z(x)^N."""
    M = semigroup_mask(Bm)
    kap, wit = kappa_table(M, Bm)
    Dpos = [d for d in range(1, Bm + 1) if kap[d] is not None]
    S = sorted({a + b for a in M for b in M})
    tested = fails = 0
    for Nn in range(2, Nmax + 1):
        for k in (1, 2):
            if 2 * k > Nn:
                continue
            for ds in itertools.combinations(Dpos, k):
                L = sum(kap[d] for d in ds)
                V = [v for v in itertools.product(M, repeat=Nn) if sum(v) <= L]
                VmV = {tuple(a - b for a, b in zip(u, v)) for u in V for v in V}
                VpV = {tuple(a + b for a, b in zip(u, v)) for u in V for v in V}
                word = [0] * (Nn - 2 * k) + list(ds) + [-d for d in ds]
                W = set(itertools.permutations(word))
                good = W <= VmV
                for xr in (Fraction(1, 2), Fraction(4, 5)):
                    Z = 1 + 2 * sum(xr ** kap[d] for d in Dpos)
                    Pp = sum(xr ** t for t in S)
                    good &= len(VpV) * xr ** (2 * L) <= Pp ** Nn and len(VmV) * xr ** (2 * L) <= Z ** Nn
                tested += 1
                fails += not good
    return tested, fails


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def run(root, controls):
    t0 = time.time()
    ok = True
    M = semigroup_mask(B)
    MM = sorted({a + b for a in M for b in M})
    DD = sorted({a - b for a in M for b in M})
    kap, wit = kappa_table(M, B)
    print(f"M = <4,5> cap [0,{B}]: |M| = {len(M)}, |M + M| = {len(MM)}, |M - M| = {len(DD)} (M - M = [-80, 80]: "
          f"{DD == list(range(-B, B + 1))})")
    gaps = {d: kap[d] for d in range(1, B + 1) if d not in M}
    c = (len(M), len(MM), len(DD)) == (STATED["M"], STATED["MM"], STATED["DD"]) and gaps == STATED["kappa_gaps"] \
        and all(kap[d] == d for d in range(1, B + 1) if d in M) and DD == list(range(-B, B + 1))
    print(f"  kappa(d) = d on M, and at the six gaps {gaps}: as the note states: {c}")
    ok &= c
    cnt, n0, L, Z = counts(N, kap, B, *X)
    c = n0 == STATED["n0"] and 2 * L == STATED["twoL"] and all(v > 0 for v in cnt.values()) and n0 > 0
    print(f"counts (exact integers): n_0 = {n0}, 2L = {2 * L}, smallest n_d = {min(cnt.values())} "
          f"(d = {min(cnt, key=cnt.get)}), all positive; n_0 and 2L as stated: {c}")
    ok &= c
    lnW, lnVpV, lnVmV = logs(N, M, kap, cnt, n0, L, X, Y)
    a, b, r = STATED["a"], STATED["b"], STATED["r"]
    ca, cb, cr = Decimal(a) <= lnW.lo, lnVpV.hi <= Decimal(b), lnVmV.hi <= Decimal(r)
    print(f"  log|W| >= {lnW.lo:.6f} >= a = {a}: {ca}")
    print(f"  log|V + V| <= {lnVpV.hi:.6f} <= b = {b}: {cb}")
    print(f"  log|V - V| <= {lnVmV.hi:.6f} <= r = {r}: {cr}")
    ok &= ca and cb and cr
    alpha = STATED["alpha"]
    slack = a - b - 2 - alpha * (r + 8)
    c = slack > 0 and slack == STATED["slack"]
    print(f"  (6.1) alpha (r + 8) + b + 2 < a with alpha = 10709/50000: slack {slack} = {float(slack):.5f} "
          f"(the note's 36384716.28656: {slack == STATED['slack']})")
    ok &= c
    # the same comparison from our own enclosures, with log(4 max A + 1) <= log 1296 + log|V - V|
    num = LO.subtract(LO.subtract(lnW.lo, lnVpV.hi), HI.multiply(Decimal(2), iln(I(2)).hi))
    den = HI.add(lnVmV.hi, iln(I(1296)).hi)
    ours = LO.divide(num, den)
    print(f"  our own exponent bound 1 + (log|W| - log 4 - log|V+V|) / (log 1296 + log|V-V|) >= 1 + {ours:.12f}")
    ok &= num > 0 and Fraction(ours) > alpha
    t, f = brute()
    print(f"brute force on small cases (M = <4,5> cap [0,10], N <= 4, symmetric counts on one or two positive d, "
          f"x in {{1/2, 4/5}}): W subset V - V and both tilted counting bounds hold in {t - f} of {t} cases")
    ok &= f == 0
    if controls:
        print("Red controls:")
        up = Fraction(10710, 50000)
        c1 = a - b - 2 - up * (r + 8) <= 0
        print(f"  C1 alpha = 10710/50000 (1.2142): (6.1) {'fails, as expected' if c1 else 'HOLDS (control failed)'}")
        cnt2, n02, L2, _ = counts(N, kap, B, 8148, 10000)
        c2 = (n02, 2 * L2) != (STATED["n0"], STATED["twoL"])
        print(f"  C2 x = 8148/10000 (one unit lower): n_0 = {n02}, 2L = {2 * L2}: "
              f"{'differ from the stated values, as expected' if c2 else 'UNCHANGED (control failed)'}")
        c3 = Decimal(a + 10 ** 9) > lnW.lo
        print(f"  C3 a raised by 10^9 (to {a + 10 ** 9}): log|W| >= a {'refused, as expected' if c3 else 'ACCEPTED'}")
        Ms = semigroup_mask(10)
        ks, _ = kappa_table(Ms, 10)
        Ls = ks[1]                                         # N = 3, counts n_1 = n_-1 = 1, n_0 = 1: L = kappa(1)
        Vs = [v for v in itertools.product(Ms, repeat=3) if sum(v) <= Ls - 1]
        VmVs = {tuple(p_ - q_ for p_, q_ in zip(u, v)) for u in Vs for v in Vs}
        c4 = not set(itertools.permutations([0, 1, -1])) <= VmVs
        print(f"  C4 a small case (M = <4,5> cap [0,10], N = 3, n_1 = n_-1 = n_0 = 1) with the budget lowered from L = "
              f"{Ls} to {Ls - 1}: W subset V - V {'fails, as expected' if c4 else 'STILL HOLDS (control failed)'}")
        ok &= c1 and c2 and c3 and c4
    print(f"RESULT o2026: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.1f} s)")
    return ok


if __name__ == "__main__":
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    root = args[0] if args else os.path.join(HERE, "certs")
    sys.exit(0 if run(root, "--controls" in sys.argv) else 1)
