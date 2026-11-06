"""Independent check of the C_3e lower bounds 1.27255796 (B = 80) and 1.27255802 (B = 160) in
mizopl/ghr-c3a-odd-base-finite-bound/c3e_finite_certificate (6 Nov 2026), written from Theorems C1-C3 of that folder's
README alone: its check_c3e.py was neither read nor run, and the two certificate files are compared, never used.

  ./fetch.sh && python3 -I c3e_check.py [CERTS_ROOT] [--controls]

C_3e = sup log|A - A| / log|A + A| over finite sets A of integers with |A| >= 2 (constants/3e.md, Problem 6.43 of
[GGSWT2025]). For the B, N and tilts x, y that the README states: M = <4, 5> cap [0, B], S = M + M, Q = 2B + 1 and
kappa(d) = min{a + b : a, b in M, a - b = d} for d = 1..B; counts n_d = floor(N x^kappa(d) / Z(x)) with
Z(x) = 1 + 2 sum_d x^kappa(d) (the rule of [O2026], as in our c3a replay), n_0 = N - 2 sum n_d, L = sum n_d kappa(d),
all in exact integers. Theorem C2 gives |A_N - A_N| >= D_N = N! / (n_0! prod (n_d!)^2) and Theorem C3 gives
|A_N + A_N| <= S_N(y) = P_+(y)^N / y^(2L) with P_+(y) = sum_{s in S} y^s, so C_3e >= log D_N / log S_N(y). log D_N is
bounded below (Robbins' bounds on log n!) and log S_N(y) above, in decimal interval arithmetic at 90 digits rounded
outward, from the exact rationals x, y and P_+(y). Theorems C1-C3 are also checked by brute force on small cases
(B = 12, N <= 3). Python 3.9+, standard library only."""
import decimal
import itertools
import json
import os
import sys
import time
from collections import Counter
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
PREC = 90
G = (4, 5)

# the rows as the README states them; n_0, L and the logs are ours, compared with the certificate files at the end
ROWS = {
    "80_161": dict(B=80, N=10 ** 13, x=(8149, 10000), y=(42571, 50000), claim="1.27255796",
                   stated="1.27255796845801701855", json="c3e_finite_certificate_80_161.json",
                   readme=dict(n0=1917647564302, L=38217661572815)),
    "160_321": dict(B=160, N=10 ** 15, x=(814897, 10 ** 6), y=(851418, 10 ** 6), claim="1.27255802",
                    stated="1.27255802020996480086", json="c3e_finite_certificate_160_321.json", readme=None),
}
HRY = "1.27155330316361201881"   # log(1 + sqrt 2) / log 2, [HRY1999], the page's best lower bound

# ---- interval arithmetic on decimals (as in our c3a replay) -------------------------------------------------------
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

    def __truediv__(self, o):
        o = o if isinstance(o, I) else I(o)
        assert o.lo > 0 or o.hi < 0, "division by an interval containing 0"
        ps = [(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(LO.divide(a, b) for a, b in ps), max(HI.divide(a, b) for a, b in ps))


def iln(x):
    """ln of a positive interval: decimal's ln is correctly rounded, so one unit either way encloses the true value."""
    assert x.lo > 0
    a, b = NEAR.ln(x.lo), NEAR.ln(x.hi)
    return I(NEAR.next_minus(a), NEAR.next_plus(b))


def ratio(num, den, digits=PREC + 10):
    """The positive rational num/den as an interval of width 10^-digits relative to its size, by integer division."""
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
    err = 200
    return I(LO.divide(Decimal(p - err), Decimal(s)), HI.divide(Decimal(p + err), Decimal(s)))


PI = pi_interval()


def lnfact(n):
    """ln n! for an integer n >= 0: exact for n <= 20, Robbins' bounds 1/(12n+1) < ln n! - (n ln n - n + ln(2 pi n)/2)
    < 1/(12n) above."""
    if n <= 20:
        f = 1
        for k in range(2, n + 1):
            f *= k
        return iln(I(f))
    base = I(n) * iln(I(n)) - I(n) + iln(I(2) * PI * I(n)) * I("0.5")
    return I(LO.add(base.lo, LO.divide(1, 12 * n + 1)), HI.add(base.hi, HI.divide(1, 12 * n)))


# ---- the construction ---------------------------------------------------------------------------------------------
def semigroup_mask(B):
    """M = <4, 5> cap [0, B], by reachability."""
    reach = bytearray(B + 1)
    reach[0] = 1
    for v in range(1, B + 1):
        reach[v] = 1 if any(v >= g and reach[v - g] for g in G) else 0
    return [v for v in range(B + 1) if reach[v]]


def kappa_table(M, B):
    """kappa(d) = min{a + b : a, b in M, a - b = d} for 0 <= d <= B (None if d is not in M - M), with a witness."""
    kap, wit = [None] * (B + 1), [None] * (B + 1)
    for i, a in enumerate(M):
        for b in M[:i + 1]:
            d, c = a - b, a + b
            if kap[d] is None or c < kap[d]:
                kap[d], wit[d] = c, (a, b)
    return kap, wit


def counts(N, kap, B, xn, xd):
    """n_d = floor(N x^kappa(d) / Z(x)), Z(x) = 1 + 2 sum_{d=1..B} x^kappa(d), exactly: with K = max kappa,
    n_d = floor(N xn^kappa xd^(K - kappa) / Znum), Znum = xd^K + 2 sum xn^kappa(d) xd^(K - kappa(d))."""
    K = max(kap[1:])
    Znum = xd ** K + 2 * sum(xn ** kap[d] * xd ** (K - kap[d]) for d in range(1, B + 1))
    cnt = {d: N * xn ** kap[d] * xd ** (K - kap[d]) // Znum for d in range(1, B + 1)}
    n0 = N - 2 * sum(cnt.values())
    L = sum(cnt[d] * kap[d] for d in cnt)
    return cnt, n0, L


def bound(N, S, cnt, n0, L, yn, yd):
    """(log D_N interval, log S_N(y) interval, lower bound of their ratio)."""
    lnD = lnfact(N) - lnfact(n0)
    for n, m in Counter(cnt.values()).items():
        lnD = lnD - lnfact(n) * I(2 * m)
    smax = max(S)
    P = ratio(sum(yn ** s * yd ** (smax - s) for s in S), yd ** smax)        # P_+(y), an exact rational
    lny = iln(ratio(yn, yd))
    lnS = I(N) * iln(P) - I(2 * L) * lny
    assert lnD.lo > 0 and lnS.lo > 0
    return lnD, lnS, LO.divide(lnD.lo, lnS.hi)


# ---- Theorems C1-C3 by brute force on small cases -----------------------------------------------------------------
def evalq(v, Q):
    return sum(c * Q ** i for i, c in enumerate(v))


def witness_ok(Mset, d, a, b, kap):
    """(a, b) in M^2 with a - b = d and a + b = kappa(d)."""
    return a in Mset and b in Mset and a - b == d and a + b == kap[d]


def small_cases(Q_of=lambda B: 2 * B + 1, B=15, slack=0):
    """(all hold, lines): for N = 1..3 and several L, |A + A| = |V + V| and |A - A| = |V - V| (C1); every word of W_N
    lies in V - V, so |V - V| >= D_N (C2, with nonnegative counts n_0 > 0, and L = kappa(d) + slack); |V + V| <=
    #{s in S^N : sum s <= 2L} <= S_N(y) for three y (C3)."""
    M = semigroup_mask(B)
    S = sorted({a + b for a in M for b in M})
    kap, wit = kappa_table(M, B)
    Q = Q_of(B)
    ok, lines = True, []
    for N in (1, 2, 3):
        for L in sorted({0, 5, 13, 2 * B, N * B}):
            V = [v for v in itertools.product(M, repeat=N) if sum(v) <= L]
            A = [evalq(v, Q) for v in V]
            VpV = {tuple(a + b for a, b in zip(u, v)) for u in V for v in V}
            VmV = {tuple(a - b for a, b in zip(u, v)) for u in V for v in V}
            ApA = {a + b for a in A for b in A}
            AmA = {a - b for a in A for b in A}
            c1 = len(ApA) == len(VpV) and len(AmA) == len(VmV)
            box = sum(1 for s in itertools.product(S, repeat=N) if sum(s) <= 2 * L)
            c3 = len(VpV) <= box and all(
                len(VpV) <= sum(y ** s for s in S) ** N / y ** (2 * L) * (1 + 1e-12) for y in (0.5, 0.8, 0.9))
            ok &= c1 and c3
            if not c1:
                lines.append(f"N = {N}, L = {L}: |A+A| {len(ApA)} vs |V+V| {len(VpV)}, |A-A| {len(AmA)} vs |V-V| {len(VmV)}")
    # C2: one +-d pair and n_0 = N - 2 (N = 3), for every d; L = kappa(d) exactly
    for d in range(1, B + 1):
        N, L = 3, kap[d] + slack
        V = {v for v in itertools.product(M, repeat=N) if sum(v) <= L}
        VmV = {tuple(a - b for a, b in zip(u, v)) for u in V for v in V}
        words = set(itertools.permutations((0, d, -d)))
        al, be = wit[d]
        inside = all(tuple(al if c == d else be if c == -d else 0 for c in w) in V and
                     tuple(be if c == d else al if c == -d else 0 for c in w) in V for w in words)
        ok &= inside and words <= VmV and len(VmV) >= len(words) == 6
    return ok, lines


# ---- one row --------------------------------------------------------------------------------------------------------
def check_row(name, root, controls):
    t0 = time.time()
    row = ROWS[name]
    B, N = row["B"], row["N"]
    (xn, xd), (yn, yd) = row["x"], row["y"]
    Q = 2 * B + 1
    print(f"Row {name}: B = {B}, Q = 2B + 1 = {Q}, N = {N:.0e}, x = {xn}/{xd}, y = {yn}/{yd}")
    ok = True

    def check(what, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  check {what}: {'PASS' if cond else 'FAIL'}{' (' + detail + ')' if detail else ''}")

    cert = None
    p = os.path.join(root, "mizopl", row["json"]) if root else None
    if p and os.path.exists(p):
        cert = json.load(open(p))      # compared only; nothing below is taken from it
    M = semigroup_mask(B)
    S = sorted({a + b for a in M for b in M})
    kap, wit = kappa_table(M, B)
    check("hypotheses: B >= 12, 0 and B in M, M - M = [-B, B] (every d = 1..B has a witness)",
          B >= 12 and M[0] == 0 and M[-1] == B and all(kap[d] is not None for d in range(1, B + 1)),
          f"|M| = {len(M)}, gaps below 12: {[v for v in range(12) if v not in M]}, |M + M| = {len(S)}, max {S[-1]}")
    Mset = set(M)
    check("every witness (alpha_d, beta_d) lies in M^2 with difference d and sum kappa(d)",
          all(witness_ok(Mset, d, *wit[d], kap) for d in range(1, B + 1)), f"kappa from {min(kap[1:])} to {max(kap[1:])}")
    if cert:
        check("|M|, |M + M|, max(M + M), the number of positive differences and |D| as the certificate states",
              (len(M), len(S), S[-1], B, 2 * B + 1) == (cert["mask_size"], cert["sum_support_size"],
                                                        cert["sum_support_max"], cert["positive_differences_count"],
                                                        cert["difference_alphabet_size"]))
    cnt, n0, L = counts(N, kap, B, xn, xd)
    check("every count n_d is positive and n_0 = N - 2 sum n_d > 0 (Theorem C2's hypotheses)",
          min(cnt.values()) > 0 and n0 > 0, f"n_0 = {n0}, smallest n_d = {min(cnt.values())} (d = "
          f"{min(cnt, key=cnt.get)}), L = {L}")
    if row["readme"]:
        check("n_0 and 2L as the README states", (n0, 2 * L) == (row["readme"]["n0"], 2 * row["readme"]["L"]))
    if cert:
        check("n_0, L and 2L as the certificate states", (n0, L, 2 * L) == (cert["n_0"], cert["L"], cert["two_L"]))
    lnD, lnS, lower = bound(N, S, cnt, n0, L, yn, yd)
    print(f"  log D_N >= {LO.plus(lnD.lo):.30f}" + (f"  (certificate: {cert['log_D_N_lower_bound']})" if cert else ""))
    print(f"  log S_N(y) <= {HI.plus(lnS.hi):.30f}" + (f"  (certificate: {cert['log_S_N_upper_bound']})" if cert else ""))
    print(f"  C_3e >= log D_N / log S_N(y) >= {LO.plus(lower):.30f}")
    print(f"           README, pure-integer bounds: {row['stated']}")
    check(f"C_3e > {row['claim']}, the printed bound", lower > Decimal(row["claim"]))
    check(f"and above [HRY1999]'s {HRY[:12]}", lower > Decimal(HRY), f"margin {LO.subtract(lower, Decimal(HRY)):.8f}")
    print(f"  (our enclosures: log D_N in [lo, lo + {HI.subtract(lnD.hi, lnD.lo):.1e}], log S_N(y) in "
          f"[hi - {HI.subtract(lnS.hi, lnS.lo):.1e}, hi])")
    if cert:
        # the certificate's two logs, against our enclosures, allowing half a unit in the last digit it prints
        for k, enc, low in (("log_D_N_lower_bound", lnD, True), ("log_S_N_upper_bound", lnS, False)):
            theirs = Decimal(cert[k])
            half = LO.divide(Decimal(1).scaleb(theirs.as_tuple().exponent), 2)
            near, far = (enc.lo, enc.hi) if low else (enc.hi, enc.lo)      # our end on its side, and the other end
            if (theirs <= near) if low else (theirs >= near):
                cond, detail = True, f"valid as printed, {abs(HI.subtract(near, theirs)):.2e} looser than ours"
            elif (theirs <= far) if low else (theirs >= far):
                cond, detail = True, f"inside our enclosure of width {HI.subtract(enc.hi, enc.lo):.1e}"
            else:
                past = abs(HI.subtract(theirs, far))
                cond = (LO.subtract(theirs, half) <= far) if low else (HI.add(theirs, half) >= far)
                detail = (f"as printed it is {past:.2e} {'above' if low else 'below'} the exact value, less than half a "
                          f"unit ({half:.0e}) of its last digit: rounded to nearest, not {'down' if low else 'up'}")
            check(f"the certificate's {k.replace('_', ' ')} is {'a lower' if low else 'an upper'} bound up to the rounding "
                  f"of its printed digits", cond, detail)
    if controls:
        print("Red controls (the same check functions, on altered data):")
        bad = {k: v + v // 4 for k, v in cnt.items()}
        refused = not (N - 2 * sum(bad.values()) > 0)
        print(f"  R1 every count raised by 25% (n_0 = {N - 2 * sum(bad.values())}): "
              f"{'refused by the count check, as expected' if refused else 'ACCEPTED (control failed)'}")
        ok &= refused
        refused = not (lower > Decimal(row["claim"]) + Decimal("0.000001"))
        print(f"  R2 a claim 1e-6 above the printed one is not implied by our bound: "
              f"{'refused, as expected' if refused else 'ACCEPTED (control failed)'}")
        ok &= refused
        d = next(d for d in range(1, B + 1) if wit[d][0] + 1 in Mset and wit[d][1] + 1 in Mset)
        a, b = wit[d]
        refused = not witness_ok(Mset, d, a + 1, b + 1, kap)
        print(f"  R3 the witness for d = {d} shifted to ({a + 1}, {b + 1}) (same difference, cost kappa + 2): "
              f"{'refused by the witness check, as expected' if refused else 'ACCEPTED (control failed)'}")
        ok &= refused
    print(f"RESULT: {name} {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.1f} s)")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    root = args[0] if args else os.path.join(HERE, "certs")
    controls = "--controls" in sys.argv
    t0 = time.time()
    good, lines = small_cases()
    print(f"Theorems C1-C3 by brute force (B = 15, Q = 31, N = 1..3, up to five L each; C2 for every d): "
          f"{'PASS' if good else 'FAIL'} ({time.time() - t0:.1f} s)")
    short = {B: [d for d in range(1, B + 1) if kappa_table(semigroup_mask(B), B)[0][d] is None] for B in range(12, 16)}
    print("  note: M - M = [-B, B] needs B >= 15, not B >= 12 as the README's section 2 says: the d missing from M - M are "
          + ", ".join(f"B = {B}: {v or 'none'}" for B, v in short.items()) + " (B = 80 and 160 are unaffected)")
    for line in lines:
        print("  " + line)
    ok = good
    if controls:
        bad, lines = small_cases(Q_of=lambda B: 2 * B)
        print(f"  R0 the same with Q = 2B = 30 (digit differences can reach Q, so evaluation can merge two "
              f"differences): {'refused by the brute force, as expected' if not bad else 'ACCEPTED (control failed)'}")
        for line in lines[:2]:
            print("     " + line)
        ok &= not bad
        bad, _ = small_cases(slack=-1)
        print(f"  R0b the same with L = kappa(d) - 1 in Theorem C2's small cases (the words' preimages no longer fit "
              f"under L): {'refused by the brute force, as expected' if not bad else 'ACCEPTED (control failed)'}")
        ok &= not bad
    for name in ROWS:
        ok &= check_row(name, root, controls)
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
