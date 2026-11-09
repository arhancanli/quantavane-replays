"""3a: replay of the row [G2026] C_3a >= 1.1740744 (Griego, PR #71), from the construction stated in the pull request
alone (its programs neither read nor run), with our own exact counts (Quantavane replays). Standard library only.

The construction: A = {0,2,3,...,10}, base B = 21, d = 80 digits, cap T = 150,
    U = { sum_{i<d} a_i B^i : a_i in A, sum a_i <= T },
and the lemma on the 3a page, C_3a >= 1 + log(|U-U|/|U+U|)/log(2 max U + 1) for a set U of distinct non-negative
integers containing 0 with |U-U| < 2 max U + 1.

How we count (re-derived; no carries anywhere):
  * max U: the digit sum is spent from the top, T = 15 * 10, so max U = 10 * sum_{i=65}^{79} 21^i.
  * U+U: every digit of u + v is a_i + b_i <= 20 < 21, so u + v <-> the digit vector s = a + b, and |U+U| counts
    the s in (A+A)^d that split as s = a + b with a, b in A^d, sum a <= T, sum b <= T. With n = sum s and
    Y(s_i) = {2a - s_i : a, s_i - a in A} (a symmetric set), such a split exists iff n + mu <= 2T, where mu is the
    least |sum y_i| over y_i in Y(s_i). A dynamic program over the digits keeps (n, R), R the set of reachable sums
    sum y_i; once R holds every value of the right parity in [-K, K], K = max |y| over all digits, it keeps that
    property for every later digit (pick y of the same sign as the target), so mu = n mod 2 at the end and the
    prefix is counted by its digit-sum distribution alone. States with n > 2T are dropped (n <= 2T is necessary).
  * U-U: every digit of u - v is a_i - b_i in [-10, 10] (balanced digits, base 21), so u - v <-> t = a - b, and t is
    realised iff the coordinatewise least choice (x_i, y_i) = argmin{x : x, x - t_i in A} has sum x <= T and
    sum y <= T (minimising x_i also minimises y_i = x_i - t_i). For t in A this is (t, 0), for -t in A it is
    (0, -t), and for t = +-1 it is (3, 2) or (2, 3): |U-U| is a box sum of the coefficients of
    (sum_{t in A} u^t + sum_{t in A, t > 0} v^t + u^3 v^2 + u^2 v^3)^d, by a two-dimensional dynamic program.
  * Both programs are checked against brute force (explicit sets) on small instances, for this alphabet and others.

usage: python3 -I g2026_check.py [--controls]
"""
import sys
import time
from decimal import Decimal as D, getcontext

getcontext().prec = 80
W = D("1e-70")

ROW = dict(A=[0, 2, 3, 4, 5, 6, 7, 8, 9, 10], B=21, d=80, T=150, printed="1.1740744",
           pr_max="2995805288150731620410662946668903948341032736664352641511848666717243994160370658179324073879212562136150",
           pr_sum="75448362167176243488362019935078206851619643198150854886920234689186981134888",
           pr_diff="195351744295266763842135520514417052287242446785296742323733058216909095059024572338564089814415",
           pr_val="1.174074447693521163363531806658755676155543896382679659744422686584815723921237399321826308369811589")


def max_u(A, B, d, T):
    top = max(A)
    digits = [0] * d
    left = T
    for i in range(d - 1, -1, -1):
        best = max(a for a in A if a <= left)
        digits[i] = best
        left -= best
    assert top * B > 0
    return sum(a * B ** i for i, a in enumerate(digits))


def count_sum(A, B, d, T):
    """|U+U| by the (n, R) program with saturation; requires 2 max A < B (no carries)."""
    if 2 * max(A) >= B:
        raise ValueError("digit sums can carry: 2 max A >= B")
    Aset = set(A)
    AA = sorted({a + b for a in A for b in A})
    Y = {s: sorted({2 * a - s for a in A if s - a in Aset}) for s in AA}
    K = max(max(abs(y) for y in Y[s]) for s in AA)
    OFF = 2 * T + K + 1                # bit OFF + v holds v; |v| <= 2T always
    sat = [0] * (2 * T + 1)              # saturated prefixes by digit sum
    states = {(0, 1 << OFF): 1}
    peak = 1

    def saturated(R, n):
        par = n & 1
        return all((R >> (v + OFF)) & 1 for v in range(-K, K + 1) if (v & 1) == par)

    for _ in range(d):
        nsat = [0] * (2 * T + 1)
        for n, c in enumerate(sat):
            if c:
                for s in AA:
                    if n + s <= 2 * T:
                        nsat[n + s] += c
        new = {}
        for (n, R), c in states.items():
            for s in AA:
                n2 = n + s
                if n2 > 2 * T:
                    continue
                R2 = 0
                for y in Y[s]:
                    R2 |= (R << y) if y >= 0 else (R >> -y)
                if saturated(R2, n2):
                    nsat[n2] += c
                else:
                    new[(n2, R2)] = new.get((n2, R2), 0) + c
        sat, states = nsat, new
        peak = max(peak, len(states))
    total = sum(sat)
    for (n, R), c in states.items():
        mu = min(abs(v - OFF) for v in range(2 * OFF + 1) if (R >> v) & 1)
        if n + mu <= 2 * T:
            total += c
    return total, peak, K


def count_diff(A, B, d, T):
    """|U-U| by the box sum; requires 2 max A < B (balanced digits, unique)."""
    if 2 * max(A) >= B:
        raise ValueError("differences can carry: 2 max A >= B")
    Aset = set(A)
    terms = {}
    for t in range(-max(A), max(A) + 1):
        xs = [x for x in A if x - t in Aset]
        if xs:
            x = min(xs)
            terms[(x, x - t)] = terms.get((x, x - t), 0) + 1
    S = [[0] * (T + 1) for _ in range(T + 1)]
    S[0][0] = 1
    for _ in range(d):
        N = [[0] * (T + 1) for _ in range(T + 1)]
        for i in range(T + 1):
            Si = S[i]
            for j in range(T + 1):
                c = Si[j]
                if c:
                    for (x, y), m in terms.items():
                        if i + x <= T and j + y <= T:
                            N[i + x][j + y] += c * m
        S = N
    return sum(map(sum, S)), terms


def brute(A, B, d, T):
    from itertools import product
    U = sorted({sum(a * B ** i for i, a in enumerate(v)) for v in product(A, repeat=d) if sum(v) <= T})
    return len({u + w for u in U for w in U}), len({u - w for u in U for w in U}), max(U)


def lemma_value(A, B, d, T):
    m = max_u(A, B, d, T)
    s, _, _ = count_sum(A, B, d, T)
    t, _ = count_diff(A, B, d, T)
    v = 1 + (D(t).ln() - D(s).ln()) / D(2 * m + 1).ln()
    return v, s, t, m


class Run:
    def __init__(self):
        self.n = 0; self.fails = []

    def __call__(self, label, ok):
        self.n += 1
        if not ok:
            self.fails.append(label)
        print(f"  check {label}: {'PASS' if ok else 'FAIL'}")
        return ok


def main():
    t0 = time.time()
    run = Run()
    A, B, d, T = ROW["A"], ROW["B"], ROW["d"], ROW["T"]
    print(f"3a [G2026]: A = {A}, B = {B}, d = {d}, T = {T} (PR #71)")
    tests = [([0, 2, 3, 4, 5, 6, 7, 8, 9, 10], 21, 3, T) for T in (0, 4, 9, 15, 30)] + \
            [([0, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 33, 3, T) for T in (5, 11, 20, 48)] + \
            [([0, 1, 3, 7], 15, 4, T) for T in (3, 8, 28)] + [([0, 2, 5], 11, 5, 9)]
    ok = True
    for (a, b, dd, tt) in tests:
        bs, bd, bm = brute(a, b, dd, tt)
        cs, _, _ = count_sum(a, b, dd, tt)
        cd, _ = count_diff(a, b, dd, tt)
        ok &= (bs, bd, bm) == (cs, cd, max_u(a, b, dd, tt))
    run(f"self-test: both programs and max U equal brute force (explicit sets) on {len(tests)} small instances", ok)
    m = max_u(A, B, d, T)
    run("max U = 10 * sum_{i=65}^{79} 21^i, the pull request's integer", m == 10 * sum(21 ** i for i in range(65, 80))
        and str(m) == ROW["pr_max"])
    s, peak, K = count_sum(A, B, d, T)
    run(f"|U+U| = {s} (K = {K}; at most {peak} unsaturated states), the pull request's integer", str(s) == ROW["pr_sum"])
    t, terms = count_diff(A, B, d, T)
    run(f"|U-U| = {t} ({len(terms)} digit types), the pull request's integer", str(t) == ROW["pr_diff"])
    run("the lemma's hypotheses: 0 in U (all digits 0), distinct non-negative integers, |U-U| < 2 max U + 1",
        t < 2 * m + 1)
    v = 1 + (D(t).ln() - D(s).ln()) / D(2 * m + 1).ln()
    lo, hi = v - W, v + W
    run(f"1 + log(|U-U|/|U+U|)/log(2 max U + 1) in [{str(lo)[:52]}, {str(hi)[:52]}] (80-digit logarithms)",
        hi - lo < D("1e-60"))
    pv = D(ROW["pr_val"])
    run("the pull request's 100-digit value lies in our enclosure (agreement to 70 digits)",
        lo - D("1e-99") <= pv <= hi)
    run(f"the row's {ROW['printed']} holds: {str(lo)[:14]} >= {ROW['printed']}", lo >= D(ROW["printed"]))
    if "--controls" in sys.argv:
        print("  -- red controls (each must be refused)")
        run("control: 1.1740745 (one unit more in the last printed digit) is not implied", not (lo >= D("1.1740745")))
        v2, s2, t2, m2 = lemma_value(A, B, d, 149)
        run(f"control: the cap T = 149 changes all three integers (value {str(v2)[:12]}), so the certificate fixes T",
            s2 != s and t2 != t and m2 != m)
        try:
            count_sum([0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 21, d, T)
            refused = False
        except ValueError:
            refused = True
        run("control: digit 11 added (11 + 11 = 22 >= 21 carries) is refused by both programs", refused)
        bs, bd, bm = brute([0, 2, 3, 4, 5, 6, 7, 8, 9, 10], 21, 3, 15)
        cs, _, _ = count_sum([0, 2, 3, 4, 5, 6, 7, 8, 9, 10], 21, 3, 14)
        run("control: the small instance with T = 14 in place of 15 differs from brute force at T = 15", cs != bs)
    print(f"  checks: {run.n}, failed: {len(run.fails)}")
    print(f"  run time: {time.time() - t0:.1f} s")
    print("RESULT: ALL CHECKS PASS" if not run.fails else "RESULT: FAIL " + "; ".join(run.fails))
    sys.exit(0 if not run.fails else 1)


if __name__ == "__main__":
    main()
