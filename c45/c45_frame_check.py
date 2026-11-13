#!/usr/bin/env python3
"""Replays of the C45 rows with bases in the frame T0 = 12252240 (PR #290 [S2026], 54 primes; PR #281 [P2026b],
[P2026c], 57 and 59 primes), from the obstruction inequality on the 45a page and the prime sets in the rows. Our own
code (c45_frame.c and this file, with the group laws and cluster distributions of c45_check.py); no program of any
certificate was read or run. The claimant's histogram file, if present, is read as data, only to compare it with ours.

Steps. (1) Self-tests: the pipeline against brute force from the definition (every a mod M0) on small base sets that
use every group. (2) The row's base: c45_frame computes the histogram of nu0 for nu0 <= C exactly, the weight of the
classes with nu0 > C, and the exact first and second moments of nu0 over all classes; we check the total (= M0) and
both moments against the closed forms sum nu0 = T0 prod (q - 1) and sum nu0^2 = T0 sum_d prod_q (q - 1 or q - 2).
(3) The cut: (C + 1) kappa H_min >= 1/2, so every class above C contributes 1/(2M) for every cluster multiplier.
(4) S = sum_nu delta(nu) E_H[min(1/2, kappa nu H)] in exact rationals (a sorted sweep over the law of H), with our
lambda = sum_{k <= 120} 1/(k 2^k) < log 2, and also with lambda = 2 sum_{j < 40} 1/((2j+1) 3^(2j+1)) < log 2 to compare
with a printed S.

Usage: python3 -I c45_frame_check.py selftest BINDIR
       python3 -I c45_frame_check.py ROW BINDIR STATEDIR [--threads N] [--passes N] [--certs DIR] [--controls]
       (threads 2 and passes 4 by default; a finished or stopped run in STATEDIR with the same settings is resumed)
       rows: s2026 (PR #290), p2026b, p2026c (PR #281)
"""
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import time
from fractions import Fraction
from math import gcd, prod

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from c45_check import order2, is_prime, group_vectors, cluster_dist, product_dist, lam, dec   # noqa: E402

T0 = 12252240
A_SLOTS = (("G1", 12), ("G9", 36), ("G5", 60), ("G16", 48))      # see t = r mod 12
S_SLOTS = (("G7", 28), ("G13", 52), ("G11", 44), ("G17", 68))    # orders prime to 3: see only s = r mod 4

CL290 = {23: [47, 178481], 29: [233, 1103, 2089], 37: [223], 43: [431, 9719], 47: [2351, 4513], 53: [6361],
         73: [439], 79: [2687], 83: [167], 113: [3391], 131: [263], 179: [359, 1433], 191: [383], 233: [1399],
         239: [479, 1913, 5737, 176383], 251: [503]}
BASE290 = [3, 5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 61, 73, 257, 43, 127, 241, 113, 109, 97, 151,
           8191, 53, 157, 23, 89, 397, 683, 137]
ROWS = {
    "s2026": dict(base=BASE290, clusters=CL290, printed="0.489481743486002", cut=1538836,
                  where="PR #290, row [S2026] (29 base primes and 25 cluster primes)",
                  s_printed="0.48948174348600104752999373053459010942",
                  hist=("c45-54/seed_histogram.txt.gz", "e011429b8f49cd362db1f8a3273fe5ceeff97d1067af1e97bfbbdf2bddfaf91c")),
    "p2026b": dict(base=BASE290 + [953, 1613, 2113], clusters=CL290, printed="0.489473818147233", cut=1535541,
                   where="PR #281, row [P2026b] (32 base primes and 25 cluster primes)",
                   s_printed="0.48947381814723275919828409928792897470", hist=None),
    "p2026c": dict(base=BASE290 + [953, 1613, 2113, 331, 673], clusters=CL290, printed="0.489440909905826",
                   cut=1528627, where="PR #281, row [P2026c] (34 base primes and 25 cluster primes)",
                   s_printed="0.48944090990582567951473617427079243852", hist=None),
}


def assign(base):
    groups = {name: [] for name, _ in A_SLOTS + S_SLOTS}
    for q in base:
        m = order2(q)
        for name, L in A_SLOTS + S_SLOTS:
            if L % m == 0:
                groups[name].append(q)
                break
        else:
            raise ValueError(f"ord_{q}(2) = {m} divides none of 12, 36, 60, 48, 28, 52, 44, 68")
    return groups


def frame_input(base):
    groups = assign(base)
    lines = []
    for name, L in A_SLOTS:
        d = group_vectors(groups[name], L)
        lines.append(f"{name} {len(d)}")
        lines += [" ".join(map(str, v)) + f" {w}" for v, w in sorted(d.items())]
    lines.append(f"NB {len(S_SLOTS)}")
    for name, L in S_SLOTS:
        d = group_vectors(groups[name], L)
        lines.append(f"{name} {len(d)}")
        lines += [" ".join(map(str, v)) + f" {w}" for v, w in sorted(d.items())]
    return "\n".join(lines) + "\n", groups


def run_frame(base, cut, bindir, statedir, threads=1, sig=200, log=print, passes=1, stop_after=None):
    text, groups = frame_input(base)
    os.makedirs(statedir, exist_ok=True)
    t0 = time.time()
    cmd = [os.path.join(bindir, "c45_frame"), str(cut), str(threads), statedir, str(sig), str(passes)]
    if stop_after is not None:
        cmd.append(str(stop_after))
    run = subprocess.run(cmd, input=text, capture_output=True, text=True)
    if stop_after is not None and run.returncode == 9:
        return None
    if run.returncode not in (0, 1):
        raise RuntimeError(f"c45_frame failed ({run.returncode}): {run.stderr[-2000:]}")
    hist, extra = {}, {}
    for line in run.stdout.split("\n"):
        if not line.strip():
            continue
        if line.startswith("#"):
            k, v = line[1:].split()
            extra[k] = int(v)
        else:
            nu, c = line.split()
            hist[int(nu)] = int(c)
    for line in run.stderr.strip().split("\n"):
        if line.startswith(("B law", "canonical", "stream done", "done", "resumed")):
            log(f"  c45_frame: {line}")
    return hist, extra, groups, time.time() - t0


def closed_moments(base):
    """sum nu0 and sum nu0^2 over a mod M0, from the definition (frame T0): sum over d of prod_q (q - 1 or q - 2)"""
    ms = [order2(q) for q in base]
    m1 = T0 * prod(q - 1 for q in base)
    # sum_{d mod T0} f(gcd(d, T0)) = sum_{g | T0} f(g) phi(T0 / g)
    divs = [g for g in range(1, T0 + 1) if T0 % g == 0]

    def phi(n):
        r, x, p = n, n, 2
        while p * p <= x:
            if x % p == 0:
                while x % p == 0:
                    x //= p
                r -= r // p
            p += 1
        if x > 1:
            r -= r // x
        return r
    m2 = T0 * sum(phi(T0 // g) * prod((q - 1) if g % m == 0 else (q - 2) for q, m in zip(base, ms)) for g in divs)
    return m1, m2


def brute_hist(base):
    """from the definition: every a mod M0, nu(a) = #{r mod T : gcd(a - 2^r, M0) = 1}, scaled to Z/T0"""
    M0 = prod(base)
    T = 1
    for q in base:
        T = T * order2(q) // gcd(T, order2(q))
    pw = [pow(2, r, M0) for r in range(T)]
    hist = {}
    for a in range(M0):
        nu = sum(1 for r in range(T) if gcd((a - pw[r]) % M0, M0) == 1) * (T0 // T)
        hist[nu] = hist.get(nu, 0) + 1
    return hist


def evaluate(hist, above, base, clusters, lam_value, cut):
    """exact S = E[min(1/2, kappa nu0 H)] with every class above the cut capped (checked); a sorted sweep over H"""
    allq = base + [q for qs in clusters.values() for q in qs]
    M0 = prod(base)
    T = T0 * prod(clusters)
    kappa = Fraction(prod(allq), prod(q - 1 for q in allq)) / (T * lam_value)
    P, Q = kappa.numerator, kappa.denominator
    H = product_dist([cluster_dist(c, qs) for c, qs in sorted(clusters.items())])
    hs = sorted(H.items())
    WH = sum(H.values())
    hmin = hs[0][0]
    cut_ok = 2 * P * (cut + 1) * hmin >= Q          # every nu > cut: kappa nu h >= 1/2 for every h
    W, X = [0], [0]
    for h, w in hs:
        W.append(W[-1] + w)
        X.append(X[-1] + w * h)
    num = 0
    k = len(hs)
    for nu in sorted(hist):
        c = hist[nu]
        if nu == 0:
            continue                                # the term is 0
        while k > 0 and 2 * P * nu * hs[k - 1][0] >= Q:   # capped: kappa nu h >= 1/2
            k -= 1
        num += c * (2 * P * nu * X[k] + Q * (WH - W[k]))
    num += above * Q * WH
    return Fraction(num, 2 * Q * M0 * WH), cut_ok, len(hs)


def lam_atanh():
    """2 sum_{j < 40} 1/((2j+1) 3^(2j+1)), a partial sum of 2 atanh(1/3) = log 2"""
    return 2 * sum(Fraction(1, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(40))


def read_their_hist(path, sha):
    raw = gzip.open(path, "rb").read() if path.endswith(".gz") else open(path, "rb").read()
    got = hashlib.sha256(raw).hexdigest()
    h, agg = {}, None
    for line in raw.decode().split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            parts = line.lstrip("#").split()
            nums = [int(x) for x in parts if x.isdigit()]
            agg = nums
            continue
        nu, c = line.split()
        h[int(nu)] = int(c)
    return h, agg, got == sha, got


def selftest(bindir, statedir):
    ok = True
    print("Self-test (frame T0 = 12252240): c45_frame against brute force over every a mod M0, from the definition")
    # G1's four primes cover Z/12 densely, so nu0 takes values below the cut; each set adds other groups
    sets = [[3, 5, 7, 13, 23], [3, 5, 7, 13, 53], [3, 5, 7, 13, 29], [3, 5, 7, 13, 19], [3, 5, 7, 13, 31],
            [3, 5, 7, 13, 17], [3, 5, 7, 13, 97], [3, 5, 13, 683], [3, 13, 23, 29], [3, 13, 53, 17], [3, 7, 13, 137]]
    for i, base in enumerate(sets):
        t0 = time.time()
        hb = brute_hist(base)
        vals = sorted(hb)
        cut = min(max(vals[len(vals) // 2], 1), 6000000)      # about half the values below the cut
        h, extra, groups, _ = run_frame(base, cut, bindir, os.path.join(statedir, f"self{i}"), threads=2, log=lambda s: None)
        shutil.rmtree(os.path.join(statedir, f"self{i}"))         # its saved state holds C + 1 bins of 32 bytes
        below = {nu: c for nu, c in hb.items() if nu <= cut}
        above = sum(c for nu, c in hb.items() if nu > cut)
        m1, m2 = closed_moments(base)
        good = len(below) > 0 and h == below and extra["above"] == above and extra["below"] + extra["above"] == prod(base) \
            and extra["moment1"] == m1 and extra["moment2"] == m2 and extra["partition_failures"] == 0
        used = ", ".join(f"{k} {v}" for k, v in groups.items() if v)
        print(f"  base {base} ({used}; M0 = {prod(base)}): {len(h)} bins <= C ({len(hb)} values in all) and the weight above C as brute force; "
              f"moments exact: {'identical' if good else 'DIFFERENT'} ({time.time() - t0:.0f} s)")
        ok &= good
    # the other code paths: every V in the flushed table (SIG = 2), and a resumed run from a saved state
    base = sets[1]
    hb = brute_hist(base)
    vals = sorted(hb)
    cut = vals[len(vals) // 2]
    ref = run_frame(base, cut, bindir, os.path.join(statedir, "path_a"), threads=3, sig=1000, log=lambda s: None)[:2]
    alt = run_frame(base, cut, bindir, os.path.join(statedir, "path_b"), threads=2, sig=2, log=lambda s: None)[:2]
    again = run_frame(base, cut, bindir, os.path.join(statedir, "path_a"), threads=1, sig=1000, log=lambda s: None)[:2]
    split = run_frame(base, cut, bindir, os.path.join(statedir, "path_c"), threads=2, sig=1000, passes=4, log=lambda s: None)[:2]
    stopped = run_frame(base, cut, bindir, os.path.join(statedir, "path_d"), threads=2, sig=1000, passes=8, stop_after=2, log=lambda s: None)
    resumed = run_frame(base, cut, bindir, os.path.join(statedir, "path_d"), threads=2, sig=1000, passes=8, log=lambda s: None)[:2]
    for sub in ("path_a", "path_b", "path_c", "path_d"):
        shutil.rmtree(os.path.join(statedir, sub))
    good = ref == alt == again == split == resumed and stopped is None and len(ref[0]) > 0
    print(f"  base {base}: every V through the flushed table (SIG = 2), a run resumed from the saved state of a finished "
          f"run, the low vectors split into 4 passes, and an 8-pass run stopped after pass 2 and resumed give the same "
          f"histogram and totals: {good}")
    ok &= good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}")
    return ok


def check_row(key, bindir, statedir, threads, certs, controls, passes=4):
    row = ROWS[key]
    t0 = time.time()
    ok = True
    log = print
    base, clusters = row["base"], row["clusters"]
    for q in base:
        if not (is_prime(q) and T0 % order2(q) == 0):
            raise ValueError(f"{q}: not a prime with order dividing T0")
    for c, qs in clusters.items():
        for q in qs:
            if order2(q) != c or gcd(c, T0) != 1 or not is_prime(c):
                raise ValueError(f"cluster {c}: {q} has order {order2(q)}")
    log(f"{row['where']}: C45 < {row['printed']}")
    groups = assign(base)
    log("  base groups: " + "; ".join(f"{k} {v}" for k, v in groups.items()))
    hist, extra, _, secs = run_frame(base, row["cut"], bindir, statedir, threads=threads, passes=passes)
    M0 = prod(base)
    m1, m2 = closed_moments(base)
    good = extra["below"] + extra["above"] == M0 and sum(hist.values()) == extra["below"] \
        and extra["moment1"] == m1 and extra["moment2"] == m2 and extra["partition_failures"] == 0
    log(f"  histogram: {len(hist):,} bins with nu0 <= {row['cut']:,} ({extra['pairs']:,} (V, B) pairs, {secs:.0f} s); "
        f"total weight = M0, sum nu0 = T0 prod (q - 1), sum nu0^2 = T0 sum_d prod_q (q - 1 or q - 2), every V's "
        f"pairs partitioned: {'all exact' if good else 'MISMATCH'}")
    ok &= good
    if certs and row["hist"]:
        path = os.path.join(certs, row["hist"][0])
        if os.path.exists(path):
            theirs, agg, same_sha, got = read_their_hist(path, row["hist"][1])
            same = theirs == hist
            log(f"  the claimant's {row['hist'][0]} (sha256 of the text {got[:8]}..., as stated: {same_sha}; "
                f"{len(theirs):,} bins), read as data: bins identical to ours: {same}")
            if agg and len(agg) == 4:
                # '# above C count sum_nu sum_nu^2': ours above C = closed-form totals minus our bins
                s1 = m1 - sum(nu * c for nu, c in hist.items())
                s2 = m2 - sum(nu * nu * c for nu, c in hist.items())
                good = agg[0] == row["cut"] and agg[1:] == [extra["above"], s1, s2]
                log(f"    its aggregate line above C (count, sum nu0, sum nu0^2) equals ours (closed-form totals minus our "
                    f"bins): {good}")
                ok &= good
            ok &= same and same_sha
    lv = lam()
    S, cut_ok, nh = evaluate(hist, extra["above"], base, clusters, lv, row["cut"])
    log(f"  cut: (C + 1) kappa H_min >= 1/2, so every class above C gives 1/(2M) for every cluster multiplier: {cut_ok}")
    ok &= cut_ok
    b = Fraction(row["printed"])
    nd = len(row["printed"].split(".")[1])
    good = S <= b and S > b - Fraction(1, 10 ** nd)
    log(f"  lambda = sum_(k <= 120) 1/(k 2^k) < log 2; {nh:,} distinct cluster multipliers H")
    log(f"  S = {dec(S, 40)}... (exact rational)")
    log(f"  S <= {row['printed']}, and the printed value is S rounded up at its last digit: {good}")
    ok &= good
    S2, _, _ = evaluate(hist, extra["above"], base, clusters, lam_atanh(), row["cut"])
    pr = row["s_printed"]
    same = dec(S2, len(pr) - 2) == pr
    log(f"  with lambda = 2 sum_(j < 40) 1/((2j+1) 3^(2j+1)): S = {dec(S2, 45)}...; the printed {pr}... agrees in "
        f"every printed digit: {same}")
    ok &= same
    if controls:
        log("Red controls:")
        c0 = max(clusters, key=lambda c: len(clusters[c]))
        cl = {c: (qs[:-1] if c == c0 else qs) for c, qs in clusters.items()}
        Sc, _, _ = evaluate(hist, extra["above"], base, cl, lv, row["cut"])
        good = Sc > b
        log(f"  without the cluster prime {clusters[c0][-1]} (order {c0}): S = {dec(Sc, 15)}..., above {row['printed']}: "
            f"{good}")
        ok &= good
        try:
            wrong = {c: list(qs) for c, qs in clusters.items()}
            wrong[min(wrong)] = wrong[min(wrong)] + [53]
            evaluate(hist, extra["above"], base, wrong, lv, row["cut"])
            log("  53 placed in a cluster: NOT REFUSED")
            ok = False
        except ValueError:
            log(f"  53 (order 52) placed in the cluster of order {min(clusters)}: refused, as expected")
        try:
            assign(base + [337])           # ord_337(2) = 21 divides none of the eight group moduli
            log("  337 (order 21) added to the base: NOT REFUSED")
            ok = False
        except ValueError:
            log("  337 (order 21) added to the base: refused (21 divides none of the group moduli), as expected")
        # 14 Nov 2026: this control first used 2731, whose order 26 divides the modulus 52, so assign() accepted
        # it, rightly, and the run printed FAIL for a wrong control; 2731 is now the positive side of the control.
        if 2731 not in assign(base + [2731])["G13"]:
            log("  2731 (order 26, a divisor of 52) added to the base: not placed in the group of modulus 52")
            ok = False
        _, cut_low, _ = evaluate(hist, extra["above"], base, clusters, lv, row["cut"] - 2000)
        good = not cut_low
        log(f"  a cut 2000 lower: the classes above it are not all capped, refused: {good}")
        ok &= good
        bad = dict(hist)
        k = sorted(bad)[len(bad) // 2]
        bad[k] -= 1
        bad[k + 1] = bad.get(k + 1, 0) + 1
        Sb, _, _ = evaluate(bad, extra["above"], base, clusters, lv, row["cut"])
        good = Sb != S
        log(f"  one class moved from nu0 = {k} to {k + 1}: S changes ({dec(Sb, 30)}...): {good}")
        ok &= good
    log(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = sys.argv[1:]
    threads = int(opts[opts.index("--threads") + 1]) if "--threads" in opts else 2
    passes = int(opts[opts.index("--passes") + 1]) if "--passes" in opts else 4
    certs = opts[opts.index("--certs") + 1] if "--certs" in opts else None
    args = [a for a in args if a not in (str(threads), str(passes), certs)]
    if len(args) < 2:
        print(__doc__)
        return 2
    if args[0] == "selftest":
        return 0 if selftest(args[1], args[2] if len(args) > 2 else os.path.join(args[1], "state_selftest")) else 1
    if args[0] not in ROWS or len(args) < 3:
        print(__doc__)
        return 2
    return 0 if check_row(args[0], args[1], args[2], threads, certs, "--controls" in opts, passes) else 1


if __name__ == "__main__":
    sys.exit(main())
