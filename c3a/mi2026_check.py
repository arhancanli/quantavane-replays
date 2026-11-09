"""3a: replay of the row [MI2026] C_3a >= 1.1835129324 (Mosaic Intelligence, PR #95), the largest 3a value certified by
exact finite counting, from the construction alone: A = {0,3,4,6,7,...,16}, base B = 33, d = 420 digits, digit sum at
most T = 1392 (the parameters of the certificate's counts file, gh3a_b33_d420_T1392_counts.json, extracted by itself
from the Zenodo archive 10.5281/zenodo.20794135 and read as data; the archive's two programs were neither read nor
run). Its three integers are compared with ours only after our computation (Quantavane replays).

Our counts use the same programs as g2026_check.py (loaded from this folder): no carries anywhere (2 * 16 < 33);
|U+U| by the dynamic program on (digit sum, reachable signed split sums) with the saturation lemma (here K = 16);
|U-U| by the coordinatewise least split ((t, 0), (0, -t), and (4,3), (6,4), (8,3) for t = 1, 2, 5 with their
mirrors) and a box sum over [0, 1392]^2 of the coefficients of the 420th power, computed by diff2d.c modulo primes
below 2^61 whose product exceeds the a priori bound 33^420 >= |U-U| (35 primes), put together by the Chinese
remainder theorem; two more primes check the result. diff2d is first checked against the exact Python box sum on
small instances.

usage: python3 -I mi2026_check.py CERTS_ROOT BINDIR [--controls]
"""
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
from decimal import Decimal as D

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("g2026_check", os.path.join(HERE, "g2026_check.py"))
G = importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)
COUNTS_SHA = "7acecfea1298e444c7d1045f0019d1bdfcc43ecca67f4b8b174f4e26655545f1"


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    dd, s = n - 1, 0
    while dd % 2 == 0:
        dd //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):   # deterministic below 3.3e24
        x = pow(a, dd, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_below(top, k):
    out, n = [], top - 1
    while len(out) < k:
        if is_prime(n):
            out.append(n)
        n -= 1
    return out


def diff_terms(A):
    Aset = set(A)
    terms = {}
    for t in range(-max(A), max(A) + 1):
        xs = [x for x in A if x - t in Aset]
        if xs:
            x = min(xs)
            terms[(x, x - t)] = terms.get((x, x - t), 0) + 1
    return terms


def diff_mod(bindir, A, d, T, primes):
    terms = diff_terms(A)
    inp = "\n".join(f"{x} {y} {m}" for (x, y), m in sorted(terms.items())) + "\n"
    run = subprocess.run([os.path.join(bindir, "diff2d"), str(d), str(T)] + [str(p) for p in primes],
                         input=inp, capture_output=True, text=True, check=True)
    res = {}
    for line in run.stdout.split("\n"):
        if line.strip():
            p, r = line.split()
            res[int(p)] = int(r)
    return res


def crt(res):
    x, m = 0, 1
    for p, r in res.items():
        t = ((r - x) * pow(m, -1, p)) % p
        x, m = x + m * t, m * p
    return x, m


class Run:
    def __init__(self):
        self.n = 0; self.fails = []

    def __call__(self, label, ok):
        self.n += 1
        if not ok:
            self.fails.append(label)
        print(f"  check {label}: {'PASS' if ok else 'FAIL'}", flush=True)
        return ok


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    if len(args) != 2:
        print(__doc__); sys.exit(2)
    root, bindir = args
    t0 = time.time()
    run = Run()
    raw = open(os.path.join(root, "mi3a", "gh3a_b33_d420_T1392_counts.json"), "rb").read()
    run("counts file sha256 = 7acecfea...45f1 (c3a_pr95/gh3a_b33_d420_T1392_counts.json of Zenodo 20794135)",
        hashlib.sha256(raw).hexdigest() == COUNTS_SHA)
    cert = json.loads(raw)
    con = cert["construction"]
    A, B, d, T = [int(a) for a in con["digit_set"]], int(con["base"]), int(con["digit_count"]), int(con["sum_cap"])
    print(f"3a [MI2026]: A = {A}, B = {B}, d = {d}, T = {T} (PR #95)", flush=True)
    run("the construction as the 3a page describes it: base 33, digits in A, digit sum capped (2 max A < B: no carries)",
        B == 33 and 2 * max(A) < B and 0 in A)
    # diff2d against the exact Python box sum on small instances
    ok = True
    nterms = len(diff_terms(A))
    bound = nterms ** d                    # |U-U| <= (number of digit differences)^d
    ps = []
    mod_needed = 1
    for p in primes_below(1 << 61, 80):
        if mod_needed > bound:
            break
        ps.append(p); mod_needed *= p
    nmain = len(ps)
    ps += primes_below(min(ps) - 1, 2)      # two check primes
    for (a, b, dd, tt) in [(A, B, 6, 30), (A, B, 9, 41), ([0, 2, 3, 4, 5, 6, 7, 8, 9, 10], 21, 7, 25), ([0, 1, 3, 7], 15, 8, 20)]:
        exact, _ = G.count_diff(a, b, dd, tt)
        res = diff_mod(bindir, a, dd, tt, ps[:3])
        ok &= all(r == exact % p for p, r in res.items())
    run("diff2d equals the exact Python box sum on four small instances (three primes each)", ok)
    m = G.max_u(A, B, d, T)
    run(f"max U by the greedy digits: 87 top digits 16 (T = 1392 = 87 * 16)",
        m == 16 * sum(B ** i for i in range(d - 87, d)))
    t1 = time.time()
    res = diff_mod(bindir, A, d, T, ps)
    main_res = {p: res[p] for p in ps[:nmain]}
    D_, mod = crt(main_res)
    check = all(D_ % p == res[p] for p in ps[nmain:])
    run(f"|U-U| = {D_} (CRT over {nmain} primes below 2^61, product {mod.bit_length()} bits > {nterms}^{d} >= |U-U|; "
        f"two more primes agree; {time.time() - t1:.0f} s)", check and mod > bound)
    t2 = time.time()
    S_, peak, K = G.count_sum(A, B, d, T)
    run(f"|U+U| = {S_} (K = {K}; at most {peak} unsaturated states; {time.time() - t2:.0f} s)", S_ > 0)
    print("  -- only now compared with the certificate's integers", flush=True)
    run("max U equals the certificate's max_u", str(m) == str(cert["max_u"]))
    run("|U+U| equals the certificate's sum_count", str(S_) == str(cert["sum_count"]))
    run("|U-U| equals the certificate's diff_count", str(D_) == str(cert["diff_count"]))
    run("the lemma's hypotheses: 0 in U, |U-U| < 2 max U + 1", D_ < 2 * m + 1)
    v = 1 + (D(D_).ln() - D(S_).ln()) / D(2 * m + 1).ln()
    lo, hi = v - G.W, v + G.W
    run(f"1 + log(|U-U|/|U+U|)/log(2 max U + 1) in [{str(lo)[:50]}, {str(hi)[:50]}]", hi - lo < D("1e-60"))
    run("the certificate's 40-digit floor 1.1835129324218615106564747894020784326947 lies in our enclosure's floor",
        D(cert["value_lo_40"]) <= lo and lo - D(cert["value_lo_40"]) < D("1e-40"))
    run("the row's 1.1835129324 holds", lo >= D("1.1835129324"))
    if "--controls" in sys.argv:
        print("  -- red controls (each must be refused)", flush=True)
        run("control: 1.1835129325 (one unit more in the row's last digit) is not implied", not (lo >= D("1.1835129325")))
        r2 = diff_mod(bindir, A, d, T - 1, ps[:1])
        run("control: the cap T = 1391 gives a different |U-U| (residue modulo one prime)", r2[ps[0]] != D_ % ps[0])
        try:
            G.count_sum(A + [17], B, d, T); refused = False
        except ValueError:
            refused = True
        run("control: a digit 17 added (17 + 17 = 34 >= 33 carries) is refused", refused)
    print(f"  checks: {run.n}, failed: {len(run.fails)}")
    print(f"  run time: {time.time() - t0:.0f} s")
    print("RESULT: ALL CHECKS PASS" if not run.fails else "RESULT: FAIL " + "; ".join(run.fails))
    sys.exit(0 if not run.fails else 1)


if __name__ == "__main__":
    main()
