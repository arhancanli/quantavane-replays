#!/usr/bin/env python3
"""Exact replay of weight-class step graphon bounds on the Ramsey multiplicity constants (89 pages: C_89c for K_6).

A row is a step graphon on the 2^n equal blocks of Z_2^n whose red density between blocks x and y is n_w / D with
w = |x xor y| (n_0 is the density inside a block). Its monochromatic K_t density is computed exactly:

  1. Two C programs written apart, c89_multiset.c (multisets of column types) and c89_refine.c (vertices placed one
     at a time on cells of coordinates, with pruning), count the (2^n)^(t-1) vertex tuples behind each red and blue
     monomial prod_k f_k^(e_k) and prod_k (1 - f_k)^(e_k); the two count tables must agree exactly and each total must
     be 2^(n(t-1)).
  2. The density sum_e [red_e prod f^e + blue_e prod (1-f)^e] / 2^(n(t-1)) is formed in exact rationals (fractions).
  3. It is compared with the fraction the row states, and the printed bound is tested exactly.

Usage:
  python3 -I c89_check.py selftest            brute force on small cases, f = 1/2, and Franek's 2002 value
  python3 -I c89_check.py ROW [CERTS] [--controls] [--workers=N]
Python 3.9+, standard library only, and a C compiler (cc).
"""
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from fractions import Fraction
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")

# Rows. "file" is the data file fetch.sh downloads (None when the PR text is the only public source of the data).
ROWS = {
    # PR #279 (JeetuSK0808), [SK2026] on the proposed 89c page: construction.json @ 98b6ee3
    "sk2026": dict(
        t=6, n=10, D=2 ** 20,
        nw=[1048576, 461327, 0, 1048576, 1048576, 13579, 0, 1048576, 1048576, 0, 799568],
        file=("sk2026/construction.json", "f90cc19090067c5a480257e802a3e47bd55cb75f6b71d3a545419e91cd11ac49"),
        claimed=Fraction(
            180244428715156181493112659521741405180580225001966676111887908798425648907356856409,
            3978585891278293137243057985174566720803649206378781739523711815145275976100267004264448),
        printed="0.7422548616", scale=Fraction(1, 2 ** 14), pr="#279"),
    # PR #331 (Doldrums), [SCIEL2026]: the eleven n_w as stated in the PR and on its 89c page @ b77711a; the Zenodo
    # record 10.5281/zenodo.23974391 that holds construction.json is restricted, so no file is fetched for it
    "sciel2026": dict(
        t=6, n=10, D=2 ** 32,
        nw=[2 ** 32, 1917548695, 0, 2 ** 32, 2 ** 32, 0, 81376308, 2 ** 32, 2 ** 32, 0, 0],
        file=None,
        claimed=Fraction(
            67432845284234356713012367110897548622422530589772294031336898680659579871445572986246952507783885816983659288359622409120016266716137,
            1488565707357402911845015158554633286356257506687627387456491927921949262056238946972039271861787782268441644249633132407878864135402815488),
        printed="0.7422042116", scale=Fraction(1, 2 ** 14), pr="#331"),
}
# Franek (Combinatorica 22, 2002): red iff |x xor y| in {1, 3, 4, 7, 8, 10}, red inside the blocks; 51162598917 / 2^50
FRANEK = dict(t=6, n=10, D=1, nw=[1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1], value=Fraction(51162598917, 2 ** 50))


def build():
    os.makedirs(BIN, exist_ok=True)
    for name in ("c89_multiset", "c89_refine"):
        src, exe = os.path.join(HERE, name + ".c"), os.path.join(BIN, name)
        if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
            subprocess.run(["cc", "-O2", "-std=c99", "-o", exe, src], check=True)


def classes(nw, D):
    """R / B / fractional index per weight; refuses anything outside [0, D]."""
    if D <= 0 or any(not isinstance(x, int) for x in nw + [D]):
        raise ValueError("densities must be integers n_w over a positive integer D")
    codes, fracs = [], []
    for w, x in enumerate(nw):
        if x < 0 or x > D:
            raise ValueError(f"n_{w} = {x} is outside [0, {D}]")
        if x == D:
            codes.append("R")
        elif x == 0:
            codes.append("B")
        else:
            if x not in fracs:
                fracs.append(x)
            codes.append(str(fracs.index(x)))
    if len(fracs) > 4:
        raise ValueError("at most four fractional densities")
    return ",".join(codes), [Fraction(x, D) for x in fracs]


def parse(out):
    tab, total = {}, None
    for line in out.split("\n"):
        p = line.split()
        if not p:
            continue
        if p[0] == "TOTAL":
            total = int(p[1])
        else:
            key = (p[0], tuple(int(x) for x in p[1:5]))
            tab[key] = tab.get(key, 0) + int(p[5])
    return tab, total


def run_multiset(t, n, codes, workers):
    procs = [subprocess.Popen([os.path.join(BIN, "c89_multiset"), str(t), str(n), codes, str(i), str(workers)],
                              stdout=subprocess.PIPE, text=True) for i in range(workers)]
    tab, total = {}, 0
    for p in procs:
        out, _ = p.communicate()
        if p.returncode:
            raise RuntimeError("c89_multiset failed")
        a, tt = parse(out)
        for k, v in a.items():
            tab[k] = tab.get(k, 0) + v
        total += tt
    return tab, total


def run_refine(t, n, codes):
    r = subprocess.run([os.path.join(BIN, "c89_refine"), str(t), str(n), codes], capture_output=True, text=True,
                       check=True)
    return parse(r.stdout)


def density(tab, fracs, t, n):
    s = Fraction(0)
    for (colour, ex), cnt in tab.items():
        term = Fraction(cnt)
        for k, e in enumerate(ex):
            if e:
                term *= (fracs[k] if colour == "R" else 1 - fracs[k]) ** e
        s += term
    return s / 2 ** (n * (t - 1))


def exact(t, n, D, nw, workers, log=print):
    """Both programs, tables compared, totals checked; returns the exact density."""
    codes, fracs = classes(nw, D)
    t0 = time.time()
    a, ta = run_multiset(t, n, codes, workers)
    t1 = time.time()
    b, tb = run_refine(t, n, codes)
    t2 = time.time()
    full = 2 ** (n * (t - 1))
    log(f"  classes {codes}; method A {len(a)} cells, total {ta} ({t1 - t0:.0f} s); "
        f"method B {len(b)} cells, total {tb} ({t2 - t1:.0f} s); 2^(n(t-1)) = {full}")
    if ta != full or tb != full:
        raise AssertionError("a total is not 2^(n(t-1))")
    if a != b:
        raise AssertionError("the two count tables differ")
    return density(a, fracs, t, n)


def brute(t, n, D, nw):
    """Direct enumeration of every tuple (small n only), same density by a third route."""
    f = [Fraction(x, D) for x in nw]
    pairs = [(i, j) for i in range(t) for j in range(i + 1, t)]
    s = Fraction(0)
    for xs in product(range(2 ** n), repeat=t - 1):
        v = (0,) + xs
        r = b = Fraction(1)
        for i, j in pairs:
            q = f[bin(v[i] ^ v[j]).count("1")]
            r *= q
            b *= 1 - q
        s += r + b
    return s / 2 ** (n * (t - 1))


def selftest(workers):
    ok = True
    rng = random.Random(20261105)
    print("1. Brute force on small cases (every tuple enumerated in Python) against both programs:")
    for t, n in ((4, 3), (4, 4), (5, 3), (6, 2), (6, 3)):
        for trial in range(2):
            D = 7
            pool = [0, D, rng.randint(1, D - 1), rng.randint(1, D - 1)]
            nw = [rng.choice(pool) for _ in range(n + 1)]
            v = exact(t, n, D, nw, workers, log=lambda s: None)
            w = brute(t, n, D, nw)
            good = v == w
            ok &= good
            print(f"  t={t} n={n} n_w/D={nw}/{D}: {'agree' if good else 'DIFFER'} ({float(v):.12g})")
    print("2. f = 1/2 everywhere gives 2^(1 - C(t,2)) exactly:")
    for t, n in ((6, 6), (5, 8)):
        v = exact(t, n, 2, [1] * (n + 1), workers, log=lambda s: None)
        good = v == Fraction(2, 2 ** (t * (t - 1) // 2))
        ok &= good
        print(f"  t={t} n={n}: {v} {'= 2^(1-C(t,2))' if good else 'WRONG'}")
    print("3. Franek 2002 (red iff |x xor y| in {1,3,4,7,8,10}, red inside blocks), stated 51162598917/2^50:")
    v = exact(FRANEK["t"], FRANEK["n"], FRANEK["D"], FRANEK["nw"], workers)
    good = v == FRANEK["value"]
    ok &= good
    print(f"  {v} = {float(v / Fraction(1, 2 ** 14)):.12f} * 2^-14: {'reproduced exactly' if good else 'DIFFERS'}")
    print("RESULT: selftest", "PASS" if ok else "FAIL")
    return ok


def check_row(key, root, controls, workers):
    r = ROWS[key]
    t0 = time.time()
    print(f"Row {key} ({r['pr']}): K_{r['t']} on Z_2^{r['n']}, n_w / {r['D']} = {r['nw']}")
    if r["file"]:
        path, sha = r["file"]
        data = open(os.path.join(root, path), "rb").read()
        h = hashlib.sha256(data).hexdigest()
        print(f"  {path}: sha256 {h} ({'as pinned' if h == sha else 'NOT AS PINNED'})")
        assert h == sha, "file differs from the pinned bytes"
        c = json.loads(data)
        assert c["D"] == r["D"] and c["n"] == r["nw"], "the file's D and n differ from the row"
        print("  the file's D and n are the row's")
    else:
        print("  data: the eleven n_w as stated in the pull request (its certificate record is restricted)")
    v = exact(r["t"], r["n"], r["D"], r["nw"], workers)
    print(f"  exact density = {v.numerator}/{v.denominator}")
    print(f"               = {float(v):.15e} = {float(v / r['scale']):.15f} * 2^-14")
    same = v == r["claimed"]
    print(f"  stated fraction: {'IDENTICAL' if same else 'DIFFERENT'}")
    bound = Fraction(r["printed"]) * r["scale"]
    holds = v <= bound
    print(f"  printed bound {r['printed']} * 2^-14: {'implied' if holds else 'NOT implied'}, margin "
          f"{float((bound - v) / r['scale']):.3e} * 2^-14")
    ok = same and holds
    if controls:
        print("Red controls:")
        last = Fraction(10) ** -len(r["printed"].split(".")[1])
        lower = (Fraction(r["printed"]) - last) * r["scale"]
        c1 = not v <= lower
        print(f"  C1 a bound one last-digit unit lower ({float(lower / r['scale']):.10f} * 2^-14): "
              f"{'not implied, as expected' if c1 else 'IMPLIED (control failed)'}")
        bad = list(r["nw"])
        bad[1] = r["D"] + 1
        try:
            classes(bad, r["D"])
            c2 = False
        except ValueError as e:
            c2 = True
            print(f"  C2 n_1 = D + 1 (a density above 1): refused ({e})")
        if not c2:
            print("  C2 n_1 = D + 1: ACCEPTED (control failed)")
        bump = list(r["nw"])
        bump[1] += 1
        v2 = exact(r["t"], r["n"], r["D"], bump, workers, log=lambda s: None)
        c3 = v2 != r["claimed"] and v2 != v
        print(f"  C3 n_1 raised by one unit: density {float(v2 / r['scale']):.15f} * 2^-14, "
              f"{'differs from the stated fraction, as expected' if c3 else 'SAME (control failed)'}")
        ok &= c1 and c2 and c3
    print(f"RESULT: {key} {'PASS' if ok else 'FAIL'}")
    print(f"run time {time.time() - t0:.0f} s")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    workers = next((int(a.split("=")[1]) for a in sys.argv[1:] if a.startswith("--workers=")),
                   min(8, os.cpu_count() or 1))
    if not args:
        print(__doc__)
        sys.exit(2)
    build()
    if args[0] == "selftest":
        ok = selftest(workers)
    else:
        root = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
        ok = check_row(args[0], root, "--controls" in sys.argv, workers)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
