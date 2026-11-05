#!/usr/bin/env python3
"""Exact replay of PR #307's upper bound on C_89b (Ramsey multiplicity constant of K_5): a step graphon on 256 equal
blocks indexed by F_2^8 with red densities 0, p = 2648/4041 and 1.

  1. The frozen full matrix (certificates/c5_cayley256_fine.json.gz) is read and checked: equal block weights, every
     entry 0, the one numerator of p, or the denominator; symmetric; R[x][y] = R[0][x xor y] for every pair (so the
     translations of F_2^8, which act transitively on the blocks, preserve it); and R[0][z] as the two-mask definition
     file says.
  2. c89b_count.c counts the 5-tuples of blocks behind each red monomial p^k and blue monomial (1 - p)^k twice: with
     the first vertex in block 0 (root), and over every ordered tuple with no symmetry assumed (all); the second must
     be exactly 256 times the first.
  3. The density sum_k [red_k p^k + blue_k (1 - p)^k] / 256^4 is formed in exact rationals, as a polynomial in p and
     at p = 2648/4041 (and at 19/29 for the simpler alternative), and compared with the stated values.

Usage: python3 -I c89b_check.py selftest | pr307 [CERTS] [--controls] [--workers=N]
Python 3.9+, standard library only, and a C compiler (cc).
"""
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from fractions import Fraction
from itertools import product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
EXE = os.path.join(BIN, "c89b_count")

PR307 = dict(
    matrix=("pr307/c5_cayley256_fine.json.gz", "78978a9f0cfc8dcf56ca371994b155f07768964fad62fe3e037d457b9aef0125"),
    definition=("pr307/c5_cayley256_definition.json",
                "52ba8e51bba96da34526b3a777ad545c5a37ca3a10937c5c84d03273f5ea0a85"),
    t=5,
    claimed=Fraction(15058298592194934908767847, 9131939205947984819238469632),
    printed="0.001648970526",
    poly=([8734, -5530, 4320, -410, 395, -60, 10], 4194304),
    simple=(Fraction(19, 29), Fraction(4113966854539, 2494869834563584)),
)


def build():
    os.makedirs(BIN, exist_ok=True)
    src = os.path.join(HERE, "c89b_count.c")
    if not os.path.exists(EXE) or os.path.getmtime(EXE) < os.path.getmtime(src):
        subprocess.run(["cc", "-O2", "-std=c99", "-o", EXE, src], check=True)


def codes_of(R, D):
    """0 / 1 / 2 codes and p; refuses entries outside [0, D] or more than one fractional value."""
    fr = sorted({x for row in R for x in row} - {0, D})
    if any(x < 0 or x > D for row in R for x in row):
        raise ValueError("an entry is outside [0, D]")
    if len(fr) > 1:
        raise ValueError(f"{len(fr)} different fractional densities (this counter takes one)")
    p = Fraction(fr[0], D) if fr else Fraction(0)
    return [[0 if x == 0 else 2 if x == D else 1 for x in row] for row in R], p


def count(t, codes, mode, workers=1):
    text = f"{len(codes)}\n" + "\n".join(" ".join(map(str, row)) for row in codes) + "\n"
    if mode == "root":
        runs = [[EXE, str(t), "root"]]
    else:
        runs = [[EXE, str(t), "all", str(i), str(workers)] for i in range(workers)]
    procs = [subprocess.Popen(r, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for r in runs]
    tab = {}
    for pr in procs:
        out, err = pr.communicate(text)
        if pr.returncode:
            raise ValueError(err.strip())
        for line in out.split("\n"):
            if line.strip():
                c, k, n = line.split()
                tab[(c, int(k))] = tab.get((c, int(k)), 0) + int(n)
    return tab


def density(tab, p, m, t, root=True):
    s = sum(n * (p if c == "R" else 1 - p) ** k for (c, k), n in tab.items())
    return Fraction(s) / Fraction(m) ** (t - 1 if root else t)


def poly(tab, m, t):
    """Coefficients of the density as a polynomial in p (exact rationals)."""
    deg = max(k for _, k in tab)
    co = [Fraction(0)] * (deg + 1)
    for (c, k), n in tab.items():
        if c == "R":
            co[k] += n
        else:
            for j in range(k + 1):
                co[j] += n * comb(k, j) * (-1) ** j
    return [x / Fraction(m) ** (t - 1) for x in co]


def cayley_ok(R):
    m = len(R)
    return all(R[x][y] == R[0][x ^ y] for x in range(m) for y in range(m))


def brute(R, D, t):
    m = len(R)
    s = Fraction(0)
    for xs in product(range(m), repeat=t):
        r = b = Fraction(1)
        for i in range(t):
            for j in range(i + 1, t):
                q = Fraction(R[xs[i]][xs[j]], D)
                r *= q
                b *= 1 - q
        s += r + b
    return s / Fraction(m) ** t


def selftest(workers):
    ok = True
    rng = random.Random(307)
    print("1. Random small matrices, every ordered tuple enumerated in Python, against the C counts:")
    for trial in range(6):
        n = 3 if trial < 4 else 2
        m, t, D = 2 ** n, (5 if trial % 2 else 4), 9
        num = rng.randint(1, D - 1)
        q = [rng.choice([0, num, D]) for _ in range(m)]
        R = [[q[x ^ y] for y in range(m)] for x in range(m)]   # a Cayley matrix on F_2^n
        codes, p = codes_of(R, D)
        a, b = count(t, codes, "root"), count(t, codes, "all", workers=3)
        same = all(b.get(k, 0) == m * a.get(k, 0) for k in set(a) | set(b))
        v, w = density(a, p, m, t), brute(R, D, t)
        good = same and v == w
        ok &= good
        print(f"  Cayley m={m} t={t} q={q}/{D}: all = m x root {same}, brute force {'agrees' if v == w else 'DIFFERS'}")
    for trial in range(3):
        m, t, D = 5, 4, 7
        num = rng.randint(1, D - 1)
        R = [[0] * m for _ in range(m)]
        for i in range(m):
            for j in range(i, m):
                R[i][j] = R[j][i] = rng.choice([0, num, D])
        codes, p = codes_of(R, D)
        v, w = density(count(t, codes, "all", workers=2), p, m, t, root=False), brute(R, D, t)
        ok &= v == w
        print(f"  arbitrary symmetric m={m} t={t}: all-tuple count against brute force "
              f"{'agrees' if v == w else 'DIFFERS'}")
    print("2. f = 1/2 everywhere (m = 4, one block pair class) gives 2^(1 - C(t,2)):")
    for t in (4, 5):
        R = [[1] * 4 for _ in range(4)]
        codes, p = codes_of(R, 2)
        v = density(count(t, codes, "root"), p, 4, t)
        good = v == Fraction(2, 2 ** (t * (t - 1) // 2))
        ok &= good
        print(f"  t={t}: {v} {'as expected' if good else 'WRONG'}")
    print("RESULT: selftest", "PASS" if ok else "FAIL")
    return ok


def check_pr307(root, controls, workers):
    c = PR307
    t0 = time.time()
    raw = open(os.path.join(root, c["matrix"][0]), "rb").read()
    h = hashlib.sha256(raw).hexdigest()
    print(f"Row pr307 (PR #307, rpinho, round 4): {c['matrix'][0]} sha256 {h} "
          f"({'as pinned' if h == c['matrix'][1] else 'NOT AS PINNED'})")
    assert h == c["matrix"][1]
    d = json.loads(gzip.decompress(raw))
    D, Rm, wts = d["edge_probability_denominator"], d["red_probability_numerators"], d["block_weights"]
    m = len(Rm)
    assert m == 256 and all(len(row) == m for row in Rm) and len(wts) == m, "not a 256 x 256 matrix"
    assert len(set(wts)) == 1, "block weights are not equal"
    codes, p = codes_of(Rm, D)
    print(f"  {m} blocks of equal weight; denominator {D}; entries 0, {p.numerator * D // p.denominator}, {D}: "
          f"p = {p}")
    sym = all(Rm[i][j] == Rm[j][i] for i in range(m) for j in range(m))
    cay = cayley_ok(Rm)
    print(f"  symmetric: {sym}; R[x][y] = R[0][x xor y] for all 65,536 pairs: {cay}")
    raw2 = open(os.path.join(root, c["definition"][0]), "rb").read()
    h2 = hashlib.sha256(raw2).hexdigest()
    df = json.loads(raw2)
    H, Fm = int(df["hard_red_mask_hex"], 16), int(df["fractional_mask_hex"], 16)
    pd = Fraction(df["probability_numerator"], df["probability_denominator"])
    want = [D if H >> z & 1 else (p.numerator * D // p.denominator if Fm >> z & 1 else 0) for z in range(m)]
    defn = h2 == c["definition"][1] and H & Fm == 0 and pd == p and Rm[0] == want
    print(f"  two-mask definition ({c['definition'][0]}, sha256 {h2[:8]}...): masks disjoint, p and row 0 match: "
          f"{defn} ({bin(H).count('1')} red, {bin(Fm).count('1')} fractional, {m - bin(H | Fm).count('1')} zero)")
    ta = time.time()
    a = count(c["t"], codes, "root")
    tb = time.time()
    b = count(c["t"], codes, "all", workers)
    tc = time.time()
    same = all(b.get(k, 0) == m * a.get(k, 0) for k in set(a) | set(b))
    print(f"  counts, root: {sorted(a.items())} ({tb - ta:.0f} s)")
    print(f"  every ordered tuple, no symmetry used: exactly {m} x root in every cell: {same} ({tc - tb:.0f} s)")
    v = density(a, p, m, c["t"])
    print(f"  exact density = {v.numerator}/{v.denominator} = {float(v):.16f}")
    same_f = v == c["claimed"]
    print(f"  stated fraction: {'IDENTICAL' if same_f else 'DIFFERENT'}")
    bound = Fraction(c["printed"])
    holds = v < bound
    print(f"  printed bound {c['printed']}: {'implied' if holds else 'NOT implied'} (margin {float(bound - v):.3e})")
    co, den = c["poly"]
    pc = poly(a, m, c["t"])
    same_p = [x * den for x in pc] == [Fraction(x) for x in co] + [Fraction(0)] * (len(pc) - len(co))
    print(f"  density as a polynomial in p, times {den}: {[str(x * den) for x in pc]}: "
          f"{'the stated polynomial' if same_p else 'NOT the stated polynomial'}")
    ps, fs = c["simple"]
    same_s = density(a, ps, m, c["t"]) == fs
    print(f"  at p = {ps} (the simpler alternative): {'the stated ' + str(fs) if same_s else 'DIFFERENT'}")
    ok = sym and cay and defn and same and same_f and holds and same_p and same_s
    if controls:
        print("Red controls:")
        last = Fraction(1, 10 ** len(c["printed"].split(".")[1]))
        c1 = not v < bound - last
        print(f"  C1 a bound one last-digit unit lower ({float(bound - last):.12f}): "
              f"{'not implied, as expected' if c1 else 'IMPLIED (control failed)'}")
        bad = [row[:] for row in Rm]
        bad[1][2] = (bad[1][2] + 2648) % (D + 1) if bad[1][2] != 2648 else 0
        c2 = not cayley_ok(bad)
        print(f"  C2 one entry R[1][2] changed: translation invariance {'fails, as expected' if c2 else 'HOLDS (control failed)'}")
        try:
            codes_of([row[:] for row in Rm[:1]] + [[1] + Rm[1][1:]] + Rm[2:], D)
            c3 = False
        except ValueError as e:
            c3 = True
            print(f"  C3 a second fractional value (R[1][0] = 1/{D}): refused ({e})")
        if not c3:
            print("  C3 a second fractional value: ACCEPTED (control failed)")
        asym = [row[:] for row in codes]
        asym[1][2] = (asym[1][2] + 1) % 3
        try:
            count(c["t"], asym, "root")
            c4 = False
        except ValueError as e:
            c4 = True
            print(f"  C4 an asymmetric matrix: refused by the counter ({e})")
        if not c4:
            print("  C4 an asymmetric matrix: ACCEPTED (control failed)")
        ok &= c1 and c2 and c3 and c4
    print(f"RESULT: pr307 {'PASS' if ok else 'FAIL'}")
    print(f"run time {time.time() - t0:.0f} s")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    workers = next((int(a.split("=")[1]) for a in sys.argv[1:] if a.startswith("--workers=")),
                   min(8, os.cpu_count() or 1))
    if not args or args[0] not in ("selftest", "pr307"):
        print(__doc__)
        sys.exit(2)
    build()
    if args[0] == "selftest":
        ok = selftest(workers)
    else:
        ok = check_pr307(args[1] if len(args) > 1 else os.path.join(HERE, "certs"), "--controls" in sys.argv,
                         workers)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
