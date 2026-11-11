#!/usr/bin/env python3
"""Exact replay of the 89a upper-bound row [PPSS2024] (teorth/optimizationproblems, main 87b4ebb; added in PR #255):
C_89 <= 4551721 * 2^-24 * 3^-2 < 0.030145, "the sequence of blow-ups of a 2-coloured Cayley graph on 768 vertices"
(Parczyk, Pokutta, Spiegel and Szabo, Theorem 1.1).  The graph is read from the authors' data (graph6, Zenodo record
6602512); none of their programs was read or run.

C_89 is the limit of the least density of monochromatic K_4 in 2-colourings of K_n.  For the blow-ups C[m] of a graph
C on n vertices (a vertex replaced by an independent set of m vertices), the density of K_4 plus independent 4-sets
tends to (hom(K_4, C) + hom(K_4, C-bar with a loop at every vertex)) / n^4 (the paper's blow-up lemma, after Thomason;
read, not checked), so the row's bound is that number.  k4hom.c (our own program) computes it twice:
  route A: hom(K_4, X) = sum over roots a of tr((X D_a)^3), D_a = diag(X_a.), for X = C and X = looped complement;
  route B: clique counts c_1, ..., c_4 of C and of its complement, hom = 24 c_4 for C (no loops) and
           c_1 + 14 c_2 + 36 c_3 + 24 c_4 for the looped complement (surjections of 4 points onto k points);
the routes must agree; the per-root traces are reported (equal for every root, as for a vertex-transitive graph).
Integers only (64-bit in C; the total is below 768^4 < 2^39); the density is compared with 4551721/(2^24 * 9) exactly.
Red controls (--controls): one edge toggled changes the total (both routes still agree); a G(768, 1/2) random graph
(seeded) gives about 1/32 > 0.030145; counting the complement without loops gives a different number.
usage: python3 -I ppss_check.py CERTS_ROOT BIN_DIR [--controls]
"""
import os
import random
import subprocess
import sys
from fractions import Fraction

ROOT, BIN = sys.argv[1], sys.argv[2]
CONTROLS = "--controls" in sys.argv


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def run(path, *extra):
    out = subprocess.run([os.path.join(BIN, "k4hom"), path, *extra], capture_output=True, text=True, check=True).stdout
    d = {}
    for line in out.splitlines():
        t = line.split()
        d[t[0]] = t[1:]
    n, edges = int(d["n"][0]), int(d["n"][2])
    Ar, Ab = int(d["routeA"][1]), int(d["routeA"][3])
    br = [int(x) for x in d["routeB"][1:5]]
    bb = [int(x) for x in d["routeB"][6:10]]
    roots = [int(x) for x in d["roots"][1:3] + d["roots"][4:6]]
    return n, edges, Ar, Ab, br, bb, roots


def to_graph6(n, adj):
    bits = [1 if adj[i][j] else 0 for j in range(1, n) for i in range(j)]
    bits += [0] * (-len(bits) % 6)
    head = bytes([126, 63 + (n >> 12 & 63), 63 + (n >> 6 & 63), 63 + (n & 63)]) if n > 62 else bytes([63 + n])
    body = bytes(63 + int("".join(map(str, bits[k:k + 6])), 2) for k in range(0, len(bits), 6))
    return head + body


def main():
    ok = True
    tex = open(os.path.join(ROOT, "ppss", "main.tex"), encoding="utf-8").read()
    need("c_4 \\leq 4551721 \\cdot 2^{-24} \\cdot 3^{-2} < 0.030145" in tex, "the theorem's statement was not found")
    g6 = os.path.join(ROOT, "ppss", "c4.graph6.txt")
    n, edges, Ar, Ab, br, bb, roots = run(g6)
    hom_r_B = 24 * br[3]
    hom_b_B = bb[0] + 14 * bb[1] + 36 * bb[2] + 24 * bb[3]
    total = Ar + Ab
    dens = Fraction(total, n ** 4)
    stated = Fraction(4551721, 2 ** 24 * 9)
    good = (n == 768 and Ar == hom_r_B and Ab == hom_b_B and roots[0] == roots[1] and roots[2] == roots[3]
            and dens == stated and dens < Fraction(30145, 10 ** 6))
    print(f"the graph (graph6, {n} vertices, {edges} edges, degree {2 * edges // n})")
    print(f"  route A: hom(K_4, C) = {Ar}, hom(K_4, looped complement) = {Ab}")
    print(f"  route B: cliques of C (c_1..c_4) = {br}, 24 c_4 = {hom_r_B}; of the complement = {bb}, "
          f"c_1 + 14 c_2 + 36 c_3 + 24 c_4 = {hom_b_B}; routes agree: {Ar == hom_r_B and Ab == hom_b_B}")
    print(f"  per-root traces: C {roots[0]}..{roots[1]}, complement {roots[2]}..{roots[3]} (one value each)")
    print(f"  density = {total}/768^4 = {dens} = {float(dens):.15f}; 4551721 * 2^-24 * 3^-2 = {stated}: equal "
          f"{dens == stated}; < 0.030145: {dens < Fraction(30145, 10 ** 6)}")
    print(f"RESULT [PPSS2024]: the blow-ups of the 768-vertex graph have monochromatic K_4 density tending to exactly "
          f"4551721 * 2^-24 * 3^-2 = {float(stated):.12f} < 0.030145, so C_89 <= that: {'HOLDS' if good else 'FAILS'}")
    ok &= good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        n2, e2, Ar2, Ab2, br2, bb2, _ = run(g6, "flip", "0", "1")
        t2 = Ar2 + Ab2
        agree = Ar2 == 24 * br2[3] and Ab2 == bb2[0] + 14 * bb2[1] + 36 * bb2[2] + 24 * bb2[3]
        print(f"  edge 0-1 toggled: total {t2} (was {total}), routes agree: {agree}")
        res.append(("one edge toggled", t2 != total and agree))
        rng = random.Random(89)
        adj = [[False] * n for _ in range(n)]
        for j in range(1, n):
            for i in range(j):
                adj[i][j] = adj[j][i] = rng.random() < 0.5
        rp = os.path.join(BIN, "random768.g6")
        open(rp, "wb").write(to_graph6(n, adj))
        _, _, Ar3, Ab3, _, _, _ = run(rp)
        d3 = Fraction(Ar3 + Ab3, n ** 4)
        print(f"  G(768, 1/2), seed 89: density {float(d3):.6f}")
        res.append(("a random graph", d3 > Fraction(30145, 10 ** 6)))
        unl = Ar + 24 * bb[3]
        print(f"  complement without loops: {unl} != {total}")
        res.append(("no loops in the complement", unl != total and Fraction(unl, n ** 4) != stated))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
