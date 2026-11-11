#!/usr/bin/env python3
"""Exact replay of two 28a upper-bound rows (teorth/optimizationproblems, main 87b4ebb), on a G_2(4) graph that we
build from the construction of Jenrich's paper (a second, independent route to the graph of borsuk_check.py, which
uses Brouwer's model with a diagonal Hermitian form):

  [Bon2014]  C_28 <= 65: "Two-distance counterexample (416 points on S^64 in R^65); cannot be partitioned into 83
             smaller-diameter sets (so needs >= 84)" (Bondarenko, arXiv 1305.2584v3, Theorem 1).
  [JB2014]   C_28 <= 64: "A 352-point two-distance subset giving a counterexample in R^64; cannot be partitioned into
             70 smaller-diameter sets (so needs >= 71)" (Jenrich, arXiv 1308.0206v6; Jenrich and Brouwer).

Construction (Jenrich, section 4, as text): F_16 = F_2[a]/(a^4 + a + 1); PG(2, 16) with the Hermitian form
h(x, y) = x1 conj(y3) + x2 conj(y2) + x3 conj(y1), conj(z) = z^4; vertices: the orthogonal bases {a, b, c} of
non-isotropic points; the triangle of a basis: the isotropic points on its three lines (the line ab is the polar of c);
two vertices adjacent when their triangles share exactly 3 isotropic points.

Computed exactly (integers, bit sets, rationals):
  1. 273 points, 65 isotropic, 416 bases, 15 isotropic points per triangle; strongly regular (416, 100, 36, 20) on
     every pair; Bondarenko's formulas then give f = 65, beta = 26, p = 1/5, q = -1/15 (computed from his equations);
  2. clique number 5: a 5-clique, and no 6-clique by an exhaustive bit-set search;
  3. [Bon2014]: with y_i the columns of A + 4I and z_i = y_i - mean, <z_i, z_j> = (A^2)_ij + 8 A_ij + 16 [i = j] - 26
     = 90, 18, -6 (equal, adjacent, not), i.e. x_i = z_i/sqrt 90 has inner products 1, 1/5, -1/15 as stated; the Gram
     matrix G = 96 I + 24 A - 6 J satisfies G^2 = 576 G (from the verified identities A^2 = 80 I + 16 A + 20 J and
     A J = 100 J), so G/576 is an idempotent of trace 416 * 90/576 = 65: the 416 points span exactly R^65 and lie on a
     sphere; also rank 65 modulo 2^61 - 1 (second route). Squared distances 144 (adjacent) and 192 (not): the diameter
     is attained exactly on non-adjacent pairs, so a part of smaller diameter is a clique, at most 5 points:
     at least ceil(416/5) = 84 > 66 parts;
  4. [JB2014]: for an isotropic q0 (all 65 are checked), B = vertices whose triangle contains q0, C = the rest; |B| = 96
     in three components B_1, B_2, B_3 of 32; claim (1): |N(i) n B_h| = 20 (i in B_h), 0 (i in B - B_h), 8 (i in C)
     for every vertex i and h; p = 1_{B_2} - 1_{B_3} satisfies A p = 20 p (so p lies in the column space of A + 4I),
     <p, y_i> = 0 on C u B_1 and +-24 on B_2, B_3. So the 352 points y_i (i in C u B_1) span a linear space of
     dimension at most 65 inside the hyperplane sum = 104, i.e. an affine space of dimension at most 64; the Gram
     matrix of their differences has rank 64 modulo 2^61 - 1, so exactly 64. Distances are again 144 and 192 with the
     diameter on non-adjacent pairs; parts of smaller diameter are cliques: at least ceil(352/5) = 71 > 65 parts.
  The statements quoted above are checked to be in the two LaTeX sources (pinned by ./fetch_bonjb.sh).
Red controls (--controls): the adjacency rule |T n T'| = 2 is not SRG(416, 100, 36, 20); a non-isotropic q0 does not
give the 96/320 split; p replaced by 1_{B_2} + 1_{B_3} is not an eigenvector for 20 and is not orthogonal to the y_i
on C; with an added edge completing a 6-clique the part count drops to ceil(416/6) = 70.
usage: python3 -I bon_jb_check.py CERTS_ROOT [--controls]
"""
import itertools
import os
import sys
from fractions import Fraction

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv
PRIME = 2 ** 61 - 1


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ---------------------------------------------------------------- F_16
def gf_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 0b10011            # a^4 = a + 1
    return r


MUL = [[gf_mul(a, b) for b in range(16)] for a in range(16)]
INV = [0] + [next(b for b in range(1, 16) if MUL[a][b] == 1) for a in range(1, 16)]


def conj(z):
    return MUL[MUL[z][z]][MUL[z][z]]      # z^4


def h(x, y):
    return MUL[x[0]][conj(y[2])] ^ MUL[x[1]][conj(y[1])] ^ MUL[x[2]][conj(y[0])]


def points():
    pts = []
    for v in itertools.product(range(16), repeat=3):
        if v == (0, 0, 0):
            continue
        lead = next(c for c in v if c)
        if lead == 1:
            pts.append(v)
    return pts


def build(adj_size=3):
    P = points()
    iso = [p for p in P if h(p, p) == 0]
    non = [p for p in P if h(p, p) != 0]
    iso_idx = {p: i for i, p in enumerate(iso)}
    bases = set()
    for a in non:
        perp = [b for b in non if b != a and h(a, b) == 0]
        for b, c in itertools.combinations(perp, 2):
            if h(b, c) == 0:
                bases.add(frozenset((a, b, c)))
    bases = sorted(bases, key=lambda s: sorted(s))
    tri = []
    for B in bases:
        pts = set()
        for c in B:
            pts |= {iso_idx[q] for q in iso if h(q, c) == 0}       # the line through the other two is the polar of c
        tri.append(frozenset(pts))
    n = len(bases)
    nb = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if len(tri[i] & tri[j]) == adj_size:
                nb[i] |= 1 << j
                nb[j] |= 1 << i
    return P, iso, non, bases, tri, nb


def popcount(x):
    return bin(x).count("1")


def srg(nb):
    n = len(nb)
    ks = {popcount(x) for x in nb}
    if len(ks) != 1:
        return None
    k = ks.pop()
    lam, mu = set(), set()
    for i in range(n):
        for j in range(i + 1, n):
            c = popcount(nb[i] & nb[j])
            (lam if (nb[i] >> j) & 1 else mu).add(c)
    if len(lam) != 1 or len(mu) != 1:
        return None
    return (n, k, lam.pop(), mu.pop())


def max_clique(nb, cap=7):
    n = len(nb)
    best = [0, []]

    def rec(cl, cand):
        if len(cl) > best[0]:
            best[0], best[1] = len(cl), list(cl)
            if best[0] >= cap:
                return
        while cand:
            if len(cl) + popcount(cand) <= best[0]:
                return
            v = cand.bit_length() - 1
            cand &= ~(1 << v)
            cl.append(v)
            rec(cl, cand & nb[v])
            cl.pop()

    rec([], (1 << n) - 1)
    return best[0], best[1]


def rank_mod(rows):
    rows = [list(r) for r in rows]
    rank, ncol = 0, len(rows[0])
    for c in range(ncol):
        piv = next((r for r in range(rank, len(rows)) if rows[r][c] % PRIME), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][c], PRIME - 2, PRIME)
        pr = [(x * inv) % PRIME for x in rows[rank]]
        rows[rank] = pr
        for r in range(len(rows)):
            if r != rank and rows[r][c] % PRIME:
                f = rows[r][c]
                rows[r] = [(x - f * y) % PRIME for x, y in zip(rows[r], pr)]
        rank += 1
    return rank


def comps(nb, S):
    S = set(S)
    out = []
    while S:
        v = S.pop()
        comp, stack = {v}, [v]
        while stack:
            u = stack.pop()
            for w in list(S):
                if (nb[u] >> w) & 1:
                    S.discard(w)
                    comp.add(w)
                    stack.append(w)
        out.append(sorted(comp))
    return out


def check_sources(root):
    bon = open(os.path.join(root, "bon2014", "borsuk.tex"), encoding="utf-8").read()
    jb = open(os.path.join(root, "jb2014", "BORSUK64.tex"), encoding="latin-1").read()
    for s in ("$416$ points on the unit sphere $S^{64}\\subset\\mathbb{R}^{65}$ which cannot be partitioned into $83$ "
              "parts of smaller diameter", "$p=1/5$ and $q=-1/15$", "Hence $b(65)\\ge b_2(65)\\ge 84$"):
        need(s in bon, f"not in Bondarenko's source: {s}")
    for s in ("a_1 \\overline {b_3} + a_2 \\overline {b_2} + a_3 \\overline {b_1}", "2\\times{(100-36-1+(4-1)^2)} = 144",
              "2\\times{(100-20+4^2)} = 192", "contains 352 vectors", "a division into less than 71",
              "20 & \\mbox{if $i \\in B_h$}"):
        need(s in jb, f"not in Jenrich's source: {s}")


def main():
    ok = True
    check_sources(ROOT)
    print("the quoted statements are in both LaTeX sources: True")
    P, iso, non, bases, tri, nb = build()
    n = len(bases)
    par = srg(nb)
    sizes = {len(t) for t in tri}
    print(f"PG(2,16) with h(x,y) = x1 y3^4 + x2 y2^4 + x3 y1^4: {len(P)} points, {len(iso)} isotropic, {len(non)} "
          f"non-isotropic; {n} orthogonal bases; triangle sizes {sorted(sizes)}; strongly regular {par}")
    good_graph = (len(P), len(iso), n, sizes, par) == (273, 65, 416, {15}, (416, 100, 36, 20))
    v, k, lam, mu = 416, 100, 36, 20
    s = Fraction(lam - mu, 2) - Fraction(isqrt_exact((lam - mu) ** 2 + 4 * (k - mu)), 2)
    disc = isqrt_exact((lam - mu) ** 2 + 4 * (k - mu))
    f = Fraction(v - 1, 2) - Fraction(2 * k + (v - 1) * (lam - mu), 2 * disc)
    beta = Fraction(s * s + k + k * (lam - 2 * s) + (v - k - 1) * mu, v)
    p = (lam - 2 * s - beta) / (s * s + k - beta)
    q = (mu - beta) / (s * s + k - beta)
    print(f"  Bondarenko's formulas: s = {s}, f = {f}, beta = {beta}, p = {p}, q = {q}")
    good_graph &= (s, f, beta, p, q) == (-4, 65, 26, Fraction(1, 5), Fraction(-1, 15))
    omega, cl = max_clique(nb)
    print(f"  clique number {omega} (a clique {cl}; no larger by exhaustive search)")
    good_graph &= omega == 5
    # [Bon2014]
    G = {"I": 96, "A": 24, "J": -6}
    # products in the algebra spanned by I, A, J: A^2 = 80 I + 16 A + 20 J (k - mu = 80, lam - mu = 16, mu = 20),
    # AJ = JA = 100 J, J^2 = 416 J
    def mul_alg(X, Y):
        out = {"I": 0, "A": 0, "J": 0}
        for a, ca in X.items():
            for b, cb in Y.items():
                c = ca * cb
                if a == "I":
                    out[b] += c
                elif b == "I":
                    out[a] += c
                elif a == "A" and b == "A":
                    out["I"] += 80 * c
                    out["A"] += 16 * c
                    out["J"] += 20 * c
                elif {a, b} == {"A", "J"}:
                    out["J"] += 100 * c
                else:
                    out["J"] += 416 * c
        return out
    G2 = mul_alg(G, G)
    idem = all(G2[t] == 576 * G[t] for t in G)
    trace_rank = Fraction(416 * (96 - 6), 576)
    gram_ok = True
    for i in range(n):
        for j in range(i, n):
            c = popcount(nb[i] & nb[j]) if i != j else k
            a = (nb[i] >> j) & 1
            z = c + 8 * a + 16 * (i == j) - 26
            want = 90 if i == j else (18 if a else -6)
            gram_ok &= z == want
    rows = [[(96 if i == j else 0) + 24 * ((nb[i] >> j) & 1) - 6 for j in range(n)] for i in range(n)]
    r_mod = rank_mod(rows)
    d_adj, d_non = 2 * 90 - 2 * 18, 2 * 90 + 2 * 6
    parts = -(-416 // omega)
    good_bon = good_graph and idem and trace_rank == 65 and gram_ok and r_mod == 65 and d_non > d_adj and parts == 84
    print(f"[Bon2014] <z_i, z_j> = (A^2)_ij + 8 A_ij + 16 [i = j] - 26 is 90 / 18 / -6 on every pair: {gram_ok}; "
          f"G^2 = 576 G in the algebra of I, A, J: {idem}; rank = trace(G)/576 = {trace_rank}; rank mod 2^61 - 1: "
          f"{r_mod}")
    print(f"  squared distances (times 1/90): adjacent {d_adj}, non-adjacent {d_non}; diameter on non-adjacent pairs; "
          f"parts of smaller diameter are cliques of size <= {omega}: at least ceil(416/{omega}) = {parts} > 66")
    print(f"RESULT [Bon2014]: 416 points of S^64 in R^65 need at least {parts} > 66 parts of smaller diameter, so "
          f"C_28 <= 65: {'HOLDS' if good_bon else 'FAILS'}")
    ok &= good_bon
    # [JB2014]
    all_q0 = True
    first = None
    for q0 in range(len(iso)):
        B = [i for i in range(n) if q0 in tri[i]]
        C = [i for i in range(n) if q0 not in tri[i]]
        cs = sorted(comps(nb, B), key=lambda c: c[0])
        good = len(B) == 96 and len(C) == 320 and sorted(len(c) for c in cs) == [32, 32, 32]
        if good:
            masks = [sum(1 << x for x in c) for c in cs]
            for i in range(n):
                for hh, m in enumerate(masks):
                    cnt = popcount(nb[i] & m)
                    want = 20 if (m >> i) & 1 else (0 if any((mm >> i) & 1 for mm in masks) else 8)
                    good &= cnt == want
        all_q0 &= good
        if first is None:
            first = (B, C, cs)
    B, C, cs = first
    B1, B2, B3 = cs
    pv = [0] * n
    for x in B2:
        pv[x] = 1
    for x in B3:
        pv[x] = -1
    Ap = [sum(pv[j] for j in range(n) if (nb[i] >> j) & 1) for i in range(n)]
    eig = Ap == [20 * x for x in pv]
    ydot = [sum(pv[j] * (4 if i == j else ((nb[i] >> j) & 1)) for j in range(n)) for i in range(n)]
    sub = sorted(C + B1)
    orth = all(ydot[i] == 0 for i in sub) and all(ydot[i] == 24 for i in B2) and all(ydot[i] == -24 for i in B3)
    i0 = sub[0]

    def gy(i, j):
        c = popcount(nb[i] & nb[j]) if i != j else k
        return c + 8 * ((nb[i] >> j) & 1) + 16 * (i == j)

    drows = [[gy(i, j) - gy(i, i0) - gy(i0, j) + gy(i0, i0) for j in sub[1:]] for i in sub[1:]]
    r_aff = rank_mod(drows)
    parts_jb = -(-len(sub) // 5)
    good_jb = (good_graph and all_q0 and eig and orth and len(sub) == 352 and r_aff == 64 and parts_jb == 71)
    print(f"[JB2014] for all 65 isotropic q0: |B| = 96 in three components of 32, |C| = 320, claim (1) (20 / 0 / 8) for "
          f"every vertex and component: {all_q0}")
    print(f"  q0 = the first isotropic point: p = 1_B2 - 1_B3 has A p = 20 p: {eig}; <p, y_i> = 0 on C u B_1 and +-24 on "
          f"B_2, B_3: {orth}; |C u B_1| = {len(sub)}")
    print(f"  so dim span <= 65 inside the hyperplane sum = 104: affine dimension <= 64; Gram of differences rank mod "
          f"2^61 - 1: {r_aff} (so exactly 64); parts of smaller diameter are cliques (<= 5 points): at least "
          f"ceil(352/5) = {parts_jb} > 65")
    print(f"RESULT [JB2014]: 352 points in R^64 need at least {parts_jb} > 65 parts of smaller diameter, so C_28 <= 64: "
          f"{'HOLDS' if good_jb else 'FAILS'}")
    ok &= good_jb
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        _, _, _, _, _, nb2 = build(adj_size=2)
        par2 = srg(nb2)
        print(f"  adjacency at |T n T'| = 2: {par2 if par2 else 'not strongly regular'} (degree "
              f"{sorted({popcount(x) for x in nb2})})")
        res.append(("|T n T'| = 2", par2 != (416, 100, 36, 20)))
        # a non-isotropic q0: the bases whose triangle "contains" it cannot be formed (triangles hold isotropic points)
        nonq = non[0]
        Bn = [i for i in range(n) if any(h(nonq, c) == 0 for c in bases[i])]
        print(f"  non-isotropic q0 (the bases with a point orthogonal to it, the description of B for isotropic q0): "
              f"{len(Bn)} vertices, not 96")
        res.append(("non-isotropic q0", len(Bn) != 96))
        pv2 = [1 if (x in B2 or x in B3) else 0 for x in range(n)]
        Ap2 = [sum(pv2[j] for j in range(n) if (nb[i] >> j) & 1) for i in range(n)]
        yd2 = [sum(pv2[j] * (4 if i == j else ((nb[i] >> j) & 1)) for j in range(n)) for i in C[:5]]
        res.append(("1_B2 + 1_B3", Ap2 != [20 * x for x in pv2] and any(yd2)))
        print(f"  p' = 1_B2 + 1_B3: A p' = 20 p': {Ap2 == [20 * x for x in pv2]}; <p', y_i> on C: {yd2[:3]}")
        nb3 = list(nb)
        out = next(x for x in range(n) if all((nb[c] >> x) & 1 for c in cl[:4]) and x not in cl and not (nb[cl[4]] >> x) & 1)
        nb3[out] |= 1 << cl[4]
        nb3[cl[4]] |= 1 << out
        om3, _ = max_clique(nb3)
        print(f"  an added edge completing a 6-clique: clique number {om3}, parts ceil(416/{om3}) = {-(-416 // om3)}")
        res.append(("added edge", om3 == 6 and -(-416 // om3) == 70))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


def isqrt_exact(x):
    import math
    r = math.isqrt(x)
    need(r * r == x, "not a perfect square")
    return r


if __name__ == "__main__":
    main()
