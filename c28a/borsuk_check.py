"""Replay of the 28a row [Gri2026], C_28 <= 63 (Grinsztajn, PR #84): Borsuk's conjecture fails in dimension 63.

Written from the proof PDF of maaxgrin/borsuk-63-counterexample @ cdcdbeac (read as text) and the 28a page. Neither
verify_borsuk63.py, verify_borsuk63_sage.sage nor export_borsuk63_certificates.py was read or run, and the DIMACS
files are not used: we build the G_2(4) graph ourselves from Brouwer's projective model as the PDF describes it.

  F_16 = F_2[a]/(a^4 + a + 1); PG(2, 16); h(u, w) = u0 w0^4 + u1 w1^4 + u2 w2^4 (Hermitian). Vertices: unordered
  triples {a1, a2, a3} of pairwise orthogonal non-isotropic points; T(A): the isotropic points on the three lines
  spanned by pairs of A (the line through a_j, a_k is the polar line of a_i); A ~ A' iff |T(A) & T(A')| = 3.

What is computed, exactly (integers, bit sets, rationals, and Q(sqrt 222) with exact sign tests):
  1. 273 points, 65 isotropic; 416 triples; |T(A)| = 15; the graph is strongly regular (416, 100, 36, 20), checked
     on every pair; hence (standard) eigenvalues 100, 20, -4 with multiplicities 1, 65, 350, and G = 96I + 24A - 6J
     is positive semidefinite of rank 65 (we also find rank 65 modulo the prime 2^61 - 1, a lower bound for the rank
     over Q): the Gram matrix of vectors x_v in R^65 with x_u . x_v = 90, 18, -6 (u = v, u ~ v, otherwise).
  2. omega = 5: a 5-clique is found and an exhaustive bit-set search finds no 6-clique.
  3. For an isotropic q0 (we check all 65): B = triples with a point orthogonal to q0, C = the rest; |B| = 96,
     |C| = 320; B induces three components of 32; the degree data of Lemma 1(4).
  4. Lemmas 3-5 from these data: x_c . S_i = 0 (S_i = sum of x_y over B_i), Gram(S_1, S_2, S_3) has rank 2, so the
     x_c lie in the 63-dimensional W = span(S)^perp; z_b = x_b - S_1/32 is orthogonal to S_1, S_2, S_3, |z_b|^2 = 78,
     z_b . x_c = 18 or -6; with t = (sqrt 222 - 1)/13 (78 t^2 + 12 t = 102 exactly) and p = t z_b, the squared
     distances in X = {x_c} u {p} are 144, 192 (pairs in C), 192 - 48t and 192 (to p); so diam^2 = 192, attained.
     The Gram matrix of {x_c} u {z_b} (integers) has rank 63 modulo 2^61 - 1, so dim span X = 63 exactly.
  5. Lemma 6: a subset of smaller diameter is a clique of the graph (with b standing for p), so has at most 5 points;
     ceil(321 / 5) = 65 > 64 = 63 + 1.
  Anchors (the same graph, not the replayed row): the 416 points give Bondarenko's 84 > 66 in R^65; the 352 points
  C u B_3 are orthogonal to S_1 - S_2, so lie in a hyperplane: 71 > 65 in R^64 (the bound of the [JB2014] row).

Controls (each must come out as stated): the adjacency rule |T & T'| = 2 is not SRG(416, 100, 36, 20); a
non-isotropic q0 does not give the 96/320 split with the degree data; t + 1/10^6 makes the diameter come only from
p, so X minus p would have smaller diameter (the construction fails); an added edge that completes a 6-clique is
found by the clique search and gives only ceil(321/6) = 54 parts.

Usage: python3 -I borsuk_check.py [--controls]   (no download: the construction is the PDF's; under a minute)
"""
import sys
from fractions import Fraction

# ---------- F_16 ----------
EXP = [0] * 30
LOG = [0] * 16
x = 1
for i in range(15):
    EXP[i] = x
    LOG[x] = i
    x <<= 1
    if x & 16:
        x ^= 0b10011            # a^4 = a + 1
for i in range(15, 30):
    EXP[i] = EXP[i - 15]


def mul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]


def frob4(a):                    # a -> a^4
    return mul(mul(a, a), mul(a, a))


def h(u, w):
    return mul(u[0], frob4(w[0])) ^ mul(u[1], frob4(w[1])) ^ mul(u[2], frob4(w[2]))


def points():
    pts = []
    for a in range(16):
        for b in range(16):
            for c in range(16):
                v = (a, b, c)
                if v == (0, 0, 0):
                    continue
                first = next(t for t in v if t)
                if first == 1:
                    pts.append(v)
    return pts


def build(rule_size=3):
    P = points()
    iso = [i for i, u in enumerate(P) if h(u, u) == 0]
    non = [i for i, u in enumerate(P) if h(u, u) != 0]
    isoset = set(iso)
    orth = {i: set(j for j in range(len(P)) if j != i and h(P[i], P[j]) == 0) for i in range(len(P))}
    triples = set()
    for a in non:
        cand = [j for j in orth[a] if j not in isoset]
        for i, b in enumerate(cand):
            for c in cand[i + 1:]:
                if c in orth[b]:
                    triples.add(tuple(sorted((a, b, c))))
    V = sorted(triples)
    isoidx = {q: k for k, q in enumerate(iso)}
    T = []
    for A in V:
        m = 0
        for q in iso:
            if any(q in orth[a] for a in A):
                m |= 1 << isoidx[q]
        T.append(m)
    n = len(V)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if bin(T[i] & T[j]).count("1") == rule_size:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return P, iso, non, orth, V, T, adj


def srg_params(adj):
    n = len(adj)
    degs = {bin(a).count("1") for a in adj}
    lam, mu = set(), set()
    for i in range(n):
        for j in range(i + 1, n):
            c = bin(adj[i] & adj[j]).count("1")
            (lam if (adj[i] >> j) & 1 else mu).add(c)
    return n, degs, lam, mu


def max_clique(adj, stop_at=None):
    """exact maximum clique by bit-set branch and bound; returns (size, one clique)"""
    n = len(adj)
    best = [0, []]

    def expand(R, P):
        if not P:
            if len(R) > best[0]:
                best[0], best[1] = len(R), R[:]
            return
        if len(R) + bin(P).count("1") <= best[0]:
            return
        while P:
            if len(R) + bin(P).count("1") <= best[0]:
                return
            v = P.bit_length() - 1
            P &= ~(1 << v)
            expand(R + [v], P & adj[v])
            if stop_at is not None and best[0] >= stop_at:
                return
    for v in range(n):
        later = adj[v] & ~((1 << (v + 1)) - 1)
        expand([v], later)
        if stop_at is not None and best[0] >= stop_at:
            break
    return best[0], best[1]


def rank_mod(M, p=(1 << 61) - 1):
    M = [[x % p for x in row] for row in M]
    r = 0
    rows, cols = len(M), len(M[0])
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


def components(nodes, adj):
    left = set(nodes)
    out = []
    while left:
        s = left.pop()
        comp, stack = {s}, [s]
        while stack:
            u = stack.pop()
            for w in list(left):
                if (adj[u] >> w) & 1:
                    left.discard(w)
                    comp.add(w)
                    stack.append(w)
        out.append(sorted(comp))
    return out


def gram(u, w, adj):
    if u == w:
        return 90
    return 18 if (adj[u] >> w) & 1 else -6


def split(q0, V, P, orth, non, adj):
    nonset = set(non)
    B = [i for i, A in enumerate(V) if any(a in orth[q0] and a in nonset for a in A)]
    Bs = set(B)
    C = [i for i in range(len(V)) if i not in Bs]
    comps = components(B, adj) if B else []
    return B, C, comps


def degree_data(B, C, comps, adj):
    """the facts of Lemma 1(4): returns (ok, text)"""
    def nb(u, S):
        return sum((adj[u] >> w) & 1 for w in S)
    ok = len(B) == 96 and len(C) == 320 and sorted(len(c) for c in comps) == [32, 32, 32]
    if not ok:
        return False, f"|B| = {len(B)}, |C| = {len(C)}, components {sorted(len(c) for c in comps)}"
    for i, Bi in enumerate(comps):
        for u in Bi:
            if nb(u, Bi) != 20 or nb(u, C) != 80:
                return False, "a B vertex has the wrong degrees"
            for j, Bj in enumerate(comps):
                if j != i and nb(u, Bj) != 0:
                    return False, "an edge between two components"
    for c in C:
        if any(nb(c, Bi) != 8 for Bi in comps) or nb(c, C) != 76:
            return False, "a C vertex has the wrong degrees"
    return True, "|B| = 96 = 3 x 32 (three components), |C| = 320; degrees 20 inside B_i, 0 between, 80 to C; 8 to each B_i and 76 inside C"


# ---------- Q(sqrt 222): numbers a + b sqrt(222), a, b rational ----------
D = 222


class Q222:
    def __init__(self, a, b=0):
        self.a, self.b = Fraction(a), Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, Q222) else Q222(o)
        return Q222(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        o = o if isinstance(o, Q222) else Q222(o)
        return Q222(self.a - o.a, self.b - o.b)

    def __mul__(self, o):
        o = o if isinstance(o, Q222) else Q222(o)
        return Q222(self.a * o.a + D * self.b * o.b, self.a * o.b + self.b * o.a)
    __rmul__ = __mul__

    def sign(self):
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if (a > 0) == (b > 0):
            return 1 if a > 0 else -1
        # opposite signs: compare a^2 with 222 b^2
        if a * a > D * b * b:
            return 1 if a > 0 else -1
        if a * a < D * b * b:
            return 1 if b > 0 else -1
        return 0

    def __repr__(self):
        if self.b == 0:
            return f"{self.a}"
        sgn = "-" if self.b < 0 else "+"
        return f"{self.a} {sgn} {abs(self.b)} sqrt(222)"


def lemma4_5(B, C, comps, adj, t):
    """inner products of z_b and the squared distances in X, for the first vertex b of B_1 and the number t"""
    B1, B2, B3 = comps
    b = B1[0]
    S = [Bi for Bi in comps]
    # x_c . S_i, S_i . S_j, x_b . S_i from the Gram rule (exact integers)
    xS = {c: [sum(gram(c, y, adj) for y in Si) for Si in S] for c in C}
    SS = [[sum(gram(u, w, adj) for u in Si for w in Sj) for Sj in S] for Si in S]
    xbS = [sum(gram(b, y, adj) for y in Si) for Si in S]
    zS = [Fraction(xbS[i]) - Fraction(SS[0][i], 32) for i in range(3)]          # z_b . S_i
    zz = Fraction(90) - Fraction(2, 32) * xbS[0] + Fraction(SS[0][0], 32 * 32)    # |z_b|^2
    zx = {c: Fraction(gram(b, c, adj)) - Fraction(xS[c][0], 32) for c in C}       # z_b . x_c
    pp = t * t * zz                                                                # |p|^2
    dist_p = {c: pp + 90 - 2 * t * zx[c] for c in C}
    return b, xS, SS, xbS, zS, zz, zx, pp, dist_p


def main():
    controls = "--controls" in sys.argv[1:]
    ok = True
    P, iso, non, orth, V, T, adj = build()
    tsz = {bin(m).count("1") for m in T}
    print(f"PG(2,16): {len(P)} points, {len(iso)} isotropic, {len(non)} non-isotropic; vertices (orthogonal "
          f"non-isotropic triples): {len(V)}; |T(A)| = {sorted(tsz)}")
    ok &= len(P) == 273 and len(iso) == 65 and len(V) == 416 and tsz == {15}
    n, degs, lam, mu = srg_params(adj)
    srg = (n, degs, lam, mu) == (416, {100}, {36}, {20})
    print(f"strongly regular (416, 100, 36, 20), every pair checked: {srg} (degrees {degs}, lambda {lam}, mu {mu})")
    ok &= srg
    # Gram matrix G = 96 I + 24 A - 6 J: rank mod 2^61 - 1
    G = [[gram(u, w, adj) for w in range(n)] for u in range(n)]
    rG = rank_mod(G)
    print(f"G = 96I + 24A - 6J (entries 90, 18, -6): rank {rG} modulo 2^61 - 1 (so >= 65 over Q; the SRG eigenvalues "
          f"100, 20, -4 with multiplicities 1, 65, 350 give exactly 65 and G >= 0)")
    ok &= rG == 65
    w, K = max_clique(adj)
    isclique = all((adj[a] >> b) & 1 for i, a in enumerate(K) for b in K[i + 1:])
    print(f"clique number: {w} (a maximum clique {sorted(K)}, a clique: {isclique}); no 6-clique")
    ok &= w == 5 and isclique
    # the split, for every isotropic q0
    good_q0 = 0
    first = None
    for q0 in iso:
        B, C, comps = split(q0, V, P, orth, non, adj)
        okd, txt = degree_data(B, C, comps, adj)
        if okd:
            good_q0 += 1
            if first is None:
                first = (q0, B, C, comps, txt)
    print(f"Lemma 1(3)-(4) degree data hold for {good_q0} of the 65 isotropic points q0")
    ok &= good_q0 == 65
    q0, B, C, comps, txt = first
    print(f"  for q0 = {P[q0]}: {txt}")
    t = Q222(Fraction(-1, 13), Fraction(1, 13))
    rel = t * t * 78 + t * 12
    print(f"t = (sqrt 222 - 1)/13 > 0: {t.sign() == 1}; 78 t^2 + 12 t = {rel.a} + {rel.b} sqrt(222) (must be 102)")
    ok &= rel.a == 102 and rel.b == 0 and t.sign() == 1
    b, xS, SS, xbS, zS, zz, zx, pp, dist_p = lemma4_5(B, C, comps, adj, t)
    orthC = all(v == [0, 0, 0] for v in xS.values())
    det = (SS[0][0] * (SS[1][1] * SS[2][2] - SS[1][2] * SS[2][1]) - SS[0][1] * (SS[1][0] * SS[2][2] - SS[1][2] * SS[2][0])
           + SS[0][2] * (SS[1][0] * SS[2][1] - SS[1][1] * SS[2][0]))
    minor = SS[0][0] * SS[1][1] - SS[0][1] * SS[1][0]
    print(f"Lemma 3: x_c . S_i = 0 for all c in C and i = 1, 2, 3: {orthC}; Gram(S) = {SS}: det {det}, a 2 x 2 minor "
          f"{minor}, so rank 2 and dim W = 63")
    ok &= orthC and det == 0 and minor != 0
    zx_ok = all(zx[c] == (18 if (adj[b] >> c) & 1 else -6) for c in C)
    print(f"Lemma 4 (b = vertex {b} of B_1): z_b . S_i = {[str(v) for v in zS]}, |z_b|^2 = {zz}, z_b . x_c = 18 or -6 as "
          f"b ~ c or not: {zx_ok}")
    ok &= zS == [0, 0, 0] and zz == 78 and zx_ok
    # Gram of {x_c} u {z_b}: integers; rank modulo 2^61 - 1
    Xg = [[gram(u, v, adj) for v in C] + [int(zx[u])] for u in C] + [[int(zx[v]) for v in C] + [int(zz)]]
    rX = rank_mod(Xg)
    print(f"Gram matrix of {{x_c}} u {{z_b}} (321 x 321, integers): rank {rX} modulo 2^61 - 1, so dim span X = 63 "
          f"(with Lemma 3)")
    ok &= rX == 63
    # Lemma 5: squared distances
    vals = {}
    for c in C:
        d = dist_p[c]
        key = "b~c" if (adj[b] >> c) & 1 else "b!~c"
        vals.setdefault(key, set()).add((d.a, d.b))
    d_non = vals["b!~c"]
    d_adj = vals["b~c"]
    nonc = sum(1 for c in C if not (adj[b] >> c) & 1)
    l5 = d_non == {(Fraction(192), Fraction(0))} and len(d_adj) == 1
    da = Q222(*next(iter(d_adj)))
    below = (Q222(192) - da).sign() == 1
    print(f"Lemma 5: |p - x_c|^2 = 192 for the {nonc} c not adjacent to b, = {da} (= 192 - 48t, below 192: {below}) for "
          f"the {len(C) - nonc} neighbours; pairs inside C: 144 or 192; |p|^2 = 78 t^2 != 90: "
          f"{(pp - 90).sign() != 0}; so diam(X)^2 = 192, attained")
    ok &= l5 and below and (pp - 90).sign() != 0 and nonc > 0
    # Lemma 6: smaller-diameter subsets are cliques (with b for p): omega = 5 bounds them
    parts = -(-321 // w)
    print(f"Lemma 6: a smaller-diameter subset is a clique of the graph (b standing for p), so at most {w} points; "
          f"ceil(321/{w}) = {parts} > 64 = 63 + 1: {parts > 64}")
    ok &= parts > 64
    # anchors
    print("anchors (same graph, not the replayed row):")
    print(f"  all 416 points (rank 65): ceil(416/5) = {-(-416 // 5)} > 66, Bondarenko's dimension 65")
    B1, B2, B3 = comps
    pts352 = C + B3
    dif = all(sum(gram(u, y, adj) for y in B1) == sum(gram(u, y, adj) for y in B2) for u in pts352)
    g352 = [[gram(u, v, adj) for v in pts352] for u in pts352]
    r352 = rank_mod(g352)
    print(f"  C u B_3 (352 points) orthogonal to S_1 - S_2: {dif}; Gram rank {r352} mod 2^61 - 1; ceil(352/5) = "
          f"{-(-352 // 5)} > 65, the bound 64 of the [JB2014] row")
    ok &= dif and r352 == 64
    cok = True
    if controls:
        print("\ncontrols (each must come out as stated):")
        sizes = sorted({bin(T[i] & T[j]).count("1") for i in range(n) for j in range(i + 1, n)})
        _, _, _, _, V1, T1, adj1 = build(rule_size=2)
        p1 = srg_params(adj1)
        c = p1 != (416, {100}, {36}, {20})
        print(f"[ctl] the sizes |T(A) & T(A')| that occur are {sizes}; the rule |T & T'| = 2 instead of 3 gives degrees "
              f"{sorted(p1[1])}, lambda {sorted(p1[2])}, mu {sorted(p1[3])}: not SRG(416,100,36,20): {c}")
        cok &= c
        q_bad = non[0]
        Bx, Cx, compsx = split(q_bad, V, P, orth, non, adj)
        okx, txtx = degree_data(Bx, Cx, compsx, adj)
        c = not okx
        print(f"[ctl] a non-isotropic q0 = {P[q_bad]}: {txtx}; the degree data fail: {c}")
        cok &= c
        t2 = t + Fraction(1, 10 ** 6)
        _, _, _, _, _, _, _, _, dist2 = lemma4_5(B, C, comps, adj, t2)
        big = [c_ for c_ in C if (dist2[c_] - 192).sign() > 0]
        c = len(big) > 0
        print(f"[ctl] t + 1/10^6: |p - x_c|^2 > 192 for {len(big)} points, so the diameter comes only from p and "
              f"X minus p (320 points) has smaller diameter: the construction fails: {c}")
        cok &= c
        Kc = sorted(K)
        wv = None
        for v in range(n):
            if v in Kc:
                continue
            hits = [u for u in Kc if (adj[v] >> u) & 1]
            if len(hits) == 4:
                wv = (v, next(u for u in Kc if u not in hits))
                break
        adj2 = adj[:]
        adj2[wv[0]] |= 1 << wv[1]
        adj2[wv[1]] |= 1 << wv[0]
        w2, K2 = max_clique(adj2, stop_at=6)
        c = w2 == 6
        print(f"[ctl] one added edge {wv} completing a 6-clique: the search finds clique number {w2}, so only "
              f"ceil(321/6) = {-(-321 // 6)} parts: {c}")
        cok &= c
        print(f"controls all as expected: {cok}")
    res = ok and cok
    print(f"\nRESULT: {'PASS' if res else 'FAIL'} (321 points in R^63, diameter^2 192, every smaller-diameter subset at "
          f"most 5 points: at least 65 > 64 parts, so Borsuk fails in dimension 63)")
    return 0 if res else 1


if __name__ == "__main__":
    sys.exit(main())
