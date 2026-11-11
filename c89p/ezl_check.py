#!/usr/bin/env python3
"""Exact replay of two more 89a upper-bound rows (teorth/optimizationproblems, main 87b4ebb; both added in PR #255),
from the constructions as Even-Zohar and Linial describe them (arXiv 1312.1205v2, section 2, read as text):

  [Tho1997]  0.030291: Thomason's product M_4 (x) K_4 (x) G_18, G_18 = (K_3 (x) K_3) (.) K_2, for which EZL state
             R(K_4 + A_4, .) = 3769/124416 ~ 1/33.0103 (and that Thomason showed C_89 strictly below this value).
  [EZL2015]  0.030285: the nested blow-up, R(K_4 + A_4, M_4 (x) K_4 (x) (K_3 (x) K_3)^(.)n) -> 1411/46592 ~ 1/33.0205.

Definitions (EZL): G (x) H, the tensor product, has (g,h) ~ (g',h') iff exactly one of g ~ g', h ~ h' holds (loops
count, none here); G (.) H, the composition, has (g,h) ~ (g',h') iff g ~ g', or g = g' and h ~ h'; R(K_4 + A_4, G) is
the probability that four vertices sampled with replacement span a K_4 or an independent 4-set (a repeated vertex is a
non-edge), the limit density of monochromatic K_4 in the blow-ups of G; M_4 is a perfect matching on 4 vertices.

Computed (exact rationals, plus ../c89p/k4hom.c on explicit graphs):
  1. labelled 4-point profiles r(G) (64 patterns on the 6 pairs) by enumerating all 4-tuples, for K_4, M_4, K_3 (x) K_3
     (checked to be the Paley graph of order 9, SRG(9,4,1,2)) and G_18; the tensor product as the XOR convolution of
     profiles (EZL's Corollary 1, after Thomason) and, independently, the composition by our own rule (classes of equal
     outer coordinates, inner coordinates independent), checked against direct enumeration on G_18 = (K_3 (x) K_3) (.) K_2;
  2. R for M_4 (x) K_4 (x) K_3 (x) K_3 (EZL: 11411/373248) and for M_4 (x) K_4 (x) G_18 (3769/124416), each also by
     k4hom on the explicit graph (144 and 288 vertices: hom(K_4) + hom(K_4, looped complement), two routes inside);
  3. the nested blow-up H^(.)n, H = K_3 (x) K_3: its profile r_n satisfies r_n = Phi(r_{n-1}), where the outer tuple
     with all four coordinates equal (probability 9^-3) reproduces the inner pattern and every other term uses marginals
     on at most three points; so the limit is solved level by level (pairs, triples, quadruples) exactly, and R for
     M_4 (x) K_4 (x) H^(.)infinity is compared with 1411/46592; R for n = 1..6 is listed (it decreases to the limit).
Notes: 1411/46592 = 0.03028416...: the row's 0.030285 is it rounded up (valid). 3769/124416 = 0.03029353...: the
[Tho1997] row prints 0.030291, below this; EZL report Thomason's bound as strictly below 3769/124416, the referee's
perturbation improving it by less than 10^-7, so 0.030291 is not reached by the construction as described (Thomason's
paper itself we could not read).
Red controls (--controls): K_2 in G_18 replaced by its complement (no edge) changes the value; the composition order
reversed (K_2 (.) (K_3 (x) K_3)) gives a different value; a claimed 0.030284 is below 1411/46592.
usage: python3 -I ezl_check.py CERTS_ROOT BIN_DIR [--controls]
"""
import itertools
import os
import subprocess
import sys
from fractions import Fraction

ROOT, BIN = sys.argv[1], sys.argv[2]
CONTROLS = "--controls" in sys.argv
PAIRS = list(itertools.combinations(range(4), 2))
BIT = {p: 1 << i for i, p in enumerate(PAIRS)}


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


class Graph:
    def __init__(self, n, adj):
        self.n = n
        self.adj = adj          # adj[u][v] in {0, 1}; adj[v][v] = 1 means a loop


def complete(n):
    return Graph(n, [[int(i != j) for j in range(n)] for i in range(n)])


def matching4():
    return Graph(4, [[int({i, j} in ({0, 1}, {2, 3})) for j in range(4)] for i in range(4)])


def edge2():
    return complete(2)


def tensor(G, H):
    V = [(g, h) for g in range(G.n) for h in range(H.n)]
    return Graph(len(V), [[G.adj[a][c] ^ H.adj[b][d] for (c, d) in V] for (a, b) in V])


def compose(G, H):
    V = [(g, h) for g in range(G.n) for h in range(H.n)]
    return Graph(len(V), [[1 if (a != c and G.adj[a][c]) else (H.adj[b][d] if a == c else 0) for (c, d) in V]
                          for (a, b) in V])


def profile(G):
    n = G.n
    cnt = [0] * 64
    for t in itertools.product(range(n), repeat=4):
        m = 0
        for (i, j), b in BIT.items():
            if G.adj[t[i]][t[j]]:
                m |= b
        cnt[m] += 1
    return [Fraction(c, n ** 4) for c in cnt]


def xor_conv(r, s):
    out = [Fraction(0)] * 64
    for a in range(64):
        if r[a]:
            for b in range(64):
                if s[b]:
                    out[a ^ b] += r[a] * s[b]
    return out


def R(r):
    return r[63] + r[0]


def marginal(r, S):
    """Distribution of the pattern on the pairs inside S (a tuple of indices), as a dict {mask restricted: prob}."""
    keep = sum(BIT[p] for p in PAIRS if p[0] in S and p[1] in S)
    out = {}
    for m in range(64):
        if r[m]:
            out[m & keep] = out.get(m & keep, 0) + r[m]
    return out


def relabel(mask, S):
    """A pattern on the pairs of S (indices of 0..3) expressed with S's elements replaced by 0..|S|-1 in order."""
    pos = {v: i for i, v in enumerate(S)}
    out = 0
    for (i, j), b in BIT.items():
        if mask & b:
            out |= BIT[(pos[i], pos[j])]
    return out


def outer_structures(G):
    """{(partition as tuple of classes, cross pattern mask): probability} over 4-tuples of G sampled with replacement."""
    n = G.n
    out = {}
    for t in itertools.product(range(n), repeat=4):
        classes = {}
        for i, v in enumerate(t):
            classes.setdefault(v, []).append(i)
        part = tuple(sorted(tuple(c) for c in classes.values()))
        m = 0
        for (i, j), b in BIT.items():
            if t[i] != t[j] and G.adj[t[i]][t[j]]:
                m |= b
        key = (part, m)
        out[key] = out.get(key, 0) + Fraction(1, n ** 4)
    return out


def compose_profile(Gouter, inner_marg):
    """Profile of Gouter (.) Inner, given inner_marg(k) = labelled k-point profile of Inner as {mask on 0..k-1: p}."""
    res = [Fraction(0)] * 64
    for (part, cross), w in outer_structures(Gouter).items():
        dists = [{0: Fraction(1)}]
        for S in part:
            if len(S) == 1:
                continue
            md = inner_marg(len(S))
            new = []
            for d in dists:
                nd = {}
                for m0, p0 in d.items():
                    for mk, pk in md.items():
                        # map the pattern on 0..k-1 back to the indices S
                        mm = 0
                        for (i, j), b in BIT.items():
                            if i < len(S) and j < len(S) and mk & b:
                                mm |= BIT[(S[i], S[j])]
                        nd[m0 | mm] = nd.get(m0 | mm, 0) + p0 * pk
                new.append(nd)
            dists = new
        for m, p in dists[-1].items():
            res[cross | m] += w * p
    return res


def k_marginal_from_profile(r):
    """inner_marg(k) for a 4-point profile r: the marginal on the first k points, relabelled 0..k-1."""
    def f(k):
        d = marginal(r, tuple(range(k)))
        return {relabel(m, tuple(range(k))): p for m, p in d.items()}
    return f


def nested_limit(H):
    """Labelled 4-point profile of H^(.)infinity, solved level by level: r_k = c_k / (1 - h^(1-k))."""
    h = H.n
    outer = outer_structures(H)
    sol = {1: {0: Fraction(1)}}
    for k in (2, 3, 4):
        # outer structures restricted to the first k points
        acc = {}
        for (part, cross), w in outer.items():
            sub = tuple(tuple(i for i in S if i < k) for S in part)
            sub = tuple(S for S in sub if S)
            crossk = cross & sum(BIT[p] for p in PAIRS if p[1] < k)
            key = (sub, crossk)
            acc[key] = acc.get(key, 0) + w
        c = {}
        diag = Fraction(0)
        for (part, cross), w in acc.items():
            if len(part) == 1 and len(part[0]) == k:
                diag += w           # all k equal: the inner pattern on k points again
                continue
            dists = [{0: Fraction(1)}]
            for S in part:
                if len(S) == 1:
                    continue
                md = sol[len(S)]
                nd = {}
                for m0, p0 in dists[-1].items():
                    for mk, pk in md.items():
                        mm = 0
                        for (i, j), b in BIT.items():
                            if j < len(S) and mk & b:
                                mm |= BIT[(S[i], S[j])]
                        nd[m0 | mm] = nd.get(m0 | mm, 0) + p0 * pk
                dists.append(nd)
            for m, p in dists[-1].items():
                c[cross | m] = c.get(cross | m, 0) + w * p
        need(diag == Fraction(1, h ** (k - 1)), "unexpected diagonal weight")
        sol[k] = {m: v / (1 - diag) for m, v in c.items()}
    r = [Fraction(0)] * 64
    for m, v in sol[4].items():
        r[m] += v
    return r


def to_graph6(G):
    n = G.n
    bits = [G.adj[i][j] for j in range(1, n) for i in range(j)]
    bits += [0] * (-len(bits) % 6)
    head = bytes([126, 63 + (n >> 12 & 63), 63 + (n >> 6 & 63), 63 + (n & 63)]) if n > 62 else bytes([63 + n])
    return head + bytes(63 + int("".join(map(str, bits[k:k + 6])), 2) for k in range(0, len(bits), 6))


def k4hom_density(G, name):
    path = os.path.join(BIN, name + ".g6")
    open(path, "wb").write(to_graph6(G))
    out = subprocess.run([os.path.join(BIN, "k4hom"), path], capture_output=True, text=True, check=True).stdout
    d = {l.split()[0]: l.split()[1:] for l in out.splitlines()}
    Ar, Ab = int(d["routeA"][1]), int(d["routeA"][3])
    br = [int(x) for x in d["routeB"][1:5]]
    bb = [int(x) for x in d["routeB"][6:10]]
    need(Ar == 24 * br[3] and Ab == bb[0] + 14 * bb[1] + 36 * bb[2] + 24 * bb[3], "k4hom routes disagree")
    return Fraction(Ar + Ab, G.n ** 4)


def main():
    ok = True
    tex = open(os.path.join(ROOT, "ezl", "ind.tex"), encoding="utf-8").read()
    for s in ("\\frac{11411}{373248}", "\\frac{3769}{124416}", "\\frac{1411}{46592}", "G_{18} = (K_3 \\otimes K_3) \\odot K_2"):
        need(s in tex, f"not in EZL's source: {s}")
    K3, K4, M4, K2 = complete(3), complete(4), matching4(), edge2()
    P9 = tensor(K3, K3)
    deg = {sum(row) for row in P9.adj}
    lam = {sum(P9.adj[a][x] * P9.adj[b][x] for x in range(9)) for a in range(9) for b in range(9) if a != b and P9.adj[a][b]}
    mu = {sum(P9.adj[a][x] * P9.adj[b][x] for x in range(9)) for a in range(9) for b in range(9) if a != b and not P9.adj[a][b]}
    paley = (deg, lam, mu) == ({4}, {1}, {2}) and all(P9.adj[v][v] == 0 for v in range(9))
    G18 = compose(P9, K2)
    rK4, rM4, rP9, rG18 = profile(K4), profile(M4), profile(P9), profile(G18)
    rG18c = compose_profile(P9, k_marginal_from_profile(profile(K2)))
    comp_ok = rG18c == rG18
    v_a = R(xor_conv(xor_conv(rM4, rK4), rP9))
    v_t = R(xor_conv(xor_conv(rM4, rK4), rG18))
    big_a = tensor(tensor(M4, K4), P9)
    big_t = tensor(tensor(M4, K4), G18)
    k_a, k_t = k4hom_density(big_a, "m4k4p9"), k4hom_density(big_t, "m4k4g18")
    rinf = nested_limit(P9)
    v_inf = R(xor_conv(xor_conv(rM4, rK4), rinf))
    print(f"K_3 (x) K_3: SRG(9,4,1,2), loopless (the Paley graph of order 9): {paley}; G_18 = (K_3 (x) K_3) (.) K_2 has "
          f"{G18.n} vertices; composition rule = direct enumeration on G_18: {comp_ok}")
    print(f"  R(K_4 + A_4, M_4 (x) K_4 (x) K_3 (x) K_3) = {v_a} (profiles), {k_a} (k4hom, 144 vertices); EZL: 11411/373248")
    print(f"  R(K_4 + A_4, M_4 (x) K_4 (x) G_18) = {v_t} = {float(v_t):.12f} (profiles), {k_t} (k4hom, 288 vertices); "
          f"EZL: 3769/124416")
    seq = []
    rn = rP9
    for n in range(1, 7):
        seq.append(R(xor_conv(xor_conv(rM4, rK4), rn)))
        rn = compose_profile(P9, k_marginal_from_profile(rn))
    print("  nested blow-ups H^(.)n, H = K_3 (x) K_3: R(K_4 + A_4, M_4 (x) K_4 (x) H^(.)n) for n = 1..6: "
          + ", ".join(f"{float(x):.10f}" for x in seq))
    print(f"  limit (solved level by level): {v_inf} = {float(v_inf):.12f}; EZL: 1411/46592")
    good_t = paley and comp_ok and v_a == k_a == Fraction(11411, 373248) and v_t == k_t == Fraction(3769, 124416)
    good_e = good_t and v_inf == Fraction(1411, 46592) and Fraction(1411, 46592) <= Fraction(30285, 10 ** 6)
    print(f"RESULT [EZL2015]: the nested blow-ups give C_89 <= 1411/46592 = {float(v_inf):.10f}...; the printed 0.030285 "
          f"is it rounded up: {'HOLDS' if good_e else 'FAILS'}")
    below = Fraction(30291, 10 ** 6) < Fraction(3769, 124416)
    print(f"RESULT [Tho1997]: the product M_4 (x) K_4 (x) G_18 gives exactly 3769/124416 = {float(v_t):.12f}... (two routes, "
          f"as EZL state); the printed 0.030291 is {'BELOW' if below else 'not below'} it by "
          f"{float(Fraction(3769, 124416) - Fraction(30291, 10 ** 6)):.2e}, so the construction as described supports "
          f"0.030294 (rounded up), not 0.030291: {'CONSTRUCTION HOLDS, PRINTED VALUE NOT REACHED' if good_t and below else 'CHECK'}")
    ok &= good_e and good_t
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        Kbar = Graph(2, [[0, 0], [0, 0]])
        v1 = R(xor_conv(xor_conv(rM4, rK4), profile(compose(P9, Kbar))))
        res.append(("K_2 replaced by its complement", v1 != Fraction(3769, 124416)))
        v2 = R(xor_conv(xor_conv(rM4, rK4), profile(compose(K2, P9))))
        res.append(("composition order reversed", v2 != Fraction(3769, 124416)))
        print(f"  (K_3 (x) K_3) (.) (2 isolated vertices): {float(v1):.10f}; K_2 (.) (K_3 (x) K_3): {float(v2):.10f}")
        res.append(("a claimed 0.030284", Fraction(30284, 10 ** 6) < Fraction(1411, 46592)))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
