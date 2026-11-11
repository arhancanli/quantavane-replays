"""Replay of the 9a row [BPZ2026], Theta(C_7) >= 3.258805369885... (Buys, Polak and Zuiddam, arXiv 2607.29681v1,
Section 3, "The cycle C_7"), from the paper's LaTeX source and the base data of their Lean repository
(spectra-research/shannon-capacity-lean @ aa21eeb, ShannonBounds/BaseC7Data.lean, read as data: the lists Ilist,
Xlist and pairs; no Lean proof or program was built or run).

The paper (Section 2): a valid tuple for a graph G is (I, S, f_0, f_1, X) with I, S, X subsets of V(G) and
f_0, f_1 : S -> V(G), such that I is independent, S is a subset of I, for every s in S exactly one of f_0(s), f_1(s)
equals s and N[f_i(s)] meets I only in s (i = 0, 1); f_0(S), f_1(S) and X are independent; and X misses
N[f_0(S)] and N[f_1(S)] at the same time. Its profile is (|I|, |S|, |X|, |X minus (N[f_0(S)] u N[f_1(S)])|).
Gao's product theorem (cited there) gives valid tuples for strong products with profiles combined by the star
operation; the first coordinate is the size of an independent set.

We check, with our own code:
 1. the data: 367 + 367 codes of Z_7^5 (base 7, coordinate k = (n / 7^k) % 7), eight (parent, alternative, side)
    triples; S = the parents, and for each, the side bit says which of f_0, f_1 is the identity (the other gives the
    alternative; exchanging the names of f_0 and f_1 changes neither the conditions nor the profile);
 2. every condition of the definition in C_7^5 (u ~ v iff u != v and every coordinate differs by 0 or +-1 mod 7),
    exhaustively, and the profile, which must be (367, 8, 367, 322) as the paper states;
 3. the paper's schedule pi_2 = pi_1*pi_1, pi_3 = pi_1*pi_2, pi_5 = pi_2*pi_3, pi_10, pi_20, pi_40 with the star
    formula, exactly in integers; then M = first coordinate of pi_40 is an independent set size in C_7^200, and
    floor(M^(1/200) 10^12) / 10^12 = 3.258805369885 is checked exactly by integer powers.
Anchors: the paper's example on C_5 (profile (2,1,2,0), star square (5,0,4,4)); Gao's profile (367, 8, 367, 321) with
the same schedule reproduces Gao's 3.258789153908 (our replay of [Gao2026]).
Not checked: Gao's product theorem (we replayed its k = 2 case for [Gao2026]); the paper's later sections.

Usage: python3 -I bpz_check.py CERTS_ROOT [--controls]
"""
import hashlib
import os
import re
import sys

TEX_SHA = "60bf79af0ab1148b16a81514504d1903f4e5b9add02bf9af53bd186028bbff21"
# per cycle: q, coordinates, data file and SHA-256 (aa21eeb), list names, the paper's profile, schedule, power, bound
CYCLES = {
    7: dict(k=5, file="BaseC7Data.lean", sha="bb50af06aeefdd1a2aafce35802e6f129fa16132b1a070a9d3806a0227f718ba",
            names=("Ilist", "Xlist", "pairs"), profile=(367, 8, 367, 322),
            schedule=[(2, 1, 1), (3, 1, 2), (5, 2, 3), (10, 5, 5), (20, 10, 10), (40, 20, 20)], top=40,
            bound="3.258805369885",
            tex=[r"\(\pi_1=(367,8,367,322)\). Form", r"\(\pi_2=\pi_1\star\pi_1\), \(\pi_3=\pi_1\star\pi_2\),",
                 r"\(\pi_5=\pi_2\star\pi_3\), and then",
                 r"\(\pi_{2r}=\pi_r\star\pi_r\) for \(r=5,10,20\). The resulting \(\pi_{40}\)",
                 r"\(\Theta(C_7)\geq 3.258805369885\ldots\)."]),
    11: dict(k=3, file="BaseC11Data.lean", sha="b6233477afdc35e6cd63a39232dabf1c5f238dfc139767e0bd08f2ca8305e517",
             names=("Ilist", "Xlist", "pairs"), profile=(148, 3, 148, 142),
             schedule=[(2, 1, 1), (3, 1, 2), (4, 2, 2), (7, 3, 4), (8, 4, 4), (14, 7, 7), (15, 7, 8), (29, 14, 15),
                       (58, 29, 29)], top=58, bound="5.294502522149",
             tex=[r"\(\pi_1=(148,3,148,142)\). Form", r"\(\pi_4=\pi_2\star\pi_2\), \(\pi_7=\pi_3\star\pi_4\),",
                  r"\(\pi_8=\pi_4\star\pi_4\), \(\pi_{14}=\pi_7\star\pi_7\),", r"\(\pi_{15}=\pi_7\star\pi_8\),",
                  r"\(\pi_{29}=\pi_{14}\star\pi_{15}\), and", r"\(\pi_{58}=\pi_{29}\star\pi_{29}\). The last profile lies in",
                  r"\(\Theta(C_{11})\geq 5.294502522149\ldots\)."]),
    15: dict(k=4, file="BaseC15Data.lean", sha="a09fc625d159ad083e5c9da35f0707050659133ffa78b23edb1383fe2e07efa5",
             names=("Iraw", "Xraw", "pairsRaw"), profile=(2842, 3, 2842, 2833),
             schedule=[(2 * r, r, r) for r in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)], top=1024,
             bound="7.301600534487",
             tex=[r"\(\pi_1=(2842,3,2842,2833)\). Form", r"\(r=1,2,4,\ldots,512\). The resulting \(\pi_{1024}\) lies in",
                  r"\(\Theta(C_{15})\geq 7.301600534487\ldots\)."]),
    19: dict(k=4, file="BaseC19Data.lean", sha="298ea875746240f154e4ee5461468ea60a7ba468d30cd8d1892c9ee24e8d730b",
             names=("Iraw", "Xraw", "pairsRaw"), profile=(7666, 2, 7666, 7661),
             schedule=[(2 * r, r, r) for r in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)], top=4096,
             bound="9.357192705918",
             tex=[r"\(\pi_1=(7666,2,7666,7661)\). Form", r"\(r=1,2,4,\ldots,2048\). The resulting \(\pi_{4096}\) lies in",
                  r"\(\Theta(C_{19})\geq 9.357192705918\ldots\)."]),
}
STAR_TEX = [r"(a_1-b_1)(a_2-b_2)+c_1b_2+b_1c_2,\ ", r"d_1b_2+b_1d_2,\ ", r"c_1c_2,\ ", r"d_1d_2+(c_1-d_1)(c_2-d_2)"]


def star(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return ((a1 - b1) * (a2 - b2) + c1 * b2 + b1 * c2, d1 * b2 + b1 * d2, c1 * c2, d1 * d2 + (c1 - d1) * (c2 - d2))


def chain(p1, schedule):
    pi = {1: p1}
    for k, i, j in schedule:
        pi[k] = star(pi[i], pi[j])
    return pi


def trunc_root(M, p, digits):
    """floor(M^(1/p) * 10^digits), exactly: the largest integer x with x^p <= M * 10^(digits p)"""
    target = M * 10 ** (digits * p)
    lo, hi = 0, 1
    while hi ** p <= target:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** p <= target:
            lo = mid
        else:
            hi = mid
    return lo


def parse_list(text, name):
    m = re.search(r"def " + name + r" : List (?:Code|Nat) :=\s*\[(.*?)\]", text, re.S)
    return [int(x) for x in re.findall(r"\d+", m.group(1))]


def parse_pairs(text, name):
    m = re.search(r"def " + name + r" : List \((?:Code|Nat) × (?:Code|Nat) × Bool\) :=\s*\[(.*?)\]\n", text, re.S)
    return [(int(a), int(b), c == "true") for a, b, c in re.findall(r"\((\d+), (\d+), (true|false)\)", m.group(1))]


def decode(n, k=5, q=7):
    return tuple((n // q ** i) % q for i in range(k))


def adjacent_or_equal(u, v, q=7):
    return all((a - b) % q in (0, 1, q - 1) for a, b in zip(u, v))


def closed_nbhd(v, q=7):
    """all words within +-1 (mod q) in every coordinate"""
    out = [()]
    for x in v:
        out = [w + ((x + d) % q,) for w in out for d in (-1, 0, 1)]
    return set(out)


def independent(S, q=7):
    """no two distinct words within +-1 in every coordinate (each word's closed neighbourhood looked up)"""
    S = list(S)
    Sset = set(S)
    if len(Sset) != len(S):
        return False
    for v in S:
        for w in closed_nbhd(v, q):
            if w != v and w in Sset:
                return False
    return True


def check_tuple(I, S, f0, f1, X, q=7):
    """the conditions of the paper's definition; returns (ok, list of failures, profile)"""
    fails = []
    Iset = set(I)
    if len(Iset) != len(I):
        fails.append("I has repeated words")
    if not independent(I, q):
        fails.append("I is not independent")
    if not set(S) <= Iset:
        fails.append("S is not a subset of I")
    for s in S:
        if (f0[s] == s) == (f1[s] == s):
            fails.append(f"not exactly one of f_0(s), f_1(s) equals s at {s}")
        for f in (f0, f1):
            meet = closed_nbhd(f[s], q) & Iset
            if meet != {s}:
                fails.append(f"N[f(s)] meets I in {len(meet)} words at {s}")
    F0 = {f0[s] for s in S}
    F1 = {f1[s] for s in S}
    for name, T in (("f_0(S)", F0), ("f_1(S)", F1), ("X", X)):
        if not independent(T, q):
            fails.append(f"{name} is not independent")
    if len(set(X)) != len(X):
        fails.append("X has repeated words")
    N0 = set().union(*(closed_nbhd(v, q) for v in F0)) if F0 else set()
    N1 = set().union(*(closed_nbhd(v, q) for v in F1)) if F1 else set()
    both = [x for x in X if x in N0 and x in N1]
    if both:
        fails.append(f"{len(both)} words of X lie in N[f_0(S)] and N[f_1(S)]")
    d = len([x for x in X if x not in N0 and x not in N1])
    h0 = len([x for x in X if x in N0])
    h1 = len([x for x in X if x in N1])
    return not fails, fails, (len(I), len(S), len(X), d), (h0, h1)


def load_cycle(base, n, C):
    raw = open(os.path.join(base, C["file"]), "rb").read()
    h = hashlib.sha256(raw).hexdigest()
    text = raw.decode("utf-8")
    nI, nX, nP = C["names"]
    Ic, Xc, pairs = parse_list(text, nI), parse_list(text, nX), parse_pairs(text, nP)
    q, k = n, C["k"]
    in_range = all(0 <= v < q ** k for v in Ic + Xc + [x for pp in pairs for x in pp[:2]])
    I = [decode(v, k, q) for v in Ic]
    X = [decode(v, k, q) for v in Xc]
    S, f0, f1 = [], {}, {}
    for par, alt, side in pairs:
        sw = decode(par, k, q)
        S.append(sw)
        a = decode(alt, k, q)
        f1[sw], f0[sw] = (sw, a) if side else (a, sw)
    return h, in_range, I, S, f0, f1, X, len(Ic), len(Xc), len(pairs)


def main():
    root = sys.argv[1]
    controls = "--controls" in sys.argv[2:]
    ok = True
    base = os.path.join(root, "bpz2026")
    tex = open(os.path.join(base, "main3.tex"), "rb").read()
    th = hashlib.sha256(tex).hexdigest()
    print("9a [BPZ2026] (Buys, Polak and Zuiddam, arXiv 2607.29681v1): the C_7 row Theta(C_7) >= 3.258805369885..., "
          "and the page's note on C_11, C_15 and C_19")
    print(f"bpz2026/main3.tex sha256 {th}: as pinned {th == TEX_SHA}")
    ok &= th == TEX_SHA
    t = tex.decode("utf-8")
    for s_ in STAR_TEX:
        v = s_ in t
        print(f"  our transcription appears verbatim: {v}   {s_}")
        ok &= v
    ex_ok = star((2, 1, 2, 0), (2, 1, 2, 0)) == (5, 0, 4, 4)
    v5, _, p5, _ = check_tuple([(0,), (2,)], [(0,)], {(0,): (0,)}, {(0,): (4,)}, [(1,), (3,)], q=5)
    print(f"anchor: the paper's C_5 example is a valid tuple with profile {p5}: {v5 and p5 == (2, 1, 2, 0)}; its star "
          f"square is (5, 0, 4, 4): {ex_ok}")
    ok &= ex_ok and v5 and p5 == (2, 1, 2, 0)
    saved = {}
    for n, C in CYCLES.items():
        print(f"\nC_{n}: {C['file']} (spectra-research/shannon-capacity-lean @ aa21eeb)")
        for s_ in C["tex"]:
            v = s_ in t
            print(f"  our transcription appears verbatim: {v}   {s_}")
            ok &= v
        h, in_range, I, S, f0, f1, X, nI, nX, nP = load_cycle(base, n, C)
        print(f"  sha256 {h}: as pinned {h == C['sha']}; |I| = {nI}, |X| = {nX}, {nP} (parent, alternative, side) "
              f"triples; codes in range: {in_range}")
        ok &= h == C["sha"] and in_range
        v_ok, fails, prof, hv = check_tuple(I, S, f0, f1, X, q=n)
        print(f"  valid tuple in C_{n}^{C['k']} (every condition of Section 2): {v_ok}" + ("" if v_ok else f"; {fails[:3]}"))
        print(f"  profile {prof} (the paper: {C['profile']}): {prof == C['profile']}; footprints in X: {hv}")
        ok &= v_ok and prof == C["profile"]
        pi = chain(prof, C["schedule"])
        M = pi[C["top"]][0]
        p = C["k"] * C["top"]
        x = trunc_root(M, p, 12)
        val = f"{x // 10 ** 12}.{x % 10 ** 12:012d}"
        print(f"  pi_{C['top']}: an independent set in C_{n}^{p} of size M ({len(str(M))} digits); floor(M^(1/{p}) "
              f"10^12)/10^12 = {val}: the stated {C['bound']}: {val == C['bound']}")
        ok &= val == C["bound"]
        saved[n] = (I, S, f0, f1, X, M, p)
    gao = chain((367, 8, 367, 321), CYCLES[7]["schedule"])
    xg = trunc_root(gao[40][0], 200, 12)
    gv = f"{xg // 10 ** 12}.{xg % 10 ** 12:012d}"
    print(f"\nanchor: Gao's profile (367, 8, 367, 321) on the C_7 schedule gives {gv} (the [Gao2026] row: "
          f"3.258789153908): {gv == '3.258789153908'}")
    ok &= gv == "3.258789153908"
    cok = True
    if controls:
        print("\ncontrols (each must come out as stated):")
        I, S, f0, f1, X, M, p = saved[7]
        X2 = list(X)
        X2[0] = S[0]
        v2, f2, _, _ = check_tuple(I, S, f0, f1, X2)
        c = not v2
        print(f"[ctl] C_7: one word of X replaced by a parent (in both neighbourhoods): refused ({f2[-1] if f2 else ''}): {c}")
        cok &= c
        f0b, f1b = dict(f0), dict(f1)
        f0b[S[0]] = f1b[S[0]] = S[0]
        v3, f3, _, _ = check_tuple(I, S, f0b, f1b, X)
        c = not v3
        print(f"[ctl] C_7: both maps the identity at one parent: refused ({f3[0] if f3 else ''}): {c}")
        cok &= c
        I2 = list(I)
        I2[5] = tuple((c_ + 1) % 7 if k == 0 else c_ for k, c_ in enumerate(I2[6]))
        v4, f4, _, _ = check_tuple(I2, S, f0, f1, X)
        c = not v4
        print(f"[ctl] C_7: one word of I moved next to another: refused ({f4[0] if f4 else ''}): {c}")
        cok &= c
        I, S, f0, f1, X, M, p = saved[19]
        X3 = list(X)
        X3[1] = X3[0]
        v5b, f5, _, _ = check_tuple(I, S, f0, f1, X3, q=19)
        c = not v5b
        print(f"[ctl] C_19: a repeated word in X: refused ({f5[0] if f5 else ''}): {c}")
        cok &= c
        for n, C in CYCLES.items():
            M, p = saved[n][5], saved[n][6]
            y = trunc_root(M, p, 12) + 1
            c = y ** p > M * 10 ** (12 * p)
            print(f"[ctl] C_{n}: the stated value plus 10^-12 is not implied: {c}")
            cok &= c
        print(f"controls all as expected: {cok}")
    res = ok and cok
    print(f"\nRESULT: {'PASS' if res else 'FAIL'} (four valid base tuples with the paper's profiles; the schedules give "
          f"Theta(C_7) >= 3.258805369885, Theta(C_11) >= 5.294502522149, Theta(C_15) >= 7.301600534487, "
          f"Theta(C_19) >= 9.357192705918)")
    return 0 if res else 1


if __name__ == "__main__":
    sys.exit(main())
