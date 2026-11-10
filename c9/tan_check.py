"""Replay of the arithmetic and rule layer of [Tan2026], Theta(C_7) >= M_*^(1/500) = 3.25883262... (Tandon,
"Strengthening recursive constructions for zero-error Shannon capacity", arXiv:2608.30273v1, Theorem 2), the 9a row,
written from the paper alone (its verification package was neither read nor run); and the checks we can make on
PR #296 (LT2026, Theta(C_7) >= 3.258838705444608*) from its certificate data alone. Standard library only.

What is computed (definitions from the paper):
  * profiles pi = (a, t, s, o, h, v); sigma swaps h and v (Sec. 4.1);
  * Gao's product (Lemma 2.3): a = (a_L - t_L)(a_R - t_R) + t_L s_R + s_L t_R, t = t_L o_R + o_L t_R, and for the
    auxiliary set X_L x X_R: s = s_L s_R, o = o_L o_R + (s_L - o_L)(h_R + v_R), h = h_L o_R + o_L v_R,
    v = v_L o_R + o_L h_R (the case J_0 = J_H = J_V = X_L of Theorem 1, with q = s_L - o_L because an independent set
    X meets N(X^0) exactly in X^0);
  * the heterogeneous product (Theorem 1) from the codebook parameters (j_0, o_0, h_0, v_0), (j_H, q_H), (j_V, q_V);
  * Phi(a, t, s, o, h, v) = (a - t, o, v, h, t, t, t) in the order (B, N, A, D, O, H, V) (eq. (Gao-to-BPZ));
  * the ternary rule S_3b and the terminal code K_4a as transcribed in Appendix A, with the separation relation of
    Definition 4.2: S_3b checked admissible (Definition 4.4 (i) and (ii)) and K_4a checked pairwise separated;
  * Tan's chain of Theorem 2 (Steps 1-5, codebooks as in Table 6), from BPZ's base profile (367, 8, 367, 322, 26, 19)
    and Tan's four certified counts C1 = 27488, C2 = 12839823, C3 = 14045805, C4 = 841760069965664 (inputs here; not
    recomputed), every displayed profile compared, then n_25 = S_3b(n_6, n_11, n_8) and M_* = K_4a(n_25, x4).
For PR #296 (certificate.json, its stated M, d = 680): the two comparisons the pull request states, M^500 > M_*^680 and
M 10^(8*680) > 325883262^680, and the printed bound: 3.258838705444608^680 <= M exactly. The construction of M
itself is described in the pattern notation of Tan's code (N<...>, S<...>, P15A8, ...), which the paper does not
define; our natural reading of it is evaluated and printed, and it does NOT give the stated M (see the output).

Usage: python3 -I tan_check.py [CERTIFICATE_JSON] [--tex=TEX_DIR] [--controls]
"""
import hashlib
import json
import sys
from decimal import Decimal, getcontext

L = "BNADOHV"
SEP = {"B": "OHV", "N": "ADOHV", "A": "NDH", "D": "NAV", "O": "BN", "H": "BNA", "V": "BND"}   # Definition 4.2
S3B = {   # Appendix A.2, S_3b (45 words)
    "B": "BBB BHA NAV NDH NVD ABH AHN AHD DBV DVN DVD HND HDN HDD VNA VAB VDA",
    "N": "BVD NNN NVA ANH ADH AHN DNV DDH DHN HVD",
    "A": "NBA NAN ABB AAV AHA",
    "D": "NBD NDN DBB DAV DHA",
    "O": "",
    "H": "NNH NHN HNN HVV",
    "V": "NNV NVN VNN VVV",
}
K4A = ("BBBB BNVB BAHB BDVB BHNB BHDB BVAB NBBV NNVH NAVH NDHH NHNH NHAH NVDV ABBV ANHH AAVH ADHH AHAH AVNH AVDH "
       "DBBH DNVV DAVV DDHV DHNV DHAV DVDV HBBA HNVD HAHD HAVN HDHN HDVD HHNN HHND HHAN HHDD HVAD HVDN VBBN VBBD "
       "VNVN VNVA VAHA VDVA VHDA VVNA VVAA").split()   # Appendix A.3, K_4a (49 words)
BASE = (367, 8, 367, 322, 26, 19)          # BPZ base gadget (Table 4 of the paper; Lean-verified upstream)
C1, C2, C3, C4 = 27488, 12839823, 14045805, 841760069965664   # Tan's certified counts (Table 5), inputs here
PAPER = {   # displayed profiles of the proof of Theorem 2, for comparison
    "G10": (134753, 5152, 134689, 105709, 14490, 14490), "G10A": (134753, 5152, 134689, 105709, 12236, 16744),
    "G10D": (134753, 5152, 134689, 105709, 16744, 12236),
    "G15X": (49495055, 2504616, 49430863, 35342398, 6674251, 7414214),
    "G15het": (49495055, 2504616, 49433743, 35275258, 6703815, 7454670),
    "G15AX": (49495055, 2504616, 49430863, 35342398, 5948463, 8140002),
    "G15Ahet": (49495055, 2504616, 49432527, 35303606, 5948463, 8180458),
    "G15DX": (49495055, 2504616, 49430863, 35342398, 8140002, 5948463),
    "G15Dhet": (49495055, 2504616, 49432527, 35303606, 8180458, 5948463),
    "G30_6": (2455719231434017, 176870111100096, 2444563032022369, 1426136268314719, 509731462042710,
              508695301664940),
    "G25_8": (6682019915439, 446844487240, 6659829690543, 4104006957194, 1077319494019, 1478503239330),
    "G40_8": (331761855244358723969, 26071517201241409024, 329305968864222045505, 179001784071649220449,
              87767535795425989412, 62536648997146835644),
    "G30L": (2455716185820961, 176941111954464, 2443492469880801, 1446768903083453, 420235140891852,
             576488425905496),
    "Ghat30": (2455716185820961, 176941111954464, 2444076335041121, 1435184942161465, 420235140891852,
               588656251987804),
    "Ghat25": (6682019915439, 446844487240, 6659079476911, 4118763352946, 1061812884635, 1478503239330),
    "G55": (16478688413981213775272008847, 1375259276200664518625890664, 16304893524159967980117037071,
            8049523822718249447032963650, 3444579390015851502852512947, 4810790311425867030231560474),
}
N25 = (12651555102711866584006925012965893561076969596025225905818339799,
       3066132745665704601020173586804414424611769880543231510544034510,
       4664412231532807578550047842170656906931895714626194718133859805,
       5561607824514622663538712910059154540700557990699908637174352708, 0,
       911561283323995002333709184151818987697897558779375600768623640,
       911561283323995002333709184151818987697897558779375600768623640)
MSTAR = int("3396467291811745694340851757371833322531739075465110" "0605971152590324909277046876400440353958428359262989"
            "6215203947372766811966557767583863304708327733085591" "5999031707610213141337406106659583482352927456098895"
            "0698158729501949866486079791784996498585401881281")
CERT_SHA = "8d55f7465bdedf6b1ea9c4ad3504c1b68c93705500590ff11b6fca66cb9d1f2a"   # as PR #296 states


def sep(x, y):
    return y in SEP[x]


def sigma(g):
    a, t, s, o, h, v = g
    return (a, t, s, o, v, h)


def gao(Lg, Rg):
    aL, tL, sL, oL, hL, vL = Lg
    aR, tR, sR, oR, hR, vR = Rg
    return ((aL - tL) * (aR - tR) + tL * sR + sL * tR, tL * oR + oL * tR, sL * sR,
            oL * oR + (sL - oL) * (hR + vR), hL * oR + oL * vR, vL * oR + oL * hR)


def hetgao(Lg, Rg, cb):
    aL, tL, sL, oL, hL, vL = Lg
    aR, tR, sR, oR, hR, vR = Rg
    j0, o0, h0, v0, jH, qH, jV, qV = cb
    if j0 != o0 + h0 + v0 or not (0 <= qH <= jH and 0 <= qV <= jV):
        raise ValueError("inconsistent codebook parameters")
    return ((aL - tL) * (aR - tR) + tL * sR + sL * tR, tL * oR + oL * tR, j0 * oR + jH * hR + jV * vR,
            o0 * oR + qH * hR + qV * vR, h0 * oR + (jV - qV) * vR, v0 * oR + (jH - qH) * hR)


def aux(g):
    """J0 = the auxiliary set of a gadget with the same code, private pairs and transversals: its (s, o, h, v)"""
    return g[2:]


def own(g):
    """a one-sided codebook equal to the left gadget's own auxiliary set: (s, s - o)"""
    return (g[2], g[2] - g[3])


def phi(g):
    a, t, s, o, h, v = g
    return (a - t, o, v, h, t, t, t)


def apply_words(words, ns):
    tot = 0
    for w in words:
        p = 1
        for i, c in enumerate(w):
            p *= ns[i][L.index(c)]
        tot += p
    return tot


def rule(T, ns):
    return tuple(apply_words(T[lam].split(), ns) for lam in L)


def admissible(T, m):
    words = {lam: T[lam].split() for lam in L}
    i_ok = all(any(sep(x[k], y[k]) for k in range(m)) for lam in L for x in words[lam] for y in words[lam] if x != y)
    ii_ok = all(any(sep(x[k], y[k]) for k in range(m)) for lam in L for mu in L if sep(lam, mu)
                for x in words[lam] for y in words[mu])
    return i_ok and ii_ok


def separated(K):
    return len(set(K)) == len(K) and all(any(sep(x[k], y[k]) for k in range(len(x))) for x in K for y in K if x != y)


def chain(base=BASE, c=(C1, C2, C3, C4), T=S3B, K=K4A):
    c1, c2, c3, c4 = c
    G5 = base
    G = {"G10": gao(G5, G5), "G10A": gao(sigma(G5), G5), "G10D": gao(G5, sigma(G5))}
    Jp = (G["G10"][0], c1)                                         # J+ : |J+| = a_10, q = C1
    G["G15X"] = gao(G["G10"], G5)
    G["G15het"] = hetgao(G["G10"], G5, aux(G["G10"]) + Jp + Jp)                       # [(X10, J+, J+); X5]
    G["G15AX"] = gao(G["G10A"], G5)
    G["G15Ahet"] = hetgao(G["G10A"], G5, aux(G["G10A"]) + Jp + own(G["G10A"]))       # [(X10, J+, X10); X5]
    G["G15DX"] = gao(G["G10D"], sigma(G5))
    G["G15Dhet"] = hetgao(G["G10D"], sigma(G5), aux(G["G10D"]) + own(G["G10D"]) + Jp)  # [(X10, X10, J+); X5]
    J15 = (G["G15X"][0], c2)                                       # J15 : |J15| = a_15, q = C2
    X15 = (G["G15het"][2], c3)                                     # X15 against X15^{0,X}: q = C3
    G["G30_6"] = hetgao(G["G15X"], G["G15het"], aux(G["G15het"]) + J15 + J15)
    G["G25_8"] = hetgao(G["G15AX"], G["G10D"], aux(G["G15Ahet"]) + J15 + J15)
    G["G40_8"] = hetgao(G["G15DX"], G["G25_8"], aux(G["G15Dhet"]) + X15 + J15)
    G["G30L"] = gao(G["G15AX"], sigma(G["G15Ahet"]))
    G["Ghat30"] = hetgao(G["G15AX"], sigma(G["G15Ahet"]), aux(G["G15Ahet"]) + J15 + own(G["G15AX"]))
    G["Ghat25"] = hetgao(G["G15AX"], G["G10D"], aux(G["G15Ahet"]) + J15 + X15)
    a30pp = gao(G["G15het"], G["G15het"])[0]                       # |I30++| = |J30++|
    G["G55"] = hetgao(G["G30L"], sigma(G["Ghat25"]), aux(G["Ghat30"]) + (a30pp, c4) + (a30pp, c4))
    n6, n8, n11 = phi(G["G30_6"]), phi(G["G40_8"]), phi(sigma(G["G55"]))
    n25 = rule(T, (n6, n11, n8))
    return G, n25, apply_words(K, (n25,) * 4), a30pp


def root_digits(M, d, digits):
    """floor(M^(1/d) * 10^digits) by integer bisection"""
    lo, hi = 0, 10 ** (len(str(M)) // d + 2 + digits)
    scale = 10 ** (digits * d)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** d <= M * scale:
            lo = mid
        else:
            hi = mid
    return lo


def tex_check(tex_dir):
    """every constant typed above appears verbatim in the paper's TeX (main.tex and appendix.tex of the arXiv
    source), after removing white space and the LaTeX line-break tokens"""
    import os
    import re
    src = "".join(open(os.path.join(tex_dir, f), encoding="utf-8").read() for f in ("main.tex", "appendix.tex"))
    n = re.sub(r"\\nonumber|\\notag|\\\\|&|\{\}|\\,|\s", "", src)
    tup = lambda v: "(" + ",".join(map(str, v)) + ")"
    found = [tup(v) in n for v in PAPER.values()] + [tup(N25) in n, str(MSTAR) in n, tup(BASE) in n]
    found += [str(c) in n for c in (C1, C2, C3, C4)]
    found += [tup(w) in n for lam in L for w in S3B[lam].split()] + [tup(w) in n for w in K4A]
    return all(found), len(found)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    controls = "--controls" in sys.argv
    ok = True
    tex = [a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--tex=")]
    if tex:
        good, k = tex_check(tex[0])
        print(f"the {k} constants typed in this file (profiles, n_25, M_*, base profile, C1-C4, the 45 S_3b and 49 "
              f"K_4a words) appear verbatim in the paper's TeX: {good}")
        ok &= good
    adm, sepk = admissible(S3B, 3), separated(K4A)
    sym = all((x in SEP[y]) == (y in SEP[x]) for x in L for y in L)
    print(f"S_3b: {sum(len(S3B[x].split()) for x in L)} words, admissible (Definition 4.4 (i) and (ii)): {adm}; "
          f"K_4a: {len(K4A)} words, pairwise separated: {sepk}; the separation relation is symmetric: {sym}")
    ok &= adm and sepk and sym
    G, n25, M, a30pp = chain()
    for k, want in PAPER.items():
        same = G[k] == want
        ok &= same
        print(f"  {k:8s} {G[k]}  {'= the paper' if same else 'DIFFERS FROM THE PAPER ' + str(want)}")
    same = n25 == N25 and a30pp == 2455726444728097
    print(f"  n_25 = S_3b(n_6, n_11, n_8) equals the paper's (Appendix B.4), |I30++| = {a30pp}: {same}")
    ok &= same
    same = M == MSTAR
    r = root_digits(M, 500, 40)
    print(f"Tan2026: M_* = K_4a(n_25) has {len(str(M))} digits and {'equals' if same else 'DIFFERS FROM'} the paper's; "
          f"M_*^(1/500) = 3.{str(r)[1:]}... (40 digits, floor)")
    ok &= same and str(r).startswith("32588326203532663091215390518104754376")
    if args:
        raw = open(args[0], "rb").read()
        h = hashlib.sha256(raw).hexdigest()
        cert = json.loads(raw)
        Mc = int(cert["M"])
        d = int(cert["dimension"])
        print(f"PR #296 certificate.json sha256 {h} ({'as stated' if h == CERT_SHA else 'NOT the stated hash'}); "
              f"M has {len(str(Mc))} digits, d = {d}")
        c1 = Mc ** 500 > MSTAR ** 680
        c2 = Mc * 10 ** (8 * 680) > 325883262 ** 680
        c3 = 3258838705444608 ** 680 <= Mc * 10 ** (15 * 680)
        r = root_digits(Mc, 680, 30)
        print(f"  M^500 > M_*^680: {c1}; M 10^(8*680) > 325883262^680: {c2}; 3.258838705444608^680 <= M: {c3}; "
              f"M^(1/680) = 3.{str(r)[1:]}... (30 digits, floor)")
        ok &= h == CERT_SHA and d == 680 and c1 and c2 and c3
        # our natural reading of the certificate's gadget list (see the docstring): S<L, R, P> = HetGao(L, R) with
        # Tan's codebook triple P; N<A, B> and the outer product = Gao's ordinary product
        CB = {"P15X": ("G15X", aux(G["G15het"]) + (G["G15X"][0], C2) * 1 + (G["G15X"][0], C2)),
              "P15A8": ("G15AX", aux(G["G15Ahet"]) + (G["G15X"][0], C2) + (G["G15X"][0], C2)),
              "P15A5": ("G15AX", aux(G["G15Ahet"]) + (G["G15X"][0], C2) + (G["G15het"][2], C3)),
              "P10A": ("G10A", aux(G["G10A"]) + (G["G10"][0], C1) + own(G["G10A"]))}
        names = {"G5": BASE, "G10": G["G10"], "G15X": G["G15X"], "G15AX": G["G15AX"], "sG15AX": sigma(G["G15AX"]),
                 "sG15Ah": sigma(G["G15Ahet"])}

        def parse(s):
            if s[:2] in ("N<", "S<"):
                body, depth, parts, cur = s[2:-1], 0, [], ""
                for ch in body:
                    if ch == "," and depth == 0:
                        parts.append(cur)
                        cur = ""
                        continue
                    depth += (ch == "<") - (ch == ">")
                    cur += ch
                return (s[0], [parse(p) for p in parts + [cur]])
            return s

        def ev(t):
            if isinstance(t, str):
                return names[t]
            op, a = t
            if op == "S":
                lname, cb = CB[a[2]]
                if ev(a[0]) != G[lname]:
                    raise ValueError("codebook " + a[2] + " belongs to " + lname)
                return hetgao(ev(a[0]), ev(a[1]), cb)
            return gao(ev(a[0]), ev(a[1]))

        env = {}
        try:
            for g in cert["gadgets"]:
                if g["op"] == "sigma":
                    env[g["name"]] = sigma(env[g["of"]])
                    continue
                Rg = env[g["R"]] if g["R"] in env else names[g["R"]]
                pt = parse(g["pattern"])
                if isinstance(pt, str):
                    lname, cb = CB[pt]
                    env[g["name"]] = hetgao(G[lname], Rg, cb)
                else:
                    env[g["name"]] = gao(ev(pt), Rg)
            ins = tuple(phi(env[x[4:-1]]) for x in cert["top"]["rule_inputs"])
            Mn = apply_words(K4A, (rule(S3B, ins),) * 4)
            print(f"  our natural reading of the gadget list gives M' = {str(Mn)[:20]}... ({len(str(Mn))} digits), "
                  f"{'EQUAL to' if Mn == Mc else 'not'} the stated M ({str(Mc)[:20]}...): ratio M'/M = "
                  f"{Decimal(Mn) / Decimal(Mc):.6f}")
        except (KeyError, ValueError) as e:
            print(f"  our natural reading of the gadget list could not be evaluated: {e}")
    if controls:
        print("Red controls (each must be caught):")
        bad = dict(S3B)
        bad["H"] = bad["H"].replace("HVV", "HVA")
        r1 = not admissible(bad, 3)
        print(f"  S_3b with the word (H,V,V) of T_H changed to (H,V,A): not admissible: {r1}")
        r2 = not separated(K4A + ["BBBA"])
        print(f"  K_4a with the word (B,B,B,A) added: not separated: {r2}")
        r3 = chain(c=(C1 + 1, C2, C3, C4))[2] != MSTAR
        r4 = chain(c=(C1, C2, C3, C4 - 1))[2] != MSTAR
        r5 = chain(base=(367, 8, 367, 321, 26, 20))[2] != MSTAR
        print(f"  C1 + 1: M_* changes: {r3}; C4 - 1: M_* changes: {r4}; Gao's base profile (367, 8, 367, 321, 26, 20) "
              f"in place of BPZ's: M_* changes: {r5}")
        ok &= r1 and r2 and r3 and r4 and r5
    print(f"RESULT: {'PASS' if ok else 'FAIL'} (Tan2026 arithmetic and rules; PR #296 comparisons only)")
    return 0 if ok else 1


if __name__ == "__main__":
    getcontext().prec = 30
    sys.exit(main())
