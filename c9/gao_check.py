"""Replay of the 9a row [Gao2026], Theta(C_7) >= M_40^(1/200) = 3.258789153908... (Gao, arXiv 2607.27869v1), from
the paper's LaTeX source and the 367-point base set alone (the repository's verification script was neither read nor
run). Our own code, standard library only.

1. The base gadget (Proposition 5 of the paper; Definition 3): I_0 (367 points of Z_7^5) is independent in
   C_7^5; the eight pairs (r_j, q_j) of the paper's table are private pairs (r_j in I_0, q_j not, and the closed
   neighbourhood of q_j meets I_0 only in r_j) with 16 distinct endpoints; the transversals PH, PV (J_0 = {0,5,6})
   are independent; X = T(I_0) with one point swapped is independent, of size 367, disjoint from the endpoints, and
   splits into (o, h, v) = (321, 26, 20) with no point confusable with both transversals.
2. The product lemma (Lemma 4) at k = 2, as a test of the construction: the code I_12 built from two copies is
   written out (134,753 points of Z_7^10) for strong_check.c, which must find it independent.
3. The recursion: s, o, t computed for every k <= 40 from the product formulas over every split (they must not
   depend on the split, and must equal the closed forms of Lemma 7); M_k by the maximum over splits (the recursion
   before Proposition 8); the paper's split tree attains it; M_40 is the stated integer; M_40^(1/200) to 40 digits.

Usage: python3 -I gao_check.py CERTS_ROOT OUTFILE [--controls]   (OUTFILE receives the k = 2 code, one point a line)
"""
import hashlib
import itertools
import os
import sys
from decimal import Decimal as D, Context, ROUND_HALF_EVEN

TEX = [
    r"0&(1,3,4,4,6)&(2,3,5,4,6)\\",
    r"1&(3,4,0,3,5)&(2,4,6,3,5)\\",
    r"2&(5,3,1,3,4)&(5,3,2,3,5)\\",
    r"3&(4,4,6,1,6)&(5,4,6,0,6)\\",
    r"4&(6,0,6,4,5)&(6,1,6,5,5)\\",
    r"5&(0,3,5,6,5)&(6,3,5,0,5)\\",
    r"6&(6,4,3,4,0)&(6,4,2,4,6)\\",
    r"7&(6,4,5,3,2)&(6,5,5,3,1)",
    r"J_0=\{0,5,6\},\qquad J_1=\{1,2,3,4,7\},",
    r"T(w)=(2-w_1,w_3,w_0,2-w_2,w_4).",
    r"X=\bigl(X_0\setminus\{(2,4,6,3,5)\}\bigr) \cup\{(1,5,6,3,5)\}.",
    r"(a,t,s,o,h,v)=(367,8,367,321,26,20).",
    r"1+1=2,\quad 1+2=3,\quad 2+3=5,\quad 5+5=10,\quad 10+10=20,\quad 20+20=40.",
    r"=3.2587891539086910161967650155206769\ldots .",
]
PAIRS = [((1, 3, 4, 4, 6), (2, 3, 5, 4, 6)), ((3, 4, 0, 3, 5), (2, 4, 6, 3, 5)), ((5, 3, 1, 3, 4), (5, 3, 2, 3, 5)),
         ((4, 4, 6, 1, 6), (5, 4, 6, 0, 6)), ((6, 0, 6, 4, 5), (6, 1, 6, 5, 5)), ((0, 3, 5, 6, 5), (6, 3, 5, 0, 5)),
         ((6, 4, 3, 4, 0), (6, 4, 2, 4, 6)), ((6, 4, 5, 3, 2), (6, 5, 5, 3, 1))]
J0 = {0, 5, 6}
M40 = int("4085567963379990282748597333793399773399910167372462112610203714626938509773125898252337545387866277249")
STATED = "3.2587891539086910161967650155206769"
SPLITS = [(1, 1), (1, 2), (2, 3), (5, 5), (10, 10), (20, 20)]


def conf(x, y):
    """closed adjacency in a strong power of C_7: every coordinate differs by 0 or +-1 mod 7"""
    return all((a - b) % 7 in (0, 1, 6) for a, b in zip(x, y))


def independent(S):
    S = list(S)
    if len(set(S)) != len(S):
        return False
    return not any(conf(S[i], S[j]) for i in range(len(S)) for j in range(i + 1, len(S)))


def T(w):
    return ((2 - w[1]) % 7, w[3], w[0], (2 - w[2]) % 7, w[4])


def gadget(I0, pairs, j0, swap=True):
    """Return (ok, details) for the base gadget built from the data."""
    I0s = set(I0)
    R = [r for r, q in pairs]
    Qs = [q for r, q in pairs]
    ends = R + Qs
    PH = [r if j in j0 else q for j, (r, q) in enumerate(pairs)]
    PV = [q if j in j0 else r for j, (r, q) in enumerate(pairs)]
    X0 = [T(w) for w in I0]
    if swap:
        X = [x for x in X0 if x != (2, 4, 6, 3, 5)] + [(1, 5, 6, 3, 5)]
    else:
        X = X0
    d = {}
    d["I0 independent, 367 points"] = len(I0s) == len(I0) == 367 and independent(I0)
    d["private pairs: r_j in I0, q_j not, N(q_j) meets I0 only in r_j"] = all(
        r in I0s and q not in I0s and [p for p in I0 if conf(p, q)] == [r] for r, q in pairs)
    d["16 distinct endpoints"] = len(set(ends)) == 16
    d["PH and PV independent"] = independent(PH) and independent(PV)
    d["X independent, 367 points"] = len(set(X)) == len(X) == 367 and independent(X)
    d["X disjoint from the endpoints"] = not set(X) & set(ends)
    inH = [any(conf(x, p) for p in PH) for x in X]
    inV = [any(conf(x, p) for p in PV) for x in X]
    both = sum(1 for a, b in zip(inH, inV) if a and b)
    o = sum(1 for a, b in zip(inH, inV) if not a and not b)
    h = sum(1 for a, b in zip(inH, inV) if a and not b)
    v = sum(1 for a, b in zip(inH, inV) if b and not a)
    d["no point of X confusable with both transversals"] = both == 0
    d["(o, h, v) = (321, 26, 20)"] = (o, h, v) == (321, 26, 20)
    ok = all(d.values())
    return ok, d, dict(I=I0, R=R, Q=Qs, PH=PH, PV=PV, X=X, inH=inH, inV=inV)


def product_code(g1, g2):
    """I_12 = (B1 x B2) u A_H u A_V of Lemma 4."""
    B1 = [p for p in g1["I"] if p not in set(g1["R"])]
    B2 = [p for p in g2["I"] if p not in set(g2["R"])]
    code = [b1 + b2 for b1 in B1 for b2 in B2]
    for x, a, b in zip(g2["X"], g2["inH"], g2["inV"]):          # A_H: (p, x), p in L_1(x)
        L = g1["R"] if not a and not b else (g1["PH"] if a else g1["PV"])
        code += [p + x for p in L]
    for y, a, b in zip(g1["X"], g1["inH"], g1["inV"]):          # A_V: (y, q), q in K_2(y)
        K = g2["R"] if not a and not b else (g2["PV"] if a else g2["PH"])
        code += [y + q for q in K]
    return code


def recursion(out):
    base = (8, 367, 321)                                         # (t, s, o)
    tso = {1: base}
    consistent = True
    for k in range(2, 41):
        vals = set()
        for i in range(1, k):
            (t1, s1, o1), (t2, s2, o2) = tso[i], tso[k - i]
            vals.add((t1 * o2 + o1 * t2, s1 * s2, o1 * o2 + (s1 - o1) * (s2 - o2)))
        consistent &= len(vals) == 1
        tso[k] = vals.pop()
    closed = all(tso[k] == (2 * (367 ** k - 275 ** k) // 23, 367 ** k, (367 ** k + 275 ** k) // 2)
                 and (2 * (367 ** k - 275 ** k)) % 23 == 0 for k in range(1, 41))
    M = {1: 367}
    best_split = {}
    for k in range(2, 41):
        cands = {}
        for i in range(1, k):
            j = k - i
            ti, si, _ = tso[i]
            tj, sj, _ = tso[j]
            cands[i] = (M[i] - ti) * (M[j] - tj) + ti * sj + si * tj
        M[k] = max(cands.values())
        best_split[k] = cands
    tree = all(best_split[i + j][i] == M[i + j] for i, j in SPLITS)
    out("t, s, o do not depend on the split for any k <= 40: %s; they equal the closed forms of Lemma 7: %s"
        % (consistent, closed))
    out("M_k = max over splits for k <= 40; the split tree 1+1, 1+2, 2+3, 5+5, 10+10, 20+20 attains it: %s" % tree)
    out("M_2 = %d (= 359^2 + 2*8*367: %s); M_40 has %d digits; equal to the stated integer: %s"
        % (M[2], M[2] == 359 ** 2 + 2 * 8 * 367 == 134753, len(str(M[40])), M[40] == M40))
    return consistent and closed and tree and M[40] == M40, M


def root200(n, prec=70):
    c = Context(prec=prec, rounding=ROUND_HALF_EVEN)
    v = c.exp(c.divide(c.ln(D(n)), D(200)))
    return v


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    root, outfile = sys.argv[1], sys.argv[2]
    controls = "--controls" in sys.argv

    def out(x):
        print(x, flush=True)

    out("9a row [Gao2026]: Theta(C_7) >= M_40^(1/200) = 3.258789153908... (arXiv 2607.27869v1, Theorem 1)")
    status = True
    tp = os.path.join(root, "gao2026/paper.tex")
    with open(tp, "rb") as fh:
        raw = fh.read()
    out("gao2026/paper.tex sha256 %s" % hashlib.sha256(raw).hexdigest())
    text = " ".join(raw.decode("utf-8").split())
    for t in TEX:
        hit = " ".join(t.split()) in text
        status &= hit
        out("  our transcription appears verbatim: %s   %s" % (hit, t[:90]))
    rp = os.path.join(root, "gao2026/R367.txt")
    with open(rp, "rb") as fh:
        rraw = fh.read()
    out("gao2026/R367.txt sha256 %s" % hashlib.sha256(rraw).hexdigest())
    I0 = []
    for line in rraw.decode("ascii").strip().split("\n"):
        parts = line.split()
        if len(parts) != 5 or not all(p.isdigit() and 0 <= int(p) <= 6 for p in parts):
            out("malformed line: %r" % line)
            return 1
        I0.append(tuple(int(p) for p in parts))
    ok, d, g = gadget(I0, PAIRS, J0)
    for k, v in d.items():
        out("base gadget: %s: %s" % (k, v))
    status &= ok
    code = product_code(g, g)
    out("k = 2 product code (Lemma 4 on two copies): %d points, = (367-8)^2 + 2*8*367 = 134753: %s"
        % (len(code), len(code) == 134753))
    status &= len(code) == 134753 and len(set(code)) == len(code)
    with open(outfile, "w") as fh:
        for p in code:
            fh.write(" ".join(map(str, p)) + "\n")
    out("k = 2 product code written for strong_check.c (independence is checked there)")
    rok, M = recursion(out)
    status &= rok
    v = root200(M[40])
    stated_ok = D(STATED) <= v < D(STATED) + D("1e-34")
    page_ok = D("3.258789153908") <= v
    out("M_40^(1/200) = %s... (70-digit decimal ln and exp; the stated 34 decimals are a truncation: %s; the page's "
        "3.258789153908 is below it: %s)" % (str(v)[:45], stated_ok, page_ok))
    status &= stated_ok and page_ok
    # the stated value exactly: 3.258789153908^200 <= M_40 < 3.258789153909^200
    exact = 3258789153908 ** 200 <= M40 * 10 ** (12 * 200) < 3258789153909 ** 200
    out("exactly: 3.258789153908^200 <= M_40 < 3.258789153909^200: %s" % exact)
    status &= exact
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        ok2, d2, _ = gadget(I0, PAIRS, J0, swap=False)
        r1 = not ok2 and not d2["X disjoint from the endpoints"]
        out("[ctl] X = T(I_0) without the swap contains the endpoint (2,4,6,3,5): refused: %s" % r1)
        ctl &= r1
        bad_pairs = list(PAIRS)
        bad_pairs[2] = (PAIRS[2][0], (5, 3, 3, 3, 5))
        ok3, d3, _ = gadget(I0, bad_pairs, J0)
        r2 = not d3["private pairs: r_j in I0, q_j not, N(q_j) meets I0 only in r_j"]
        out("[ctl] q_2 = (5,3,3,3,5) in place of (5,3,2,3,5): no longer a private pair: %s" % r2)
        ctl &= r2
        moved = list(I0)
        moved[0] = ((moved[0][0] + 1) % 7,) + moved[0][1:]
        r3 = not independent(moved)
        out("[ctl] the first point of I_0 moved by +1 in its first coordinate: I_0 no longer independent: %s" % r3)
        ctl &= r3
        ok4, d4, _ = gadget(I0, PAIRS, {0, 5})
        r4 = not ok4
        out("[ctl] J_0 = {0, 5} (pair 6 moved to the other transversal): refused (%s): %s"
            % ([k for k, x in d4.items() if not x], r4))
        ctl &= r4
        r5 = not (3258789153909 ** 200 <= M40 * 10 ** (12 * 200))
        out("[ctl] the value 3.258789153909 is refused: %s" % r5)
        ctl &= r5
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT (Python part): %s" % ("PASS" if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
