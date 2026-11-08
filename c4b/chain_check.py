"""Replay of the 4b row [Kra2026], C_4b >= alpha_* = 0.752796455874514... (Krachun, arXiv 2608.01325v1, Theorem 1.1),
from the paper's LaTeX source alone. Our own code, standard library only.

1. The ten Paley chains of the table in the proof of Theorem 1.1, read from the source: p prime, p = 3 mod 4, t
   distinct residues, every forward difference s_b - s_a (a < b) a nonzero quadratic residue mod p (Euler's criterion
   and, separately, the set of squares).
2. alpha_* = (10 + sum 1/alpha_i)/(1 + 2 sum 1/alpha_i), alpha_i = log t_i/log p_i, in interval arithmetic of our
   own (Python decimal, 50 digits, outward rounding); the page's form sum log(p_i t_i)/log t_i over
   1 + 2 sum log p_i/log t_i gives the same number.
3. The construction on small instances, by brute force: the local block of Lemma 2.2 lifted by Lemma 2.3, for the
   chain (0, 1) at p = 3 (e = 1, L = 1, 2, 3) and the chains at p = 3 and p = 7 glued by Lemma 2.4 (e = 1, L = 1): the
   sets have the stated sizes and no two elements differ by a nonzero perfect square.
Not checked: the passage to the limit (the multiplicities e_i(U) and the rounding of N), which is analysis, not
computation.

Usage: python3 -I chain_check.py CERTS_ROOT [--controls]
"""
import hashlib
import itertools
import os
import re
import sys
from math import isqrt
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

P = 50
CF = Context(prec=P, rounding=ROUND_FLOOR)
CC = Context(prec=P, rounding=ROUND_CEILING)
CN = Context(prec=P + 5, rounding=ROUND_HALF_EVEN)


def ln_iv(n):
    r = CN.ln(D(n))
    return CF.plus(CN.next_minus(r)), CC.plus(CN.next_plus(r))


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def check_chain(p, S):
    squares = {(x * x) % p for x in range(1, p)}
    ok = is_prime(p) and p % 4 == 3 and len(set(S)) == len(S) and all(0 <= s < p for s in S)
    for a in range(len(S)):
        for b in range(a + 1, len(S)):
            d = (S[b] - S[a]) % p
            ok &= d != 0 and d in squares and pow(d, (p - 1) // 2, p) == 1
    return ok


def alpha_star(pairs):
    """enclosure of (k + sum 1/alpha_i)/(1 + 2 sum 1/alpha_i) with 1/alpha_i = log p_i/log t_i"""
    inv_lo = inv_hi = D(0)
    for p, t in pairs:
        lp, lt = ln_iv(p), ln_iv(t)
        inv_lo = CF.add(inv_lo, CF.divide(lp[0], lt[1]))
        inv_hi = CC.add(inv_hi, CC.divide(lp[1], lt[0]))
    k = len(pairs)
    # f(s) = (k + s)/(1 + 2 s) is decreasing in s for k > 1/2
    lo = CF.divide(CF.add(k, inv_hi), CC.add(1, CC.multiply(2, inv_hi)))
    hi = CC.divide(CC.add(k, inv_lo), CF.add(1, CF.multiply(2, inv_lo)))
    return lo, hi


def local_block(p, S, e):
    """Lemma 2.2: residues mod p^(2e) whose even base-p digits lie in S, with the rank h."""
    t = len(S)
    pos = {s: i for i, s in enumerate(S)}
    out = {}
    for digits in itertools.product(range(p), repeat=2 * e):
        if all(digits[2 * j] in pos for j in range(e)):
            x = sum(d * p ** k for k, d in enumerate(digits))
            out[x] = sum((t - 1 - pos[digits[2 * j]]) * t ** (e - 1 - j) for j in range(e))
    return out, p ** (2 * e), t ** e


def glue(blocks):
    """Lemma 2.4 under CRT: rank sum, H = 1 + sum (H_i - 1)."""
    Ps = [P_ for _, P_, _ in blocks]
    Ptot = 1
    for P_ in Ps:
        Ptot *= P_
    H = 1 + sum(H_ - 1 for _, _, H_ in blocks)
    out = {}
    for combo in itertools.product(*[list(b.items()) for b, _, _ in blocks]):
        x = 0
        for (r, _), P_ in zip(combo, Ps):
            m = Ptot // P_
            x = (x + r * m * pow(m, -1, P_)) % Ptot
        out[x] = sum(h for _, h in combo)
    return out, Ptot, H


def lift(C, P_, H, L):
    """Lemma 2.3: X + P^L h_L(X) over words of length L, translated by one."""
    items = list(C.items())
    A = []
    for word in itertools.product(items, repeat=L):
        X = sum(x * P_ ** j for j, (x, _) in enumerate(word))
        hL = sum(h * H ** (L - 1 - j) for j, (_, h) in enumerate(word))
        A.append(X + P_ ** L * hL + 1)
    return A


def square_difference_free(A):
    A = sorted(A)
    S = set(A)
    top = A[-1] - A[0]
    sq = [k * k for k in range(1, isqrt(top) + 1)]
    return len(S) == len(A) and not any(a + s in S for a in A for s in sq)


# ============================================================================ [Jon2026a]: eleven blocks
def q_set(m):
    return {(z * z) % m for z in range(m)} - {0}


def dag_height(S, m, Q):
    """arcs x -> y when (y - x) mod m in Q; return (acyclic, height = longest path in vertices, arcs)"""
    S = list(S)
    arcs = {x: [y for y in S if y != x and (y - x) % m in Q] for x in S}
    state, depth = {}, {}

    def visit(x):                       # longest path starting at x, counted in vertices; detects cycles
        if state.get(x) == 1:
            raise ValueError("cycle")
        if state.get(x) == 2:
            return depth[x]
        state[x] = 1
        d = 1 + max([visit(y) for y in arcs[x]], default=0)
        state[x] = 2
        depth[x] = d
        return d
    try:
        H = max(visit(x) for x in S)
    except ValueError:
        return False, None, arcs
    return True, H, arcs


def composite_block(m, S, h0, e=1):  # S: the support (the keys of h0)
    """Lemma A of [Jon2026a]: residues mod m^(2e) whose even base-m digits lie in S; h = sum h0(x_2j) H0^(e-1-j)."""
    H0 = max(h0.values()) + 1
    out = {}
    for digits in itertools.product(range(m), repeat=2 * e):
        if all(digits[2 * j] in h0 for j in range(e)):
            x = sum(d * m ** k for k, d in enumerate(digits))
            out[x] = sum(h0[digits[2 * j]] * H0 ** (e - 1 - j) for j in range(e))
    return out, m ** (2 * e), H0 ** e


def main_jon(root, controls):
    def out(x):
        print(x, flush=True)

    out("4b row [Jon2026a]: C_4b >= alpha_inf = 0.753741541837329405... (Jones, fs-lower-bound @ c2c0687, README)")
    status = True
    raws = {}
    for rel in ("jon2026a/README.md", "page/4b.md", "kra2026/main.tex"):
        with open(os.path.join(root, rel), "rb") as fh:
            raws[rel] = fh.read()
        out("%s sha256 %s" % (rel, hashlib.sha256(raws[rel]).hexdigest()))
    readme = raws["jon2026a/README.md"].decode("utf-8")
    page = raws["page/4b.md"].decode("utf-8")
    # block 1 with its exhibited ranking, block 2 support, from the README
    blk1 = re.search(r"Support with its exhibited ranking.*?```\n(.*?)```", readme, re.S).group(1)
    h1 = {int(a): int(b) for a, b in re.findall(r"(\d+):\s*(\d+)", blk1)}
    blk2 = re.search(r"### Composite block 2.*?```\n(.*?)```", readme, re.S).group(1)
    S2 = [int(x) for x in blk2.replace(",", " ").split()]
    # the same sets as printed on the 4b page row
    pg1 = [int(x) for x in re.search(r"on \$235\$ the residues \$\\\\\{([0-9, ]+)\\\\\}\$", page).group(1).split(",")]
    pg2 = [int(x) for x in re.search(r"on \$299\$ the residues \$\\\\\{([0-9, ]+)\\\\\}\$", page).group(1).split(",")]
    same = sorted(h1) == sorted(pg1) and sorted(S2) == sorted(pg2)
    out("block supports read from the README equal those printed on the 4b page row: %s" % same)
    status &= same
    stated = "α∞ = 0.753741541837329405…" in readme and "0.753741541837329405\\ldots" in page
    out("the stated alpha_inf 0.753741541837329405... appears in the README and on the page: %s" % stated)
    status &= stated
    Q235, Q299 = q_set(235), q_set(299)
    out("|Q_235| = %d, |Q_299| = %d (stated 71 and 83): %s" % (len(Q235), len(Q299), (len(Q235), len(Q299)) == (71, 83)))
    status &= (len(Q235), len(Q299)) == (71, 83)
    ok1, H1, arcs1 = dag_height(h1, 235, Q235)
    rank_ok = all(h1[x] > h1[y] for x in h1 for y in arcs1[x])
    out("block 1 (m = 235): %d vertices, square-DAG: %s, height %s, the exhibited ranking drops along every arc: %s "
        "(ranks 0..%d)" % (len(h1), ok1, H1, rank_ok, max(h1.values())))
    ok2, H2, arcs2 = dag_height(S2, 299, Q299)
    out("block 2 (m = 299): %d vertices, square-DAG: %s, height %s" % (len(S2), ok2, H2))
    good = ok1 and ok2 and rank_ok and (len(h1), H1, max(h1.values()) + 1) == (17, 11, 11) and (len(S2), H2) == (19, 12)
    status &= good
    tex = raws["kra2026/main.tex"].decode("utf-8")
    rows = re.findall(r"\$(\d+)\$\s*&\s*\$(\d+)\$\s*&\s*\$\(([0-9,]+)\)\$", tex)
    chains = [(int(p), int(t), tuple(int(x) for x in s.split(","))) for p, t, s in rows if int(p) != 23]
    chok = len(chains) == 9 and all(check_chain(p, S) and len(S) == t for p, t, S in chains)
    out("Krachun's nine chains with p != 23 (from his source): valid Paley chains: %s" % chok)
    status &= chok
    import math
    coprime = all(math.gcd(a, b) == 1 for a, b in itertools.combinations([p for p, _, _ in chains] + [235, 299], 2))
    sqfree = all(235 % (q * q) and 299 % (q * q) for q in range(2, 18))
    out("the eleven moduli are pairwise coprime, and 235 = 5*47, 299 = 13*23 square-free: %s" % (coprime and sqfree))
    status &= coprime and sqfree
    blocks = [(p, t, t) for p, t, _ in chains] + [(235, 17, 11), (299, 19, 12)]
    num_lo = num_hi = den_lo = den_hi = D(0)
    for m, t, H in blocks:
        lmt, lH, lm = ln_iv(m * t), ln_iv(H), ln_iv(m)
        num_lo = CF.add(num_lo, CF.divide(lmt[0], lH[1]))
        num_hi = CC.add(num_hi, CC.divide(lmt[1], lH[0]))
        den_lo = CF.add(den_lo, CF.divide(lm[0], lH[1]))
        den_hi = CC.add(den_hi, CC.divide(lm[1], lH[0]))
    lo = CF.divide(num_lo, CC.add(1, CC.multiply(2, den_hi)))
    hi = CC.divide(num_hi, CF.add(1, CF.multiply(2, den_lo)))
    out("alpha_inf in [%s, %s]" % (str(lo)[:24], str(hi)[:24]))
    tr = D("0.753741541837329405") <= lo and hi < D("0.753741541837329406")
    out("the stated 0.753741541837329405... is its truncation: %s; the Lean-certified 0.7537 < alpha_inf: %s; the page's "
        "'approximately 0.753742' is rounded up (0.753741 would be the truncation): %s"
        % (tr, D("0.7537") < lo, hi < D("0.753742")))
    status &= tr
    # brute force on small instances: Lemma A at e = 1, Lemma B at L = 1, and one glue with Krachun's p = 3 chain
    for m, h0 in ((235, h1), (299, None)):
        if h0 is None:
            ok, H, arcs = dag_height(S2, 299, Q299)
            # reversed longest-path rank: rank(x) = longest path from x (in vertices) - 1
            h0 = {}
            for x in S2:
                h0[x] = dag_height_from(x, arcs)
        C, P_, H_ = composite_block(m, list(h0), h0)
        A = lift(C, P_, H_, 1)
        ok = len(A) == m * len(h0) and square_difference_free(A)
        status &= ok
        out("Lemma A then Lemma B at m = %d (e = 1, L = 1): %d integers in [1, %d], no nonzero square difference: %s"
            % (m, len(A), P_ * H_, ok))
    C3, P3, H3 = local_block(3, (0, 1), 1)
    C235, P235, H235 = composite_block(235, list(h1), h1)
    G, PG, HG = glue([(C3, P3, H3), (C235, P235, H235)])
    A = lift(G, PG, HG, 1)
    ok = len(A) == 6 * 235 * 17 and HG == 12 and square_difference_free(A)
    status &= ok
    out("p = 3 chain glued with the m = 235 block (P = %d, H = %d), L = 1: %d integers, no nonzero square difference: %s"
        % (PG, HG, len(A), ok))
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        # Remark 1 of the README: unit-only arcs are not enough
        m = 15
        Qfull = q_set(m)
        units = {q for q in Qfull if math.gcd(q, m) == 1}
        S15 = {0: 1, 1: 0, 6: 1}
        unit_ok = all(S15[x] > S15[y] for x in S15 for y in S15 if x != y and (y - x) % m in units)
        full_ok = all(S15[x] > S15[y] for x in S15 for y in S15 if x != y and (y - x) % m in Qfull)
        C15, P15, H15 = composite_block(15, list(S15), S15)
        A15 = lift(C15, P15, H15, 1)
        r1 = unit_ok and not full_ok and not square_difference_free(A15)
        out("[ctl] README Remark 1 (m = 15, S = {0, 1, 6}): valid for unit arcs only, refused with all of Q_15, and its "
            "lift has a square difference: %s" % r1)
        ctl &= r1
        h1b = dict(h1)
        h1b[148], h1b[73] = 1, 0
        r2 = not all(h1b[x] > h1b[y] for x in h1b for y in arcs1[x])
        out("[ctl] block 1 with the ranks of 148 and 73 swapped is refused: %s" % r2)
        ctl &= r2
        okb, _, _ = dag_height(S2 + [0], 299, Q299)
        Hb = dag_height(S2 + [0], 299, Q299)[1]
        r3 = (not okb) or Hb != 12
        out("[ctl] block 2 with the residue 0 added: no longer a height-12 square-DAG: %s" % r3)
        ctl &= r3
        lo9 = alpha_star([(p, t) for p, t, _ in chains])[1]
        r4 = lo9 < D("0.7537")
        out("[ctl] without the two composite blocks the nine chains give only %s...: %s" % (str(lo9)[:10], r4))
        ctl &= r4
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT: %s" % ("PASS (two square-DAG blocks of heights 11 and 12, nine chains, alpha_inf = "
                        "0.753741541837329405...; small instances checked by brute force)" if status else "FAIL"))
    return 0 if status else 1


def dag_height_from(x, arcs, memo=None):
    if memo is None:
        memo = {}
    if x in memo:
        return memo[x]
    v = max([dag_height_from(y, arcs, memo) + 1 for y in arcs[x]], default=0)
    memo[x] = v
    return v


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    controls = "--controls" in sys.argv
    if "--row=jon2026a" in sys.argv:
        return main_jon(root, controls)

    def out(x):
        print(x, flush=True)

    out("4b row [Kra2026]: C_4b >= alpha_* = 0.752796455874514... (Krachun, arXiv 2608.01325v1, Theorem 1.1)")
    path = os.path.join(root, "kra2026/main.tex")
    with open(path, "rb") as fh:
        raw = fh.read()
    out("kra2026/main.tex sha256 %s" % hashlib.sha256(raw).hexdigest())
    tex = raw.decode("utf-8")
    rows = re.findall(r"\$(\d+)\$\s*&\s*\$(\d+)\$\s*&\s*\$\(([0-9,]+)\)\$", tex)
    chains = [(int(p), int(t), tuple(int(x) for x in s.split(","))) for p, t, s in rows]
    ours = [(3, 2, (0, 1)), (7, 3, (0, 4, 1)), (11, 4, (0, 3, 1, 4)), (19, 5, (0, 5, 11, 9, 16)),
            (23, 5, (0, 18, 1, 3, 4)), (31, 7, (0, 25, 14, 1, 19, 8, 2)), (43, 7, (0, 31, 9, 23, 4, 40, 1)),
            (59, 9, (0, 49, 15, 7, 16, 19, 35, 36, 5)), (71, 9, (0, 8, 12, 18, 48, 27, 1, 37, 20)),
            (103, 11, (0, 79, 25, 58, 55, 81, 4, 1, 34, 83, 59))]
    status = chains == ours
    out("the table's ten chains, read from the source, equal our transcription: %s" % status)
    thm = "(p_i,t_i)=(3,2),(7,3),(11,4),(19,5),(23,5),\n (31,7),(43,7),(59,9),(71,9),(103,11)," in tex
    stated = "=0.752796455874514\\ldots." in tex
    out("the theorem's ten pairs and the stated 0.752796455874514... appear verbatim: %s" % (thm and stated))
    status &= thm and stated
    for p, t, S in chains:
        ok = check_chain(p, S) and len(S) == t
        status &= ok
        out("p = %3d, t = %2d, chain %s: prime, 3 mod 4, distinct, every forward difference a nonzero square: %s"
            % (p, t, S, ok))
    lo, hi = alpha_star([(p, t) for p, t, _ in chains])
    out("alpha_* in [%s, %s]" % (str(lo)[:22], str(hi)[:22]))
    trunc = D("0.752796455874514") <= lo and hi < D("0.752796455874515")
    page = D("0.752796") <= lo
    out("the paper's 0.752796455874514... is its truncation: %s; the page's 0.752796 is below it: %s" % (trunc, page))
    status &= trunc and page
    # small instances of the construction, by brute force
    C3, P3, H3 = local_block(3, (0, 1), 1)
    for L in (1, 2, 3):
        A = lift(C3, P3, H3, L)
        ok = len(A) == 6 ** L and max(A) <= (P3 * H3) ** L and square_difference_free(A)
        status &= ok
        out("p = 3, S = (0, 1), e = 1, L = %d: %d integers in [1, %d], no nonzero square difference: %s"
            % (L, len(A), (P3 * H3) ** L, ok))
    C7, P7, H7 = local_block(7, (0, 4, 1), 1)
    G, PG, HG = glue([(C3, P3, H3), (C7, P7, H7)])
    A = lift(G, PG, HG, 1)
    ok = len(A) == 6 * 21 and HG == 4 and max(A) <= PG * HG and square_difference_free(A)
    status &= ok
    out("p = 3 and p = 7 glued (P = 441, H = 4), L = 1: %d integers in [1, %d], no nonzero square difference: %s"
        % (len(A), PG * HG, ok))
    C9, P9, H9 = local_block(3, (0, 1), 2)
    A = lift(C9, P9, H9, 1)
    ok = len(A) == 36 and H9 == 4 and square_difference_free(A)
    status &= ok
    out("p = 3, S = (0, 1), e = 2 (P = 81, H = 4), L = 1: %d integers, no nonzero square difference: %s" % (len(A), ok))
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        r1 = not check_chain(7, (0, 1, 4))
        out("[ctl] the p = 7 chain reordered as (0, 1, 4) is refused (1 - 4 = 4 is a square, 4 - 1 = 3 is not): %s"
            % r1)
        ctl &= r1
        r2 = not check_chain(13, (0, 1, 4))
        out("[ctl] p = 13 (1 mod 4) is refused: %s" % r2)
        ctl &= r2
        r3 = not check_chain(103, (0, 79, 25, 58, 55, 81, 4, 1, 34, 83, 60))
        out("[ctl] the p = 103 chain with its last entry 59 -> 60 is refused: %s" % r3)
        ctl &= r3
        Cb, Pb, Hb = local_block(3, (0, 1), 1)
        Cb = {x: 0 for x in Cb}
        r4 = not square_difference_free(lift(Cb, Pb, 1, 2))
        out("[ctl] the lift with all ranks set to 0 (H = 1) has a square difference: %s" % r4)
        ctl &= r4
        lo2, hi2 = alpha_star([(p, t) for p, t, _ in chains if p != 103])
        r5 = hi2 < D("0.752796")
        out("[ctl] without the p = 103 chain alpha_* drops to %s...: %s" % (str(hi2)[:12], r5))
        ctl &= r5
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT: %s" % ("PASS (ten valid Paley chains; alpha_* = 0.752796455874514...; the construction checked by "
                        "brute force on small instances)" if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
