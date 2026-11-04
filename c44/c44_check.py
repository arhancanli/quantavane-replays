#!/usr/bin/env python3
"""Independent replay of the asterisked upper bound C_44 <= 4.3576549226* ([P2026], R. Pinho; PR #258 of
teorth/optimizationproblems, constant 44a) from the certificate's data files.

Written from the 44a page, the paper (paper/c44.tex at f3fb25d), Wellens's paper (Discrete Analysis 2022,
arXiv 2005.00566v2, Section 3.3: Lemmas 3.8, 3.9, Corollary 3.10, Fact 3.11, Proposition 3.12) and the
certificate's README. The claimant's programs (verify_c44_exact.py, independent/verify_c44_fraction.py,
reproduce_c44.py, reproduce.sh) were neither read nor run. Standard library only (Python 3.9+); run with
python3 -I. Every decision is an exact integer or rational comparison.

  c44_check.py CERTS_ROOT [--controls]

CERTS_ROOT holds pin/ (rpinho/c44-certificate @f3fb25d: certificates/c44_candidate.json,
sources/wellens-certificate.json, sources/wellens-result.json) and epi/ (BakeLab/Epi-results @97fa6a4:
mathematics/wellens/certificate.json and result.json), as fetch.sh downloads them.

What is computed (the definitions are the paper's and Wellens's)
  D(f) = sum over relevant i of 2^-deg_i(f), deg_i(f) = deg(f(x) - f(x^i)) = largest degree of a monomial
  containing x_i; r(f) <= D(f) 2^deg(f). The bound is U(B_H, H) + T_H at the horizon H = 96, where
    U(0, d) = 0,   U(b, d) = min(c_d, d 2^-d + max_{1<=k<=d} U(min(b - 1, B_k), k))   for 1 <= b <= B_d,
    c_d = d/2 (the influence bound) for d != 3 and c_3 = 5/4 (the paper's Lemma 2);
    B_d = x_d - 1 from Epi's exact certificate (nodes 1 = x_0 < ... < x_d, weights w_i with
    sum_i w_i x_i^j = 0 for j = 1..d, w_0 > 0, w_0 + sum_{i>=1, w_i<0} w_i > 0) for d <= 48, intersected with
    Wellens's 3b(b-1) <= 2(d^4 - d^2); for 49 <= d <= 96 the largest b with 3b(b-1) <= 2(d^4 - d^2);
    T_H = (H+1)^3 / 2^(H+1) + sum_{k>=H+2} (2k^2 - k)/2^k   (Wellens (3.15)).
  The recurrence is computed three ways: (A) row by row, exact integers at scale 2^(H+1), prefix maxima;
  (B) column by column, exact integers at scale 2^(2H+3), the max over k taken explicitly; (C) the literal
  memoised recursion in Fractions for d <= 30. A and B must agree on every U(B_d, d), C on d <= 30. T_H is
  summed term by term to k = 4000 and the rest closed by a telescoping identity checked coefficient by
  coefficient, and compared with the paper's closed form.
  The finite parts of the cap lemma are enumerated: every {-1,0,1}-valued affine function of s <= 3 sign
  variables (all 3^(2^s) truth tables, Walsh coefficients), and the ten-variable witness F with D(F) = 5/4
  (Walsh-Hadamard transform of all 1024 values).
  Third-party anchors for the recurrence code: with Wellens's caps d/2 and his Table 1 it reproduces his
  D(30^2, 30) = 4.4157... and D_inf <= 4.3935, and the 3.9375 at degree 8 quoted on the page from [Wel2019];
  with Epi's data it reproduces the paper's preceding exact value and the README's horizon-72 value.
  Controls (--controls) include the bound as an exact function of the degree-3 cap c: it is constant for
  c <= 173/128, so the claim needs only D <= 173/128 at degree 3 (still R_3 <= 10: 11/8 is too much).
  After our run, the claimant's published outputs (verification/c44_exact.log: margin, tail) were compared.

What rests on mathematics rather than computation (read, not computed): r(f) <= D(f) 2^deg f; the
restriction inequality (Wellens (3.2), Fact 3.3); Lemma 3.8 (fixing a top-degree monomial lowers block
sensitivity) and Lemma 3.9 (the recurrence); Corollary 3.10 (3.14)-(3.15) (the tail, which needs
bs <= deg^2 above the horizon); the influence bound D <= I/2 <= deg/2; Wellens's Section 4 inequality
B_d^2 - B_d <= (2/3)(d^4 - d^2); the symmetrization (Fact 3.11, Prop. 3.12) that turns Epi's identities into
bs <= x_d - 1; the paper's Lemma 2 (D <= 5/4 at degree <= 3), whose prose proof uses the published R_3 <= 10
(Tarannikov-Kirienko 2000; Krotov-Valyuzhenich 2024). Only the lemma's finite cases are computed here.
"""
import hashlib
import json
import math
import os
import re
import sys
import time
from fractions import Fraction
from itertools import product

PIN = "pin"
EPI = "epi"
F_CAND = "pin/certificates/c44_candidate.json"
F_SRC = "pin/sources/wellens-certificate.json"
F_RES = "pin/sources/wellens-result.json"
F_EPI_CERT = "epi/certificate.json"
F_EPI_RES = "epi/result.json"

# Hashes stated by others: rpinho/c44-certificate SHA256SUMS @f3fb25d, the paper (Section 4), the candidate's
# source_sha256 field, THIRD_PARTY_NOTICES.md, BakeLab's checksums.sha256 and result.json @97fa6a4.
H_CAND = "95bd09f3c2684abd46f663ce4b444ae482119a56722384f8c797939a0beb823c"
H_SRC = "ce66f1097a7284e9845113fa1dc39267932cc934fd23d9019ff13d3dcb39cf49"
H_RES = "554994f8a2f7bce361d72f8f9100890d19d7f7343f56d59a6aa61e0f76627df3"

# Values printed by others, compared with ours at the end
PAGE_CLAIM = "4.3576549226"
PAGE_EXACT = Fraction(172624496193406565421448533795, 39614081257132168796771975168)     # 44a page, H = 96
PAPER_EPI_EXACT = Fraction(172933981203227910490173314851, 39614081257132168796771975168)  # paper, cap 3/2
README_H72 = Fraction(10289221775138769568521, 2361183241434822606848)                   # README, H = 72
EPI_PRINTED = "4.3654674227"
# the claimant's outputs, verification/c44_exact.log @f3fb25d (read after our run)
LOG_MARGIN = Fraction(9881116158388503765235957, 386856262276681335905976320000000000)
LOG_TAIL = Fraction(233045, 39614081257132168796771975168)
WELLENS_TABLE1 = [1, 3, 6, 10, 15, 21, 29, 38, 47, 58, 71, 84, 99, 114]    # Wellens 2022, Table 1, d = 1..14

RAT = re.compile(r"-?[0-9]+(/[0-9]+)?$")


# ----------------------------------------------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------------------------------------------

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rat(v):
    if isinstance(v, bool):
        raise ValueError("boolean where a rational is expected")
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str) and RAT.match(v):
        return Fraction(v)
    raise ValueError("not an exact rational: %r" % (v,))


def dec(fr, digits, up=False):
    """fr rounded down (or up) to `digits` decimals, as a string."""
    q, r = divmod(fr.numerator * 10 ** digits, fr.denominator)
    if up and r:
        q += 1
    s = str(abs(q)).rjust(digits + 1, "0")
    return ("-" if q < 0 else "") + s[:-digits] + "." + s[-digits:]


def pow2_str(fr):
    """'a / 2^e' when the denominator is a power of two."""
    d = fr.denominator
    if d & (d - 1) == 0:
        return "%d / 2^%d" % (fr.numerator, d.bit_length() - 1)
    return "%d / %d" % (fr.numerator, d)


class Report:
    def __init__(self, title):
        self.fails = []
        self.t0 = time.time()
        print(title, flush=True)

    def info(self, key, val):
        print("  %s: %s" % (key, val), flush=True)

    def check(self, desc, ok, detail=""):
        print("  check %s: %s%s" % (desc, "PASS" if ok else "FAIL", (" (" + detail + ")") if detail else ""),
              flush=True)
        if not ok:
            self.fails.append(desc)
        return ok

    def done(self):
        print("  run time: %.1f s" % (time.time() - self.t0))
        print("  RESULT: %s" % ("ALL CHECKS PASS" if not self.fails else "FAILED: " + "; ".join(self.fails)))
        return not self.fails


# ----------------------------------------------------------------------------------------------------------------
# Block-sensitivity inputs
# ----------------------------------------------------------------------------------------------------------------

def parse_epi(src):
    """{d: (nodes, weights)} from Epi's certificate.json, read as data."""
    out = {}
    for key, item in src["certificates"].items():
        d = int(key)
        if str(d) != key:
            raise ValueError("bad degree key %r" % key)
        nodes = item["nodes"]
        if not all(isinstance(x, int) and not isinstance(x, bool) for x in nodes):
            raise ValueError("degree %d: nodes are not integers" % d)
        out[d] = (list(nodes), [rat(w) for w in item["weights"]])
    return out


def epi_condition(d, nodes, weights):
    """The exact conditions of the paper, Section 3, for one degree. Returns (ok, reason)."""
    if len(nodes) != d + 1 or len(weights) != d + 1:
        return False, "%d nodes and %d weights, expected %d" % (len(nodes), len(weights), d + 1)
    if nodes[0] != 1 or any(a >= b for a, b in zip(nodes, nodes[1:])):
        return False, "nodes not 1 = x_0 < x_1 < ... < x_d"
    if weights[0] <= 0:
        return False, "w_0 <= 0"
    if weights[0] + sum(w for w in weights[1:] if w < 0) <= 0:
        return False, "w_0 + sum of negative w_i <= 0"
    L = 1
    for w in weights:
        L = L * w.denominator // math.gcd(L, w.denominator)
    W = [int(w * L) for w in weights]                         # exact integers
    for j in range(1, d + 1):
        if sum(Wi * x ** j for Wi, x in zip(W, nodes)) != 0:
            return False, "moment identity fails at j = %d" % j
    return True, "ok"


def lagrange_condition(d, nodes, weights):
    """Second, different identity: w_i x_i prod_{j != i} (x_i - x_j) is one nonzero constant (with d + 1 distinct
    nodes this is equivalent to sum_i w_i x_i^j = 0 for j = 1..d; checked separately as a cross-check)."""
    vals = set()
    for i, (x, w) in enumerate(zip(nodes, weights)):
        p = 1
        for j, y in enumerate(nodes):
            if j != i:
                p *= (x - y)
        vals.add(w * x * p)
    return len(vals) == 1 and 0 not in vals


def wellens_b(d):
    """Largest integer b >= 0 with 3 b (b - 1) <= 2 (d^4 - d^2) (Wellens 2022, Section 4)."""
    rhs = 2 * (d ** 4 - d ** 2)
    b = (1 + math.isqrt(1 + 4 * rhs // 3)) // 2
    while 3 * (b + 1) * b <= rhs:
        b += 1
    while b > 0 and 3 * b * (b - 1) > rhs:
        b -= 1
    assert 3 * b * (b - 1) <= rhs < 3 * (b + 1) * b
    return b


def bs_bounds(H, epi_ok, cutoff, above="wellens"):
    """B[0..H]: B[d] = min(x_d - 1, W_d) where Epi's degree-d certificate passed (d <= cutoff), else W_d;
    above='square' replaces W_d by d^2 for d > cutoff (a control)."""
    B = [0] * (H + 1)
    src = [""] * (H + 1)
    for d in range(1, H + 1):
        W = wellens_b(d) if (above == "wellens" or d <= cutoff) else d * d
        if d <= cutoff and d in epi_ok:
            e = epi_ok[d]
            B[d] = min(e, W)
            src[d] = "epi" if e < W else ("both" if e == W else "wellens")
        else:
            B[d] = W
            src[d] = "wellens" if above == "wellens" else "square"
    return B, src


# ----------------------------------------------------------------------------------------------------------------
# The recurrence, three ways
# ----------------------------------------------------------------------------------------------------------------

def caps_list(H, cap3):
    c = [Fraction(0)] + [Fraction(d, 2) for d in range(1, H + 1)]
    if H >= 3:
        c[3] = Fraction(cap3)
    return c


def rec_rows(B, H, caps):
    """(A) Row by row. Values are integers at scale 2^E, E = H + 1. Returns (E, Uf) with Uf[d] = 2^E U(B_d, d)."""
    E = H + 1
    w = [0] + [d << (E - d) for d in range(1, H + 1)]
    cap = [0] * (H + 1)
    for d in range(1, H + 1):
        v = caps[d] * (1 << E)
        if v.denominator != 1:
            raise ValueError("cap c_%d = %s is not a multiple of 2^-%d" % (d, caps[d], E))
        cap[d] = int(v)
    prev = [0] * (H + 1)                 # prev[d] = U(min(b - 1, B_d), d)
    for b in range(1, max(B[1:]) + 1):
        new = list(prev)
        pm = 0                           # running max over k <= d of prev[k]
        for d in range(1, H + 1):
            if prev[d] > pm:
                pm = prev[d]
            if b <= B[d]:
                v = w[d] + pm
                new[d] = v if v < cap[d] else cap[d]
        prev = new
    return E, prev


def rec_cols(B, H, caps):
    """(B) Column by column, degree d = 1..H in turn, the max over k = 1..d taken explicitly for every b.
    Integers at scale 2^S, S = 2H + 3. Returns (S, Uf) with Uf[d] = 2^S U(B_d, d)."""
    S = 2 * H + 3
    col = [None] * (H + 1)
    out = [0] * (H + 1)
    for d in range(1, H + 1):
        wd = d * (1 << (S - d))
        cd = caps[d] * (1 << S)
        assert cd.denominator == 1
        cd = int(cd)
        c = [0]
        prevcols = [(col[k], B[k]) for k in range(1, d)]
        for b in range(1, B[d] + 1):
            m = c[b - 1]                                    # k = d
            for ck, Bk in prevcols:                         # k = 1..d-1
                v = ck[b - 1] if b - 1 <= Bk else ck[Bk]
                if v > m:
                    m = v
            v = wd + m
            c.append(v if v < cd else cd)
        col[d] = c
        out[d] = c[B[d]]
    return S, out


def rec_literal(B, Dmax, caps):
    """(C) The definition as written, memoised, in Fractions: U(b, d) for d <= Dmax."""
    sys.setrecursionlimit(100000)
    memo = {}

    def U(b, d):
        if b > B[d]:
            b = B[d]
        if b == 0:
            return Fraction(0)
        key = (b, d)
        if key in memo:
            return memo[key]
        # fill the column bottom-up to keep the recursion shallow
        for bb in range(1, b + 1):
            if (bb, d) in memo:
                continue
            m = max(U(bb - 1, k) for k in range(1, d + 1))
            memo[(bb, d)] = min(caps[d], Fraction(d, 2 ** d) + m)
        return memo[key]

    return [Fraction(0)] + [U(B[d], d) for d in range(1, Dmax + 1)]


def tail_closed(H):
    return Fraction((H + 1) ** 3 + 2 * H * H + 11 * H + 19, 2 ** (H + 1))


def tail_series(H, K=4000):
    """(H+1)^3/2^(H+1) + sum_{k=H+2}^{K} (2k^2 - k)/2^k + S(K+1), with S(m) = (2m^2 + 3m + 5)/2^(m-1) the exact
    value of sum_{k>=m} (2k^2 - k)/2^k (telescoping identity checked by tail_identity)."""
    s = Fraction((H + 1) ** 3, 2 ** (H + 1))
    num = 0
    for k in range(H + 2, K + 1):          # common denominator 2^K
        num += (2 * k * k - k) << (K - k)
    s += Fraction(num, 2 ** K)
    m = K + 1
    return s + Fraction(2 * m * m + 3 * m + 5, 2 ** (m - 1)), Fraction(num, 2 ** K)


def tail_identity():
    """S(m) - S(m+1) = (2m^2 - m)/2^m as polynomials: 2 (2m^2+3m+5) - (2(m+1)^2 + 3(m+1) + 5) = 2m^2 - m.
    Coefficients (constant first) compared exactly; S(m) -> 0 since S(m) <= 10 m^2 / 2^(m-1)."""
    p = [10, 6, 4]                                   # 2 (2m^2 + 3m + 5)
    q = [2 + 3 + 5, 4 + 3, 2]                        # 2(m+1)^2 + 3(m+1) + 5 = 2m^2 + 7m + 10
    diff = [a - b for a, b in zip(p, q)]
    return diff == [0, -1, 2]


def bound_at(H, epi_ok, cutoff, cap3, above="wellens", B_override=None):
    """U(B_H, H) + T_H with recurrence A; returns (bound, Uf as Fractions, B)."""
    B, _ = bs_bounds(H, epi_ok, cutoff, above)
    if B_override:
        for d, v in B_override.items():
            B[d] = v
    E, Uf = rec_rows(B, H, caps_list(H, cap3))
    U = [Fraction(v, 1 << E) for v in Uf]
    return U[H] + tail_closed(H), U, B


# ----------------------------------------------------------------------------------------------------------------
# The cap lemma: finite parts
# ----------------------------------------------------------------------------------------------------------------

def affine_forms(smax=3):
    """All A: {-1,1}^s -> {-1,0,1}, A != 0, that are affine (Walsh coefficients of degree >= 2 vanish) and depend on
    every one of the s variables. Returns {s: [(table, p)]}."""
    found = {}
    for s in range(smax + 1):
        pts = list(product((-1, 1), repeat=s))
        masks = list(range(1 << s))
        chi = [[1 if bin(S & sum(1 << j for j in range(s) if pt[j] < 0)).count("1") % 2 == 0 else -1
                for pt in pts] for S in masks]
        res = []
        for vals in product((-1, 0, 1), repeat=len(pts)):
            if not any(vals):
                continue
            c = [sum(v * x for v, x in zip(vals, chi[S])) for S in masks]
            if any(c[S] for S in masks if bin(S).count("1") >= 2):
                continue
            if any(c[1 << j] == 0 for j in range(s)):
                continue
            p = Fraction(sum(1 for v in vals if v), len(pts))
            res.append((vals, p))
        found[s] = res
    return found


def witness_F():
    """F(s,t,x,y) = ((s+t) Q(x) + (s-t) Q(y))/2, Q(s,t,u,v) = ((s+t) u + (s-t) v)/2, on {-1,1}^10.
    Variables 0..9 = s, t, x1..x4, y1..y4. Returns (values, Walsh coefficients c_S = sum_x F(x) chi_S(x))."""
    def Q(a, b, u, v):
        return Fraction((a + b) * u + (a - b) * v, 2)
    vals = []
    for m in range(1 << 10):
        z = [(-1 if (m >> i) & 1 else 1) for i in range(10)]
        s, t = z[0], z[1]
        F = Fraction(s + t, 2) * Q(*z[2:6]) + Fraction(s - t, 2) * Q(*z[6:10])
        vals.append(F)
    c = list(vals)                                   # Walsh-Hadamard transform, exact
    h = 1
    while h < len(c):
        for i in range(0, len(c), 2 * h):
            for j in range(i, i + h):
                a, b = c[j], c[j + h]
                c[j], c[j + h] = a + b, a - b
        h *= 2
    return vals, c


def witness_poly(z):
    """The cubic expansion read off the definition: (s+t)(x1x3 + x2x3 + x1x4 - x2x4)/4 + (s-t)(y1y3 + ...)/4."""
    s, t, x1, x2, x3, x4, y1, y2, y3, y4 = z
    return Fraction((s + t) * (x1 * x3 + x2 * x3 + x1 * x4 - x2 * x4) +
                    (s - t) * (y1 * y3 + y2 * y3 + y1 * y4 - y2 * y4), 4)


def lemma_checks(rep):
    forms = affine_forms(3)
    worst = Fraction(0)
    rows = []
    for s in range(4):
        ps = sorted(set(p for _, p in forms[s]))
        for p in ps:
            b = Fraction(s, 4) + p / 2 + (1 - p)
            worst = max(worst, b)
            rows.append("s=%d p=%s -> %s" % (s, p, b))
        rep.info("affine {-1,0,1}-valued forms depending on exactly %d sign variables" % s,
                 "%d (of %d truth tables)%s" % (len(forms[s]), 3 ** (2 ** s),
                                                ", p in {%s}" % ", ".join(str(p) for p in ps) if ps else ""))
    rep.check("lemma, finite part: no {-1,0,1}-valued affine form depends on 3 variables (so none on >= 3: fixing "
              "variables of such a form leaves one on 3)", len(forms[3]) == 0)
    rep.check("lemma, finite part: the forms are, up to signs, 1, x_j, (1 + x_j)/2, (x_j + x_k)/2 (2, 2 + 4 and 4 "
              "truth tables for s = 0, 1, 2)",
              [len(forms[s]) for s in range(3)] == [2, 6, 4] and
              set(p for _, p in forms[1]) == {Fraction(1), Fraction(1, 2)} and
              set(p for _, p in forms[2]) == {Fraction(1, 2)})
    rep.check("lemma, finite part: max over the cases of s/4 + p/2 + (1 - p) is 5/4 (cases %s)" % "; ".join(rows),
              worst == Fraction(5, 4))
    vals, c = witness_F()
    rep.check("witness F is Boolean (every one of the 1024 values is +1 or -1)", all(v in (1, -1) for v in vals))
    coef = {S: Fraction(c[S], 1024) for S in range(1024) if c[S]}
    degs = sorted(set(bin(S).count("1") for S in coef))
    used = 0
    for S in coef:
        used |= S
    rep.check("witness F: 16 nonzero monomials, all cubic, every one of the 10 variables in one (Walsh transform)",
              len(coef) == 16 and degs == [3] and used == (1 << 10) - 1,
              "%d monomials, degrees %s, coefficients %s" % (len(coef), degs, sorted(set(coef.values()))))
    ok_poly = True
    for m in range(1 << 10):
        z = [(-1 if (m >> i) & 1 else 1) for i in range(10)]
        if witness_poly(z) != vals[m]:
            ok_poly = False
            break
    rep.check("witness F equals its cubic expansion at all 1024 points (a second evaluation)", ok_poly)
    DF = sum(Fraction(1, 2 ** max(bin(S).count("1") for S in coef if (S >> i) & 1)) for i in range(10))
    rep.check("D(F) = 10 * 2^-3 = 5/4: the cap 5/4 is attained, so no smaller cap can hold", DF == Fraction(5, 4))
    return DF


# ----------------------------------------------------------------------------------------------------------------
# The row
# ----------------------------------------------------------------------------------------------------------------

def run(root, do_controls):
    rep = Report("ROW [P2026] (R. Pinho, PR #258): C44 <= 4.3576549226*, rpinho/c44-certificate @f3fb25d")
    # ---- files and hashes
    paths = {k: os.path.join(root, k) for k in (F_CAND, F_SRC, F_RES, F_EPI_CERT, F_EPI_RES)}
    hs = {k: sha256(p) for k, p in paths.items()}
    for k in (F_CAND, F_SRC, F_RES, F_EPI_CERT, F_EPI_RES):
        rep.info("sha256 " + k, "%s (%d bytes)" % (hs[k], os.path.getsize(paths[k])))
    rep.check("sha256 of c44_candidate.json matches SHA256SUMS @f3fb25d", hs[F_CAND] == H_CAND)
    rep.check("sha256 of wellens-certificate.json matches SHA256SUMS, the paper (Section 4), THIRD_PARTY_NOTICES "
              "and BakeLab's checksums.sha256 (ce66f109...cf49)", hs[F_SRC] == H_SRC)
    rep.check("sha256 of wellens-result.json matches SHA256SUMS and BakeLab's checksums.sha256", hs[F_RES] == H_RES)
    rep.check("the certificate and result are byte-identical to BakeLab/Epi-results @97fa6a4 "
              "mathematics/wellens/", open(paths[F_SRC], "rb").read() == open(paths[F_EPI_CERT], "rb").read() and
              open(paths[F_RES], "rb").read() == open(paths[F_EPI_RES], "rb").read())
    cand = json.load(open(paths[F_CAND]))
    src = json.load(open(paths[F_SRC]))
    res = json.load(open(paths[F_RES]))
    rep.info("c44_candidate.json", json.dumps(cand, sort_keys=True))
    rep.check("candidate's source_sha256 is the source file's hash", cand.get("source_sha256") == hs[F_SRC])
    rep.check("BakeLab's result.json lists the certificate with this hash and 108347 bytes",
              any(a.get("sha256") == hs[F_SRC] and a.get("bytes") == os.path.getsize(paths[F_SRC])
                  for a in res.get("artifacts", [])))
    H = int(cand["horizon"])
    cutoff = int(cand["degree_cutoff"])
    cap3 = rat(cand["degree_three_cap"])
    claim = Fraction(cand["claim"])
    rep.check("candidate fields: horizon 96, degree_cutoff 48, degree_three_cap 5/4, R_3 = 10, claim %s (the 44a "
              "page's bound)" % PAGE_CLAIM,
              H == 96 and cutoff == 48 and cap3 == Fraction(5, 4) and cand.get("external_theorem_R3") == 10 and
              cand["claim"] == PAGE_CLAIM and src.get("degree") == cutoff and src.get("horizon") == H)
    # ---- Epi's exact certificates
    certs = parse_epi(src)
    rep.check("Epi's file has degree-d certificates for exactly d = 1..48", sorted(certs) == list(range(1, 49)))
    epi_ok = {}
    bad = []
    lag_bad = []
    for d in sorted(certs):
        nodes, weights = certs[d]
        ok, why = epi_condition(d, nodes, weights)
        if ok:
            epi_ok[d] = nodes[-1] - 1
        else:
            bad.append("%d: %s" % (d, why))
        if not lagrange_condition(d, nodes, weights):
            lag_bad.append(d)
    rep.check("all 48 certificates: 1 = x_0 < ... < x_d integers, sum_i w_i x_i^j = 0 exactly for j = 1..d, w_0 > 0, "
              "w_0 + sum_{w_i<0} w_i > 0", not bad, "; ".join(bad))
    rep.check("all 48 certificates: w_i x_i prod_{j!=i}(x_i - x_j) is one nonzero constant (the Lagrange form)",
              not lag_bad, "fails at d = %s" % lag_bad if lag_bad else "")
    signs = all(all((w > 0) == (i % 2 == 0) for i, w in enumerate(certs[d][1])) for d in certs)
    rep.info("weights alternate in sign, w_0 > 0", signs)
    B, srcB = bs_bounds(H, epi_ok, cutoff)
    rep.info("B_d = x_d - 1 (Epi), d = 1..48", " ".join(str(epi_ok.get(d, "-")) for d in range(1, 49)))
    rep.check("Epi's B_d equal Wellens's LP Table 1 for d = 1..14 (1 3 6 10 15 21 29 38 47 58 71 84 99 114)",
              [epi_ok.get(d) for d in range(1, 15)] == WELLENS_TABLE1)
    tighter = [d for d in range(1, cutoff + 1) if srcB[d] == "epi"]
    equal = [d for d in range(1, cutoff + 1) if srcB[d] == "both"]
    looser = [d for d in range(1, cutoff + 1) if srcB[d] == "wellens"]
    rep.info("where Epi's bound is below Wellens's 3b(b-1) <= 2(d^4-d^2)", "%d degrees; equal at d = %s; Wellens "
             "smaller at d = %s" % (len(tighter), equal, looser if looser else "none"))
    rep.info("B_d used, d = 49..96 (Wellens)", " ".join(str(B[d]) for d in range(49, H + 1)))
    rep.check("B_d is strictly increasing in d = 1..96", all(B[d] < B[d + 1] for d in range(1, H)))
    rep.info("B_96, sum of B_d", "%d, %d" % (B[H], sum(B)))
    # ---- the recurrence
    caps = caps_list(H, cap3)
    t0 = time.time()
    E, UA = rec_rows(B, H, caps)
    tA = time.time() - t0
    t0 = time.time()
    S, UB = rec_cols(B, H, caps)
    tB = time.time() - t0
    UfA = [Fraction(v, 1 << E) for v in UA]
    UfB = [Fraction(v, 1 << S) for v in UB]
    rep.check("recurrence A (rows, scale 2^97) and B (columns, explicit max, scale 2^195) agree on every U(B_d, d), "
              "d = 1..96", UfA == UfB, "%.1f s and %.1f s" % (tA, tB))
    t0 = time.time()
    UC = rec_literal(B, 30, caps)
    rep.check("recurrence C (literal memoised definition, Fractions) agrees with A for d = 1..30",
              UC[1:] == UfA[1:31], "%.1f s" % (time.time() - t0))
    rep.info("U(B_d, d) for d = 1..12", ", ".join(str(UfA[d]) for d in range(1, 13)))
    rep.info("U(B_d, d) for d = 24, 48, 72, 96", ", ".join(dec(UfA[d], 15) for d in (24, 48, 72, 96)))
    mx = max(range(1, H + 1), key=lambda d: UfA[d])
    rep.check("every degree d <= 96 is covered: max_d U(B_d, d) is attained at d = 96 (so D(f) <= U(B_96, 96) for "
              "deg f <= 96 directly)", all(UfA[d] <= UfA[H] for d in range(1, H + 1)), "argmax d = %d" % mx)
    # ---- the tail
    T = tail_closed(H)
    Ts, part = tail_series(H)
    rep.check("tail identity S(m) - S(m+1) = (2m^2 - m)/2^m, coefficient by coefficient", tail_identity())
    rep.check("T_96 summed term by term to k = 4000, plus the exact remainder S(4001), equals the paper's closed form "
              "((H+1)^3 + 2H^2 + 11H + 19)/2^(H+1) = %s" % pow2_str(T), Ts == T)
    rep.check("the term-by-term sum alone (no closed form) is below T_96 by less than 2^-3900",
              Fraction(97 ** 3, 2 ** 97) + part < T <= Fraction(97 ** 3, 2 ** 97) + part + Fraction(1, 2 ** 3900))
    rep.info("T_96", "%s = %.6e" % (pow2_str(T), float(T)))
    # ---- the bound
    bound = UfA[H] + T
    rep.info("U(B_96, 96)", "%s = %s" % (pow2_str(UfA[H]), dec(UfA[H], 30)))
    rep.info("bound U(B_96, 96) + T_96", "%s = %s..." % (bound, dec(bound, 30)))
    rep.info("as a dyadic", pow2_str(bound))
    rep.check("bound equals the exact value printed on the 44a page and in the paper", bound == PAGE_EXACT,
              "172624496193406565421448533795/39614081257132168796771975168")
    rep.check("claim 4.3576549226 = 21788274613/5000000000 exactly", claim == Fraction(21788274613, 5000000000))
    rep.check("bound < 4.3576549226 (the printed bound holds)", bound < claim,
              "margin %.6e" % float(claim - bound))
    rep.check("claim - bound and T_96 equal the margin and tail in the claimant's c44_exact.log", claim - bound ==
              LOG_MARGIN and T == LOG_TAIL)
    rep.check("the printed bound is the exact value rounded up at 10 decimals", claim == Fraction(dec(bound, 10, up=True))
              and dec(bound, 10) != dec(bound, 10, up=True))
    # ---- the preceding (Epi) bound and the README's horizon-72 value
    b_epi, U_epi, _ = bound_at(H, epi_ok, cutoff, Fraction(3, 2))
    rep.info("cap 3/2 at degree 3 (Wellens's influence cap, Epi's method)", "%s = %s..." % (b_epi, dec(b_epi, 30)))
    rep.check("with cap 3/2 the recurrence gives the paper's reproduced preceding value", b_epi == PAPER_EPI_EXACT,
              "172933981203227910490173314851/39614081257132168796771975168")
    rep.check("Epi's printed 4.3654674227 holds for that value (Epi's bound reproduced)", b_epi < Fraction(EPI_PRINTED),
              "the 10-decimal round-up would be %s; Epi prints %s, %.2e above the exact value" % (
                  dec(b_epi, 10, up=True), EPI_PRINTED, float(Fraction(EPI_PRINTED) - b_epi)))
    rep.check("the gain is exactly 1/128, as the page says", b_epi - bound == Fraction(1, 128))
    diffs = [d for d in range(1, H + 1) if U_epi[d] != UfA[d]]
    rep.info("degrees whose U(B_d, d) the cap changes", "%d of 96 (d = %s ... %s); U(B_3, 3): %s -> %s" % (
        len(diffs), diffs[:3], diffs[-2:], U_epi[3], UfA[3]))
    rep.check("with cap 3/2, U(B_8, 8) = 63/16 = 3.9375: the 44a page's [Wel2019] Table 2 value (R_8 <= 256 * 3.9375 "
              "= 1008)", U_epi[8] == Fraction(63, 16))
    rep.info("by-product with cap 5/4 (not on the page)", "U(B_8, 8) = %s, so R_8 <= %s under the same premises" % (
        UfA[8], 256 * UfA[8]))
    b72, _, _ = bound_at(72, epi_ok, cutoff, cap3)
    rep.info("horizon 72", "%s = %s..." % (b72, dec(b72, 30)))
    rep.check("horizon 72 gives the README's second exact value 10289221775138769568521/2361183241434822606848",
              b72 == README_H72)
    rep.check("horizon 72 bound < 4.3576549226 too", b72 < claim)
    # ---- Wellens's own published values from his Table 1 (a third-party test of the recurrence code)
    T1 = WELLENS_TABLE1
    Bw = [0] + [T1[d - 1] if d <= 14 else d * d for d in range(1, 31)]
    E30, Uw = rec_rows(Bw, 30, caps_list(30, Fraction(3, 2)))
    u30 = Fraction(Uw[30], 1 << E30)
    rep.check("Wellens (3.18) reproduced: Table 1 for d <= 14, B_d = d^2 above, caps d/2: D(30^2, 30) = %s... "
              "(Wellens: 4.4157...), plus T_30 <= 4.4158" % dec(u30, 8),
              dec(u30, 4) == "4.4157" and u30 + tail_closed(30) <= Fraction("4.4158"))
    H2 = 60
    Bw2 = [0] + [min(T1[d - 1], wellens_b(d)) if d <= 14 else wellens_b(d) for d in range(1, H2 + 1)]
    E2, Uw2 = rec_rows(Bw2, H2, caps_list(H2, Fraction(3, 2)))
    b2 = Fraction(Uw2[H2], 1 << E2) + tail_closed(H2)
    rep.check("Wellens's Theorem 3.13 reproduced: Table 1 with his Section 4 inequality gives %s... <= 4.3935 "
              "(his D_inf <= 4.3935, rounded up at 4 decimals)" % dec(b2, 8),
              b2 <= Fraction("4.3935") and dec(b2, 4, up=True) == "4.3935")
    # ---- the lemma's finite parts
    DF = lemma_checks(rep)
    # ---- red controls
    if do_controls:
        controls(rep, epi_ok, certs, cutoff, H, cap3, claim, bound, b_epi, DF)
    rep.info("what rests on mathematics", "r(f) <= D(f) 2^deg; Wellens (3.2), Lemma 3.8, Lemma 3.9, Corollary 3.10 "
             "(3.14)-(3.15) with bs <= deg^2 above 96; D <= I/2 <= deg/2; Wellens's B^2 - B <= (2/3)(d^4 - d^2); the "
             "symmetrization giving bs <= x_d - 1; the paper's Lemma 2 (cap 5/4), whose prose proof rests on the "
             "published R_3 <= 10 (finite cases computed above)")
    return rep.done()


def controls(rep, epi_ok, certs, cutoff, H, cap3, claim, bound, b_epi, DF):
    print("  -- red controls (each must fire)", flush=True)
    low = claim - Fraction(1, 10 ** 10)
    rep.check("control: the printed bound lowered by 1e-10 (4.3576549225) is refuted (exact bound above it)",
              bound > low, "bound - 4.3576549225 = %.4e" % float(bound - low))
    rep.check("control: with the old cap 3/2 the claim 4.3576549226 is not reached (the new cap is load-bearing)",
              b_epi > claim, "%s..." % dec(b_epi, 12))
    # the bound as a function of the degree-3 cap c (exact, at dyadic c)
    thr = Fraction(173, 128)
    sweep = {}
    for cnum in [64 + j for j in range(0, 33)] + [173, 347]:
        c = Fraction(cnum, 64) if cnum <= 96 else Fraction(cnum, 128 if cnum == 173 else 256)
        sweep[c] = bound_at(H, epi_ok, cutoff, c)[0]
    for e in (20, 40):
        for c in (thr - Fraction(1, 2 ** e), thr + Fraction(1, 2 ** e)):
            sweep[c] = bound_at(H, epi_ok, cutoff, c)[0]

    def model(c):                                    # new + clamp(c - 173/128, 0, 1/128)
        return bound + min(max(c - thr, Fraction(0)), Fraction(1, 128))
    rep.check("control: over %d dyadic caps c in [1, 3/2] the exact bound is %s + clamp(c - 173/128, 0, 1/128): "
              "constant for c <= 173/128, slope 1 up to 87/64, Epi's value from 87/64 on" % (len(sweep), "new"),
              all(sweep[c] == model(c) for c in sweep), "checked at c = 173/128 +- 2^-20, +- 2^-40 too")
    rep.check("control: a cap of 11/8 (what R_3 <= 11 would give in the all-cubic case) gives exactly Epi's bound, "
              "above the claim: R_3 <= 10 is load-bearing", sweep[Fraction(11, 8)] == b_epi and b_epi > claim)
    rep.check("control: a cap of 347/256 gives an intermediate exact bound (%s), above the claim" %
              pow2_str(sweep[Fraction(347, 256)]), bound < sweep[Fraction(347, 256)] < b_epi and
              sweep[Fraction(347, 256)] > claim, "%s..." % dec(sweep[Fraction(347, 256)], 12))
    rep.check("control: caps 1 and 9/8 (false: the witness F has D(F) = %s) give the same bound as 5/4, so the "
              "computation cannot see a false cap, and a smaller cap gains nothing" % DF,
              sweep[Fraction(1)] == bound and sweep[Fraction(9, 8)] == bound and DF > Fraction(9, 8))
    # corrupt one weight of Epi's degree-17 certificate
    nodes, weights = certs[17]
    w2 = list(weights)
    w2[5] = w2[5] * (1 + Fraction(1, 10 ** 30))
    ok, why = epi_condition(17, nodes, w2)
    rep.check("control: weight w_5 of the degree-17 certificate times (1 + 1e-30) is rejected (%s)" % why, not ok)
    eo = dict(epi_ok)
    del eo[17]
    bw, _, Bw = bound_at(H, eo, cutoff, cap3)
    rep.check("control: without the degree-17 certificate, B_17 falls back from %d to Wellens's %d and the exact "
              "bound rises above the claim" % (epi_ok[17], Bw[17]), bw > claim, "%s..." % dec(bw, 15))
    nodes, weights = certs[48]
    n2 = list(nodes)
    n2[-1] += 1
    ok, why = epi_condition(48, n2, weights)
    rep.check("control: the last node of the degree-48 certificate raised by one (claiming B_48 one higher) is "
              "rejected (%s)" % why, not ok)
    up, down = [], []
    worst = (Fraction(0), 0)
    for d in range(1, H + 1):
        B0, _ = bs_bounds(H, epi_ok, cutoff)
        bu = bound_at(H, epi_ok, cutoff, cap3, B_override={d: B0[d] + 1})[0]
        if bu != bound:
            up.append(d)
            worst = max(worst, (bu - bound, d))
        if B0[d] > 1:
            bd = bound_at(H, epi_ok, cutoff, cap3, B_override={d: B0[d] - 1})[0]
            if bd != bound:
                down.append(d)
    rep.check("control: raising a single B_d by one changes the exact bound for %d of the 96 degrees" % len(up),
              len(up) > 0, "unchanged for d = %s; largest rise %.3e at d = %d" % (
                  [d for d in range(1, H + 1) if d not in up], float(worst[0]), worst[1]))
    rep.info("lowering a single B_d by one changes the bound for", "%d degrees (unchanged for d = %s)" % (
        len(down), [d for d in range(1, H + 1) if d not in down]))
    bsq, _, _ = bound_at(H, epi_ok, cutoff, cap3, above="square")
    rep.check("control: with bs <= d^2 instead of Wellens's Section 4 inequality for d = 49..96 the bound rises above "
              "the claim (that inequality is load-bearing)", bsq > claim, "%s..." % dec(bsq, 15))
    vals = []
    allok = True
    for HH in (60, 84, 95, 97, 110):
        b_, _, _ = bound_at(HH, epi_ok, cutoff, cap3)
        vals.append("H=%d: %s" % (HH, dec(b_, 15)))
        allok = allok and b_ < claim
    rep.check("other horizons, each a valid bound by the same argument, are all below the claim too", allok,
              "; ".join(vals))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        print(__doc__)
        sys.exit(2)
    ok = run(args[0], "--controls" in sys.argv)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
