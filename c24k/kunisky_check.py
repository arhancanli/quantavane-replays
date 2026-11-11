#!/usr/bin/env python3
"""Exact replay of the 24a lower-bound row [Kun2023] (teorth/optimizationproblems, main 87b4ebb): C_24 >= 1 + sqrt 2,
from the explicit matrices of Kunisky (arXiv 2111.02974v1, Theorem 1.1 and its second proof): A_0 = [1] and
    A_{k+1} = (1/sqrt 2) [[A_k, I], [A_k, -I]]      (2^k x 2^k matrices with columns of unit Euclidean norm),
for which disc(A_k) = 2^(-k/2) + sum_{a=1}^k 2^(-a/2), increasing to sum_{a>=1} 2^(-a/2) = 1 + sqrt 2.
The page's definition: disc(A) = min over x in {-1,1}^n of ||Ax||_inf; C_24 = sup_n K_n, K_n = sup disc(A) over n x n
matrices whose columns have Euclidean norm at most 1.  So each A_k gives C_24 >= disc(A_k).

Arithmetic is exact in Z[sqrt 2]: M_k = 2^(k/2) A_k has entries p + q sqrt 2 with integers p, q
(M_{k+1} = [[M_k, 2^(k/2) I], [M_k, -2^(k/2) I]]).  For k = 0, ..., 4:
  1. every column of M_k has squared norm exactly 2^k (so A_k has unit columns);
  2. min over all 2^(2^k) sign vectors x of max_i |(M_k x)_i|, exhaustively (65,536 vectors at k = 4); values are
     ranked as floats, which is exact here because two different elements p + q sqrt 2 with |p|, |q| <= 64 differ by
     at least 1/(64 + 64 sqrt 2) > 0.006, and equal absolute values are compared as integer pairs;
  3. the result equals 1 + sum_{j<k} 2^(j/2), i.e. disc(A_k) = 2^(-k/2) + sum_{a=1}^k 2^(-a/2) (Theorem 1.1), and the
     step disc(A_{k+1}) = (1 + disc(A_k))/sqrt 2 of the second proof holds at each k < 4.
Beyond k = 4 the theorem's induction (two lines, read) carries the formula; the values increase (the step is
2^(-k/2)(sqrt 2 - 1) > 0) to the sum 1/(sqrt 2 - 1) = 1 + sqrt 2, the row's bound.  At k = 1 the matrix is the 2 x 2
Hadamard matrix over sqrt 2, the matrix of the 10c row [Band2024] (replayed in ../c10c).
Red controls (--controls): a claimed disc(A_4) one unit of 1/4 larger is refused; A_2 with one entry's sign flipped has
a smaller discrepancy; the recursion with +I in both blocks gives repeated columns and a discrepancy 0 at k = 1.
usage: python3 -I kunisky_check.py CERTS_ROOT [--controls]
"""
import itertools
import math
import os
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv
KMAX = 4
S2 = math.sqrt(2)


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def pow_sqrt2(k):
    """2^(k/2) as (p, q)."""
    return (2 ** (k // 2), 0) if k % 2 == 0 else (0, 2 ** (k // 2))


def build(kmax, minus=True):
    M = [[(1, 0)]]
    out = [M]
    for k in range(kmax):
        n = len(M)
        c = pow_sqrt2(k)
        top = [row + [c if j == i else (0, 0) for j in range(n)] for i, row in enumerate(M)]
        sgn = -1 if minus else 1
        bot = [row + [(sgn * c[0], sgn * c[1]) if j == i else (0, 0) for j in range(n)] for i, row in enumerate(M)]
        M = top + bot
        out.append(M)
    return out


def col_norms(M):
    n = len(M)
    res = set()
    for j in range(n):
        p = sum(M[i][j][0] ** 2 + 2 * M[i][j][1] ** 2 for i in range(n))
        q = sum(2 * M[i][j][0] * M[i][j][1] for i in range(n))
        res.add((p, q))
    return res


def disc(M):
    """min over sign vectors of max_i |(Mx)_i|, as an exact pair (p, q) with p + q sqrt 2 >= 0."""
    n = len(M)
    rows = [[(j, e) for j, e in enumerate(r) if e != (0, 0)] for r in M]
    best, best_pair = None, None
    for x in itertools.product((1, -1), repeat=n):
        worst, worst_pair = -1.0, None
        for r in rows:
            p = sum(x[j] * e[0] for j, e in r)
            q = sum(x[j] * e[1] for j, e in r)
            v = p + q * S2
            if v < 0:
                v, p, q = -v, -p, -q
            if v > worst + 1e-9:
                worst, worst_pair = v, (p, q)
            elif abs(v - worst) <= 1e-9:
                need((p, q) == worst_pair, "two different elements of Z[sqrt 2] within 1e-9")
        if best is None or worst < best - 1e-9:
            best, best_pair = worst, worst_pair
        elif abs(worst - best) <= 1e-9:
            need(worst_pair == best_pair, "two different elements of Z[sqrt 2] within 1e-9")
    return best_pair


def formula_pair(k):
    """1 + sum_{j<k} 2^(j/2) as (p, q)."""
    p, q = 1, 0
    for j in range(k):
        a, b = pow_sqrt2(j)
        p, q = p + a, q + b
    return (p, q)


def check_source(path):
    text = open(path, encoding="utf-8").read()
    need("$K(2^k) \\geq 2^{-k/2} + \\sum_{a = 1}^k 2^{-a/2}$." in text, "Theorem 1.1 was not found")
    need("K \\geq \\sum_{a = 1}^{\\infty} 2^{-a/2} = 1 + \\sqrt{2}" in text, "the limit statement was not found")
    need("\\whbA^{T_{k + 1}} = \\frac{1}{\\sqrt{2}}\\left[\\begin{array}{rr} \\whbA^{T_k} & \\bm I_{2^{k}} \\\\ "
         "\\whbA^{T_k} & -\\bm I_{2^k} \\end{array}\\right]" in text, "the recursion was not found")


def main():
    ok = True
    check_source(os.path.join(ROOT, "kun2023", "unsatisfiable-matrices.tex"))
    print("the paper's Theorem 1.1, its limit statement and the recursion of its second proof are in the source: True")
    Ms = build(KMAX)
    prev = None
    for k, M in enumerate(Ms):
        n = len(M)
        norms = col_norms(M)
        d = disc(M)
        want = formula_pair(k)
        val = (d[0] + d[1] * S2) / 2 ** (k / 2)
        step_ok = prev is None or abs(val - (1 + prev) / S2) < 1e-12
        good = norms == {(2 ** k, 0)} and d == want and step_ok
        ok &= good
        print(f"k = {k}: {n} x {n}; column squared norms of M_k all {sorted(norms)} (= 2^k: {norms == {(2 ** k, 0)}}); "
              f"min_x max_i |(M_k x)_i| = {d[0]} + {d[1]} sqrt 2 over {2 ** n} sign vectors (formula: {want[0]} + "
              f"{want[1]} sqrt 2); disc(A_k) = {val:.12f}; step (1 + disc(A_(k-1)))/sqrt 2: {step_ok}")
        prev = val
    lim = 1 + S2
    print(f"  limit of 2^(-k/2) + sum_(a<=k) 2^(-a/2): 1/(sqrt 2 - 1) = 1 + sqrt 2 = {lim:.12f}")
    print(f"RESULT [Kun2023]: disc(A_k) = 2^(-k/2) + sum_(a=1)^k 2^(-a/2) exactly for k <= {KMAX} (exhaustive), the "
          f"paper's induction beyond, increasing to 1 + sqrt 2, so C_24 >= 1 + sqrt 2: {'HOLDS' if ok else 'FAILS'}")
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        d4 = disc(Ms[4])
        res.append(("a claimed disc(A_4) + 1/4", (d4[0] + 1, d4[1]) != formula_pair(4)))
        M2 = [list(r) for r in Ms[2]]
        M2[0][0] = (-M2[0][0][0], -M2[0][0][1])
        dflip = disc(M2)
        print(f"  A_2 with one sign flipped: min max = {dflip[0]} + {dflip[1]} sqrt 2 (A_2: 2 + 1 sqrt 2)")
        res.append(("A_2 with a flipped sign", dflip != formula_pair(2)))
        Mp = build(1, minus=False)[1]
        dp = disc(Mp)
        print(f"  recursion with +I in both blocks at k = 1: columns {[tuple(r[j] for r in Mp) for j in range(2)]}, "
              f"min max = {dp}")
        res.append(("+I in both blocks", dp == (0, 0)))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
