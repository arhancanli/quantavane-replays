#!/usr/bin/env python3
"""Exact replay of the 5a upper-bound row [HZ26] (teorth/optimizationproblems, main 87b4ebb): C_5a <= 0.9435, with
gamma_0 = sqrt(ab) = 0.94349..., from the eight-pair kernel-boundary system of Hou and Zhao (arXiv 2607.01169v3,
Theorem 2.1 and section 5).  Definitions from the paper's LaTeX (read as text):

  m, L, R >= 1; kernel vectors p^(r) in R^m, symmetric probability vectors; mixing weights lambda_r >= 0 summing to 1;
  boundary vectors w^(r) in R^(Lm), with w_j = 1 for j >= Lm;
  covering inequalities  sum_r lambda_r sum_{i<m} p_i^(r) w_{q+i}^(r) >= 1  for 0 <= q <= Lm;
  a = m sum_r lambda_r sum_i (p_i^(r))^2,   b = 1 + 2((1/m) sum_r lambda_r sum_{j<Lm} (w_j^(r))^2 - L);
  Theorem 2.1: if the covering inequalities hold and b > 0, F(N) <= N^(1/2) + sqrt(ab) N^(1/4) + O(1).

The data: R = 8, m = 32, L = 4 and the eight mixing weights 10^-8 (39490874, ...) are printed in section 5 and read
from the LaTeX; the second kernel vector is the uniform one (section 5: "the uniform kernel vector was added as the
second kernel vector"); the other seven kernels (half vectors over 10^12) and the eight boundary vectors (over 10^9)
are the integer literals of the certificate file, read with Python's ast module (nothing in that file is run); the
repair adds eta = 17/250000000000 (printed in section 5) to every boundary entry (Lemma 4.6).

Checked in exact rationals: weights and kernels nonnegative, normalized, symmetric; all 129 covering inequalities;
the last (q = 128) an equality; the least positive slack and its index against the printed fraction; a and b against
the printed fractions; b > 0; ab < (9435/10000)^2; sqrt(ab) to 30 digits against the printed 0.9434925907135450.
Not checked: Theorem 2.1 itself (mathematics, read).
Red controls (--controls): eta one unit (1/250000000000) smaller leaves some covering inequality violated; the first
weight raised by 10^-8 breaks the normalization; one kernel entry moved off its mirror breaks the symmetry; a claimed
0.9434 is refused (ab > 0.9434^2).
usage: python3 -I hz26_check.py CERTS_ROOT [--controls]
"""
import ast
import os
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv
getcontext().prec = 40


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def read_literals(path, names):
    tree = ast.parse(open(path, encoding="utf-8").read())
    out = {}
    for node in tree.body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in names):
            v = ast.literal_eval(node.value)
            need(isinstance(v, list) and all(isinstance(x, int) for x in v), f"{node.targets[0].id} is not a list of ints")
            out[node.targets[0].id] = v
    need(set(out) == set(names), f"missing literals: {set(names) - set(out)}")
    return out


def read_paper(path):
    t = open(path, encoding="utf-8").read()
    m = re.search(r"R=8,\\qquad m=32,\\qquad L=4,", t)
    need(m is not None, "R, m, L not found")
    lam = re.search(r"\(\\lambda_1,\\ldots,\\lambda_8\)=10\^\{-8\}\(&([0-9,]+),\\\\\s*&([0-9,]+)\)", t)
    need(lam is not None, "the mixing weights were not found")
    lams = [int(x) for x in (lam.group(1) + "," + lam.group(2)).split(",") if x]
    eta = re.search(r"common repair increment is \$\\eta=(\d+)/(\d+)\$", t)
    a = re.search(r"a=\\frac\{(\d+)\}\s*\{(\d+)\}", t)
    b = re.search(r"b=\\frac\{(\d+)\}\s*\{(\d+)\}", t)
    sl = re.search(r"least positive slack occurs at \$q=(\d+)\$ and is\s*\\\[\s*\\frac\{(\d+)\}\s*\{(\d+)\}", t)
    g = re.search(r"\\gamma_0=([0-9.]+)\\ldots\. \\qedhere", t)
    need(all(x is not None for x in (eta, a, b, sl, g)), "a stated value was not found")
    return (lams, Fraction(int(eta.group(1)), int(eta.group(2))), Fraction(int(a.group(1)), int(a.group(2))),
            Fraction(int(b.group(1)), int(b.group(2))), int(sl.group(1)), Fraction(int(sl.group(2)), int(sl.group(3))),
            g.group(1))


def system(lit, lams, eta, m=32, L=4):
    lam = [Fraction(x, 10 ** 8) for x in lams]
    halves = {1: lit["P1_HALF_INTS"]}
    for r in range(3, 9):
        halves[r] = lit[f"P{r}_HALF_INTS"]
    kernels = []
    for r in range(1, 9):
        if r == 2:
            kernels.append([Fraction(1, m)] * m)
        else:
            h = [Fraction(x, 10 ** 12) for x in halves[r]]
            need(len(h) == m // 2, "half kernel of the wrong length")
            kernels.append(h + h[::-1])
    w = [[Fraction(x, 10 ** 9) + eta for x in lit[f"W{r}_INTS"]] for r in range(1, 9)]
    return lam, kernels, w


def evaluate(lam, kernels, w, m=32, L=4):
    n = L * m
    ok_norm = sum(lam) == 1 and all(x >= 0 for x in lam)
    ok_kern = all(sum(p) == 1 and all(x >= 0 for x in p) and all(p[i] == p[m - 1 - i] for i in range(m)) for p in kernels)
    ok_len = all(len(wr) == n for wr in w)
    M = []
    for q in range(n + 1):
        tot = Fraction(0)
        for lr, p, wr in zip(lam, kernels, w):
            tot += lr * sum(p[i] * (wr[q + i] if q + i < n else 1) for i in range(m))
        M.append(tot)
    a = m * sum(lr * sum(x * x for x in p) for lr, p in zip(lam, kernels))
    b = 1 + 2 * (Fraction(1, m) * sum(lr * sum(x * x for x in wr) for lr, wr in zip(lam, w)) - L)
    return ok_norm, ok_kern, ok_len, M, a, b


def main():
    ok = True
    lams, eta, a_st, b_st, q_sl, sl_st, g_st = read_paper(os.path.join(ROOT, "hz26", "main830.tex"))
    names = ["P1_HALF_INTS"] + [f"P{r}_HALF_INTS" for r in range(3, 9)] + [f"W{r}_INTS" for r in range(1, 9)]
    lit = read_literals(os.path.join(ROOT, "hz26", "sidon_certificate_8kernel.py"), names)
    lam, kernels, w = system(lit, lams, eta)
    ok_norm, ok_kern, ok_len, M, a, b = evaluate(lam, kernels, w)
    viol = [q for q, v in enumerate(M) if v < 1]
    slacks = [(v - 1, q) for q, v in enumerate(M) if v > 1]
    least, q_least = min(slacks)
    eq_last = M[-1] == 1
    ab = a * b
    gamma = (Decimal(ab.numerator) / Decimal(ab.denominator)).sqrt()
    good = (ok_norm and ok_kern and ok_len and not viol and eq_last and q_least == q_sl and least == sl_st and a == a_st
            and b == b_st and b > 0 and ab < Fraction(9435, 10000) ** 2 and str(gamma).startswith(g_st))
    print(f"[HZ26] R = 8, m = 32, L = 4; mixing weights 10^-8 {lams} (sum {sum(lam)}); eta = {eta}")
    print(f"  weights nonnegative and summing to 1: {ok_norm}; the 8 kernels (the 2nd uniform) symmetric probability "
          f"vectors: {ok_kern}; boundary vectors of length 128: {ok_len}")
    print(f"  covering inequalities q = 0..128: violated {len(viol)}; q = 128 an equality: {eq_last}; least positive slack "
          f"at q = {q_least}: {least} (printed: q = {q_sl}, {sl_st}): {q_least == q_sl and least == sl_st}")
    print(f"  a = {a} (printed equal: {a == a_st})")
    print(f"  b = {b} (printed equal: {b == b_st}; b > 0)")
    print(f"  ab < (9435/10000)^2: {ab < Fraction(9435, 10000) ** 2}; sqrt(ab) = {str(gamma)[:32]} (printed "
          f"{g_st}...)")
    print(f"RESULT [HZ26]: the kernel-boundary system is feasible with a, b as stated, so C_5a <= sqrt(ab) = "
          f"{str(gamma)[:18]}... < 0.9435: {'HOLDS' if good else 'FAILS'}")
    ok &= good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        lam2, k2, w2 = system(lit, lams, eta - Fraction(1, 250000000000))
        _, _, _, M2, _, _ = evaluate(lam2, k2, w2)
        v2 = [q for q, v in enumerate(M2) if v < 1]
        print(f"  eta - 1/250000000000: violated covering inequalities at q = {v2[:6]}{'...' if len(v2) > 6 else ''}")
        res.append(("eta one unit smaller", len(v2) > 0))
        lams3 = list(lams)
        lams3[0] += 1
        lam3, k3, w3 = system(lit, lams3, eta)
        res.append(("first weight + 10^-8", not evaluate(lam3, k3, w3)[0]))
        k4 = [list(p) for p in kernels]
        k4[0][0] += Fraction(1, 10 ** 12)
        k4[0][1] -= Fraction(1, 10 ** 12)
        res.append(("an asymmetric kernel", not evaluate(lam, k4, w)[1]))
        res.append(("a claimed 0.9434", not (ab < Fraction(9434, 10000) ** 2)))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
