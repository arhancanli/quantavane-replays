#!/usr/bin/env python3
"""Exact replay of three lower-bound rows of the 3d page (teorth/optimizationproblems, main 87b4ebb), where
C(A) = log(|A+A|/|A|) / log(|A-A|/|A|) and C_3d = sup_A C(A) (the page's definition):

  [FrPi73]    the 17-element set printed in the row: |A+A| = 59, |A-A| = 55, C(A) = log(59/17)/log(55/17) = 1.059793...
  [PeWe13]    the family Q_j of Penman and Wells (Integers 13 (2013) #A57): |Q_j| = 5j+17, |Q_j+Q_j| = 32j+63,
              |Q_j-Q_j| = 26j+61 (their Corollary 13, as the row quotes it), so sup_j C(Q_j) >= lim C(Q_j) =
              log(32/5)/log(26/5) = 1.125944426...
  [GGSWT2025] AlphaEvolve's set A = A_1 u A_2, A_1 = {-159, ..., 111} and A_2 of 38 elements: C(A) ~ 1.1219.

Sumsets and difference sets are computed exactly as Python sets of integers (and again by a sorted-merge count for
the two printed sets); logarithms by Python's decimal module (correctly rounded) at 60 digits, and every comparison
with a printed decimal has a margin far above the rounding error.  The AlphaEvolve list is read from the notebook's
first code cell with the ast module as a literal: no notebook code is run.  The family Q_j is transcribed from the
paper (the four pieces of its definition) and checked against all three cardinality formulas for 1 <= j <= JMAX, with
the restricted sumset |Q_j ^+ Q_j| = 32j+56 (j >= 2) and the paper's value fhat(Q_19) = 1.028377107... as anchors.
usage: python3 -I c3d_check.py CERTS_ROOT [--controls]
"""
import ast
import json
import re
import sys
from decimal import Decimal, getcontext

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv
JMAX = 120
getcontext().prec = 60


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def sizes(A):
    A = sorted(set(A))
    S = {a + b for a in A for b in A}
    D = {a - b for a in A for b in A}
    return len(A), len(S), len(D)


def sizes_merge(A):
    """Second route: count distinct sums and differences from sorted lists of all pair values."""
    A = sorted(set(A))
    sums = sorted(a + b for i, a in enumerate(A) for b in A[i:])
    diffs = sorted(a - b for a in A for b in A)
    distinct = lambda xs: sum(1 for i, x in enumerate(xs) if i == 0 or x != xs[i - 1])
    return len(A), distinct(sums), distinct(diffs)


def C(n, s, d):
    return (Decimal(s) / n).ln() / (Decimal(d) / n).ln()


def leading(x, digits):
    """The first `digits` decimals of x > 0, truncated."""
    q = int(x * (10 ** digits))
    return f"{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


def rounded(x, digits):
    q = int((x * (10 ** digits)).to_integral_value())
    return f"{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


def frpi_set(path):
    text = open(path, encoding="utf-8").read()
    m = re.search(r"\$A = \\\\\{([0-9, ]+)\\\\\}\$, which has \$\\lvert A\+A\\rvert = (\d+)\$ and \$\\lvert A-A\\rvert = (\d+)\$", text)
    need(m is not None, "the [FrPi73] set was not found on the page")
    return [int(x) for x in m.group(1).split(",")], int(m.group(2)), int(m.group(3))


def ae_set(path):
    nb = json.load(open(path, encoding="utf-8"))
    src = "".join(nb["cells"][1]["source"])
    src = "\n".join(l for l in src.split("\n") if not l.lstrip().startswith(("%", "!")))
    node = next(n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "best_list")
    vals = ast.literal_eval(node.value)
    need(all(isinstance(v, int) for v in vals), "best_list is not a list of integers")
    return vals


def Q(j):
    """Penman and Wells: Q_j = {0,2,4,12} u {1,5,...,1+4(4j+8)} u {24,40,...,8+16j} u {4+16(j+1), 12+16(j+1),
    14+16(j+1), 16(j+2)}."""
    s = {0, 2, 4, 12}
    s |= set(range(1, 1 + 4 * (4 * j + 8) + 1, 4))
    s |= set(range(24, 8 + 16 * j + 1, 16))
    s |= {4 + 16 * (j + 1), 12 + 16 * (j + 1), 14 + 16 * (j + 1), 16 * (j + 2)}
    return s


def family(qfun, jmax, quiet=False):
    """Check the three formulas (and the restricted-sum anchor) for 1 <= j <= jmax; return the first failure or None."""
    for j in range(1, jmax + 1):
        A = qfun(j)
        n, s, d = sizes(A)
        r = len({a + b for a in A for b in A if a != b})
        if (n, s, d) != (5 * j + 17, 32 * j + 63, 26 * j + 61) or (j >= 2 and r != 32 * j + 56):
            return j, (n, s, d, r)
    return None


def main():
    ok = True
    # [FrPi73]
    A, s_claim, d_claim = frpi_set(f"{ROOT}/page/3d.md")
    n, s, d = sizes(A)
    need(sizes_merge(A) == (n, s, d), "the two routes disagree on the [FrPi73] set")
    c = C(n, s, d)
    good = (n, s, d) == (17, s_claim, d_claim) and leading(c, 6) == "1.059793"
    ok &= good
    print(f"[FrPi73]: A from the page, {n} elements; |A+A| = {s}, |A-A| = {d} (the page: {s_claim}, {d_claim}; two "
          f"routes agree); C(A) = log({s}/{n})/log({d}/{n}) = {leading(c, 20)}...")
    print(f"RESULT [FrPi73]: C_3d >= {leading(c, 12)}...; the printed 1.059793... is its leading digits: "
          f"{'HOLDS' if good else 'FAILS'}")
    # [PeWe13]
    fail = family(Q, JMAX)
    lim = (Decimal(32) / 5).ln() / (Decimal(26) / 5).ln()
    cj = {j: C(5 * j + 17, 32 * j + 63, 26 * j + 61) for j in (1, 10, 100, 10 ** 3, 10 ** 6, 10 ** 12)}
    incr = all(C(5 * j + 17, 32 * j + 63, 26 * j + 61) < C(5 * j + 22, 32 * j + 95, 26 * j + 87) for j in range(1, 2000))
    below = all(v < lim for v in cj.values())
    n10, s10, d10 = sizes(Q(10))
    f10 = Decimal(s10).ln() / Decimal(d10).ln()
    q19 = Q(19)
    r19 = len({a + b for a in q19 for b in q19 if a != b})
    fh19 = Decimal(r19).ln() / Decimal(sizes(q19)[2]).ln()
    good = fail is None and leading(lim, 9) == "1.125944426" and incr and below and leading(fh19, 9) == "1.028377107"
    ok &= good
    print(f"[PeWe13]: Q_j built from the paper's definition for 1 <= j <= {JMAX}: |Q_j| = 5j+17, |Q_j+Q_j| = 32j+63, "
          f"|Q_j-Q_j| = 26j+61 for every j{' (all hold)' if fail is None else f'; FIRST FAILURE at j = {fail[0]}: {fail[1]}'}"
          f"; anchor |Q_j ^+ Q_j| = 32j+56 for j >= 2 and fhat(Q_19) = log {r19}/log {sizes(q19)[2]} = {leading(fh19, 12)}... "
          f"(the paper: 1.028377107...)")
    print("  C(Q_j) from the formulas: " + ", ".join(f"j = {j}: {leading(v, 12)}" for j, v in cj.items())
          + f"; increasing in j for j < 2000: {incr}; each below the limit: {below}")
    print(f"  reader's note: f(Q_10) = log|Q_10+Q_10|/log|Q_10-Q_10| = log {s10}/log {d10} = {leading(f10, 12)}..., "
          f"while the page (quoting the paper's Theorem 20) prints 1.030597781")
    print(f"RESULT [PeWe13]: C_3d >= lim_j C(Q_j) = log(32/5)/log(26/5) = {leading(lim, 15)}...; the printed "
          f"1.125944426 is its leading digits: {'HOLDS' if good else 'FAILS'}")
    # [GGSWT2025]
    A = ae_set(f"{ROOT}/ae/sums_differences_problems.ipynb")
    U = sorted(set(A))
    a1 = set(range(-159, 112))
    a2 = sorted(set(U) - a1)
    n, s, d = sizes(A)
    need(sizes_merge(A) == (n, s, d), "the two routes disagree on the AlphaEvolve set")
    c = C(n, s, d)
    good = len(A) == len(U) and a1 <= set(U) and len(a2) == 38 and rounded(c, 4) == "1.1219"
    ok &= good
    print(f"[GGSWT2025]: the notebook's list has {len(A)} entries, {len(U)} distinct; it contains A_1 = "
          f"{{-159, ..., 111}} ({len(a1)} elements): {a1 <= set(U)}; the other {len(a2)} elements: {a2}")
    print(f"  |A+A| = {s}, |A-A| = {d} (two routes agree); C(A) = log({s}/{n})/log({d}/{n}) = {leading(c, 20)}...")
    print(f"RESULT [GGSWT2025]: C_3d >= {leading(c, 12)}...; the printed ~1.1219 is its value to four decimals: "
          f"{'HOLDS' if good else 'FAILS'}")
    if CONTROLS:
        print("red controls (each must be refused or fail):")
        res = []
        A0, _, _ = frpi_set(f"{ROOT}/page/3d.md")
        B = [x for x in A0 if x != 29] + [30]
        nb, sb, db = sizes(B)
        res.append(("[FrPi73] with 30 for 29", (sb, db) != (59, 55)))
        print(f"  [FrPi73] set with 30 for 29: |A+A| = {sb}, |A-A| = {db}, not (59, 55)")
        res.append(("[FrPi73] printed 1.059794", leading(C(17, 59, 55), 6) != "1.059794"))
        def Qbad(j):
            s = Q(j)
            s.discard(12)
            s.add(10)
            return s
        fb = family(Qbad, 30)
        print(f"  Q_j with 10 for 12: first formula failure at j = {fb[0] if fb else None}")
        res.append(("Q_j with 10 for 12", fb is not None))
        res.append(("[PeWe13] printed 1.125944427", leading(lim, 9) != "1.125944427"))
        Ab = [x for x in A if x != 386]
        nb, sb, db = sizes(Ab)
        cb = C(nb, sb, db)
        print(f"  AlphaEvolve set without 386: C = {leading(cb, 10)}..., not ~1.1219")
        res.append(("AlphaEvolve set without 386", rounded(cb, 4) != "1.1219"))
        for label, r in res:
            print(f"  {label}: {'as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
