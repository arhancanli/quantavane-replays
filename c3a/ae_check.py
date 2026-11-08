#!/usr/bin/env python3
"""Independent replay of the two [GGSWT2025] (AlphaEvolve) rows of 3a, C_3a >= 1.1479 (a set U_1 of 2003 integers) and
C_3a >= 1.1584 (a set U_2 of 54265 integers), from the sets alone.

The sets are the list literals `solution_1` and `solution_2` of the notebook
experiments/sums_differences_problems/sums_differences_problems.ipynb of google-deepmind/
alphaevolve_repository_of_problems at commit 8f44745 (fetched by fetch.sh). They are extracted by parsing that cell
with Python's ast module and evaluating only the two list literals; no code of the notebook is run.

The lemma on the 3a page (from [GHR2007]): a finite set U of non-negative integers containing 0 with
|U - U| < 2 max(U) + 1 gives C_3a >= 1 + log(|U - U|/|U + U|)/log(2 max(U) + 1). We check the hypotheses and count
|U + U| and |U - U| exactly: for U_1 directly with Python sets, and in C with bitsets (bin/ae_count, from
ae_count.c) for both sets; the two routes must agree on U_1 and on random subsets of U_2. The logarithms are
correctly rounded in 50-digit decimal.

  python3 -I ae_check.py NOTEBOOK BINDIR [--controls]
Python 3.9+ and a C compiler."""
import ast
import hashlib
import json
import random
import subprocess
import sys
import time
from decimal import Decimal as D, getcontext

getcontext().prec = 50
NB_SHA = None   # recorded in SHA256SUMS; checked by fetch.sh


def extract(nb_path):
    raw = open(nb_path, "rb").read()
    nb = json.loads(raw.decode("utf-8"))
    found = {}
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        try:
            tree = ast.parse("".join(cell["source"]))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.List):
                for t in node.targets:
                    if isinstance(t, ast.Name) and t.id in ("solution_1", "solution_2"):
                        found[t.id] = ast.literal_eval(node.value)
    return found["solution_1"], found["solution_2"], hashlib.sha256(raw).hexdigest()


def c_counts(U, bindir):
    inp = f"{len(U)}\n" + "\n".join(map(str, U)) + "\n"
    out = subprocess.run([f"{bindir}/ae_count"], input=inp, capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1]), int(out[2])


def py_counts(U):
    S = {a + b for a in U for b in U}
    Dd = {a - b for a in U for b in U}
    return len(S), len(Dd), max(U)


def value(ns, nd, mx):
    return 1 + (D(nd) / D(ns)).ln() / D(2 * mx + 1).ln()


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) < 2:
        print(__doc__)
        return 2
    t0 = time.time()
    U1, U2, sha = extract(args[0])
    bindir = args[1]
    ok = True
    print(f"3a rows [GGSWT2025]: the notebook (sha256 {sha[:16]}...): solution_1 has {len(U1)} integers, solution_2 "
          f"has {len(U2)}")
    for name, U, claim in (("U_1", U1, "1.1479"), ("U_2", U2, "1.1584")):
        hyp = len(set(U)) == len(U) and min(U) == 0 and all(isinstance(x, int) and x >= 0 for x in U)
        ns, nd, mx = c_counts(U, bindir)
        if name == "U_1":
            p = py_counts(U)
            same = p == (ns, nd, mx)
            print(f"  {name}: Python sets and the C bitsets agree: {same}")
            ok &= same
        else:
            rnd = random.Random(3)
            same = True
            for k in (400, 800, 1200):
                sub = sorted(rnd.sample(U, k))
                if 0 not in sub:
                    sub = [0] + sub[1:]
                same &= py_counts(sub) == c_counts(sub, bindir)
            print(f"  {name}: on three random subsets (400, 800, 1200 elements) Python and C agree: {same}")
            ok &= same
        lemma = nd < 2 * mx + 1
        v = value(ns, nd, mx)
        good = hyp and lemma and v >= D(claim) and str(v).startswith(claim)
        print(f"  {name}: distinct non-negative integers containing 0: {hyp}; |U+U| = {ns}, |U-U| = {nd}, max U = {mx}, "
              f"|U-U| < 2 max U + 1 = {2 * mx + 1}: {lemma}")
        print(f"  {name}: 1 + log(|U-U|/|U+U|)/log(2 max U + 1) = {v:.20f}..., which truncates to the row's {claim}: "
              f"{good}")
        ok &= good
    if "--controls" in sys.argv:
        print("Red controls (each must fall below the row's 1.1479):")
        ns, nd, mx = c_counts(U1, bindir)
        V = sorted(U1)[:-1]                                              # the largest element removed
        ns2, nd2, mx2 = c_counts(V, bindir)
        v2 = value(ns2, nd2, mx2)
        good = v2 < D("1.1479")
        print(f"  U_1 without its largest element: |U+U| = {ns2}, |U-U| = {nd2}, value {v2:.10f}: "
              f"{'below, as expected' if good else 'NOT BELOW'}")
        ok &= good
        W = sorted(x * 2 for x in U1)
        ns3, nd3, mx3 = c_counts(W, bindir)
        v3 = value(ns3, nd3, mx3)
        good = (ns3, nd3) == (ns, nd) and v3 < D("1.1479")
        print(f"  U_1 dilated by 2 (the same |U+U| and |U-U|, max U doubled): value {v3:.10f}: "
              f"{'below, as expected' if good else 'NOT BELOW'}")
        ok &= good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
