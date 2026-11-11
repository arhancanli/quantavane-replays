#!/usr/bin/env python3
"""Exact replay of the 10c lower-bound row [Band2024] (teorth/optimizationproblems, main 87b4ebb): C_10c >= sqrt 2 from
"the 2 by 2 sign matrix with rows (1,1) and (1,-1)", read from the page (pinned by ../c10c/fetch.sh).  As in
c10c_check.py, an n x n matrix A with entries in [-1,1] and disc(A) = d gives C_10c >= d/sqrt(n).

Checked: the matrix from the row's text; disc(A) = min over the 4 sign vectors of ||Ax||_inf = 2 exactly (integers);
2/sqrt 2 = sqrt 2, and the printed 1.414214 against sqrt 2 by squaring in integers (it is sqrt 2 rounded up, next to
the exact expression).  Anchor: the page's trivial row 1 (A = [1], disc 1).  The same matrix divided by sqrt 2 is the
k = 1 matrix of Kunisky's 24a construction (../c24k).
Red controls (--controls): the all-ones 2 x 2 matrix has discrepancy 0; a claimed discrepancy 3 is refused; the row's
matrix with a 0 in place of the -1 has discrepancy 1 (so only 1/sqrt 2).
usage: python3 -I band2024_check.py CERTS_ROOT [--controls]
"""
import itertools
import os
import re
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def disc(A):
    n = len(A[0])
    return min(max(abs(sum(a * x for a, x in zip(row, xs))) for row in A) for xs in itertools.product((1, -1), repeat=n))


def main():
    ok = True
    text = open(os.path.join(ROOT, "page", "10c.md"), encoding="utf-8").read()
    row = [l for l in text.splitlines() if "[Band2024]" in l and "2 by 2 sign matrix" in l]
    need(len(row) == 1, "the [Band2024] row was not found")
    m = re.search(r"rows \$\((-?1),\s*(-?1)\)\$ and \$\((-?1),\s*(-?1)\)\$", row[0])
    need(m is not None, "the two rows were not found")
    a, b, c, d = (int(x) for x in m.groups())
    A = [[a, b], [c, d]]
    dec = re.search(r"\\approx\s*([0-9.]+)", row[0]).group(1)
    D = disc(A)
    # printed decimal vs sqrt(2): compare squares in integers; D / sqrt(n) = 2/sqrt(2) = sqrt(2)
    num = int(dec.replace(".", ""))
    scale = 10 ** (len(dec.split(".")[1]))
    above = num * num > 2 * scale * scale
    below_next = (num - 1) ** 2 < 2 * scale * scale
    good = A == [[1, 1], [1, -1]] and D == 2 and D * D == 2 * 2 and above and below_next
    print(f"[Band2024] A = {A} (from the row); disc(A) over the 4 sign vectors = {D}; bound D/sqrt(2) = sqrt 2 "
          f"= 1.41421356237...; printed {dec}: {dec}^2 > 2: {above}, ({num - 1}/{scale})^2 < 2: {below_next} "
          f"(sqrt 2 rounded up, printed with \\approx next to the exact sqrt 2)")
    print(f"RESULT [Band2024]: the 2 x 2 sign matrix has discrepancy 2, so C_10c >= 2/sqrt 2 = sqrt 2: "
          f"{'HOLDS' if good else 'FAILS'}")
    ok &= good
    print(f"anchor (the page's trivial row): A = [1], disc = {disc([[1]])}, bound 1")
    if CONTROLS:
        print("red controls (each must fail):")
        res = [("all-ones 2 x 2", disc([[1, 1], [1, 1]]) != 2),
               ("a claimed discrepancy 3", D != 3),
               ("-1 replaced by 0", disc([[1, 1], [1, 0]]) != 2)]
        print(f"  all-ones: disc {disc([[1, 1], [1, 1]])}; with 0 for -1: disc {disc([[1, 1], [1, 0]])}")
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
