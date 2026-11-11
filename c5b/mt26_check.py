#!/usr/bin/env python3
"""Independent replay of the 5b upper-bound row [MT26] (teorth/optimizationproblems, main 87b4ebb): C_5b <= 4/7, from
the 14-point set A_base printed in Ma and Tang, arXiv 2602.23282v2, Lemma 6.1 (section 6.1), read from the LaTeX
source as text.  The paper's code repository was neither read nor run.

Definitions as in c5b_check.py (the 5b page): A is a (4,5)-set if every four of its elements determine at least five
distinct differences; h(A) is the size of a largest Sidon subset of A; f(n) = min h(A) over (4,5)-sets of size n.
The row's bound is c_* = inf_n f(n)/n <= f(14)/14 <= h(A_base)/14 = 8/14 = 4/7 (the identity c_* = inf f(n)/n is
the paper's Theorem 1.4, mathematics, read and not checked).

Computed here, with the two (4,5) routes and the two h routes of c5b_check.py (imported from this folder by path):
  1. A_base: 14 distinct integers, as printed;
  2. the (4,5) property by route (a) (all 1001 four-subsets) and route (b) (equal-difference pairs of pairs);
  3. h(A_base) = 8 by route (a) (branch and bound) and route (b) (every 8- and 9-subset, sum test);
  4. the Sidon subset the proof prints, {0,136,200,243,246,298,323,528}: inside A_base, Sidon, of size 8;
  5. 8/14 = 4/7 exactly, the row's value.
Red controls (--controls): 1056 replaced by 0 + 528 - 136 = 392 makes a parallelogram {0,136,392,528} (refused by both
(4,5) routes); h claimed 7 is refused (route b finds a Sidon 8-subset); h claimed 9 is refused (no Sidon 9-subset);
the printed subset with 528 replaced by 596 is not Sidon; 4/7 claimed from 13 of the 14 points (A_base without 1056)
is refused, since 8/13 > 4/7 and h of that set is 8.

Usage: python3 -I mt26_check.py CERTS_ROOT [--controls]     (CERTS_ROOT holds mt26/757_submit_v2.tex)
"""
import importlib.util
import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("c5b_check", os.path.join(HERE, "c5b_check.py"))
c5b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(c5b)


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def read_paper(path):
    text = open(path, encoding="utf-8").read()
    m = re.search(r"A_\{\\mathrm\{base\}\}\s*=\\\{([0-9,\s]+)\\\}\\subset\\mathbb Z", text)
    need(m is not None, "A_base was not found in the LaTeX source")
    a = [int(x) for x in m.group(1).replace(" ", "").split(",")]
    m2 = re.search(r"outputs the explicit Sidon subset\s*\\\[\s*\\\{([0-9,\s]+)\\\}\.\s*\\\]", text)
    need(m2 is not None, "the printed Sidon subset was not found")
    s = [int(x) for x in m2.group(1).replace(" ", "").split(",")]
    m3 = re.search(r"\\frac\{f\(14\)\}\{14\}\\le \\frac\{8\}\{14\}=\\frac\{4\}\{7\}", text)
    need(m3 is not None, "the final inequality f(14)/14 <= 8/14 = 4/7 was not found")
    return a, s


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    root = sys.argv[1]
    controls = "--controls" in sys.argv[2:]
    a, s = read_paper(os.path.join(root, "mt26", "757_submit_v2.tex"))
    ok = True
    print(f"A_base from the paper ({len(a)} elements): {a}")
    distinct = len(set(a)) == len(a) == 14
    print(f"  14 distinct integers: {distinct}")
    a = sorted(a)
    good, ok45a, ok45b, ha, hb_ok, lines = c5b.check_one(a, 14, stated=8, label="[MT26] A_base,")
    for line in lines:
        print(line)
    sub_ok = set(s) <= set(a) and len(s) == 8 and len(set(s)) == 8 and c5b.sidon_by_sums(sorted(s))
    print(f"  the printed Sidon subset {s}: inside A_base, 8 distinct elements, Sidon (sum test): {sub_ok}")
    ratio = Fraction(ha, 14)
    print(f"  h(A_base)/14 = {ha}/14 = {ratio} (the row: 4/7)")
    ok = distinct and good and sub_ok and ratio == Fraction(4, 7)
    print(f"RESULT [MT26]: A_base is a (4,5)-set (both routes) with h(A_base) = 8 (both routes), so "
          f"C_5b = inf f(n)/n <= 8/14 = 4/7: {'HOLDS' if ok else 'FAILS'}")
    if controls:
        print("red controls (each must be refused or fail):")
        res = []
        bad = sorted(set(a) - {1056} | {392})
        b45a, wa = c5b.is45_route_a(bad)
        b45b, wb = c5b.is45_route_b(bad)
        print(f"  1056 -> 392 (0 + 528 - 136): route a {'holds' if b45a else 'refuses at ' + str(wa)}, "
              f"route b {'holds' if b45b else 'refuses at ' + str(wb)}")
        res.append(("parallelogram {0,136,392,528}", not b45a and not b45b))
        w7, l7 = c5b.h_route_b(a, 7)
        res.append(("h claimed 7", l7 is not None))
        w9, l9 = c5b.h_route_b(a, 9)
        res.append(("h claimed 9", w9 is None))
        print(f"  h claimed 7: a Sidon 8-subset {'exists ' + str(list(l7)) if l7 else 'NOT FOUND'}; "
              f"h claimed 9: a Sidon 9-subset {'FOUND' if w9 else 'does not exist'}")
        s2 = sorted(set(s) - {528} | {596})
        res.append(("printed subset with 596 for 528", not c5b.sidon_by_sums(s2)))
        a13 = sorted(set(a) - {1056})
        h13, _ = c5b.h_route_a(a13)
        print(f"  A_base without 1056: h = {h13}, ratio {Fraction(h13, 13)} > 4/7: {Fraction(h13, 13) > Fraction(4, 7)}")
        res.append(("4/7 from 13 points", Fraction(h13, 13) > Fraction(4, 7)))
        for label, r in res:
            print(f"  {label}: {'refused, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
