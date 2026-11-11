#!/usr/bin/env python3
"""Exact replay of the 58a lower-bound row [Kan2021] (teorth/optimizationproblems, main 87b4ebb): "The A = 4 version
fails (counterexamples d = 54, 150), so any universal bound must satisfy A_Zar >= 5".

For d = 54 and d = 150, every b with 1 <= b < d and gcd(b, d) = 1: the continued fraction b/d = [0; a_1, ..., a_k]
by the Euclidean algorithm (integers).  A finite continued fraction has two expansions, [.., a_k] and
[.., a_k - 1, 1]; we accept b when either has every partial quotient <= A, i.e. a_1..a_{k-1} <= A and a_k <= A + 1
(the weaker test, so a failure is a failure under either convention), and also report the canonical test (a_k <= A).
A counterexample to the A = 4 version is a d with no accepted b.  Anchors: both d pass with A = 5; and the list of
d <= DSCAN failing A = 4 under each convention.  Written from the page's row and definition; no download.
usage: python3 -I zaremba_check.py [--controls]
"""
import sys
from math import gcd

CONTROLS = "--controls" in sys.argv
DSCAN = 3000


def cf(b, d):
    out = []
    while d:
        q, r = divmod(b, d)
        out.append(q)
        b, d = d, r
    return out[1:]              # b/d < 1, so the first quotient is 0


def ok_either(qs, A):
    return all(q <= A for q in qs[:-1]) and qs[-1] <= A + 1


def ok_canonical(qs, A):
    return all(q <= A for q in qs)


def witnesses(d, A, test):
    return [b for b in range(1, d) if gcd(b, d) == 1 and test(cf(b, d), A)]


def main():
    good = True
    for d in (54, 150):
        units = [b for b in range(1, d) if gcd(b, d) == 1]
        w4 = witnesses(d, 4, ok_either)
        w4c = witnesses(d, 4, ok_canonical)
        w5 = witnesses(d, 5, ok_either)
        best = min(max(cf(b, d)[:-1] + [cf(b, d)[-1] - 1]) for b in units)
        print(f"d = {d}: {len(units)} reduced fractions b/d; with every partial quotient <= 4: {len(w4)} (either "
              f"expansion), {len(w4c)} (canonical); with <= 5: {len(w5)} (e.g. b = {w5[0]}: {cf(w5[0], d)}); the "
              f"least possible largest partial quotient is {best}")
        good &= not w4 and not w4c and bool(w5)
    fail_e = [d for d in range(2, DSCAN + 1) if not witnesses(d, 4, ok_either)]
    fail_c = [d for d in range(2, DSCAN + 1) if not witnesses(d, 4, ok_canonical)]
    print(f"anchor: 2 <= d <= {DSCAN} with no reduced b/d, all partial quotients <= 4: {fail_e} (either expansion), "
          f"{fail_c} (canonical)")
    print(f"RESULT [Kan2021]: d = 54 and d = 150 have no reduced fraction with all partial quotients <= 4, so "
          f"A_Zar >= 5: {'HOLDS' if good else 'FAILS'}")
    if CONTROLS:
        res = [("d = 55 (has a witness with A = 4)", bool(witnesses(55, 4, ok_either))),
               ("A = 5 at d = 54 (has witnesses)", bool(witnesses(54, 5, ok_either)))]
        print("red controls (a changed value must not fail):")
        for label, r in res:
            print(f"  {label}: {'as expected' if r else 'NOT AS EXPECTED'}")
        good &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
