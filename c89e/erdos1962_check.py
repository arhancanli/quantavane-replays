#!/usr/bin/env python3
"""Exact replay of the 89a upper-bound row [Erd1962] (teorth/optimizationproblems, main 87b4ebb): C_89 <= 1/32, the
uniformly random colouring.

C_89 is the limit, as n grows, of the least fraction of monochromatic K_4 over 2-colourings of the edges of K_n. A K_4 is
monochromatic in exactly 2 of the 2^6 colourings of its six edges, so by linearity the average over all 2^C(n,2)
colourings of K_n of the number of monochromatic K_4 is C(n,4)/32, and some colouring has at most that many: the least
fraction is <= 1/32 for every n. Checked exactly:
  1. the 64 colourings of K_4: exactly 2 monochromatic;
  2. for n = 4, 5, 6, every one of the 2^C(n,2) colourings (64, 1024, 32768) enumerated: the total number of
     monochromatic K_4 equals 2^C(n,2) C(n,4)/32 exactly, and the least count is listed (it is at most the average);
  3. the step-graphon form: one class with red density 1/2 gives (1/2)^6 + (1/2)^6 = 1/32.
Red controls (--controls): counting a K_4 as monochromatic when five of its six edges agree gives a different average
(refused); a claimed 1/33 from the random colouring is refused (the average is exactly 1/32).
usage: python3 -I erdos1962_check.py [--controls]      (Python 3.9+, standard library only)
"""
import itertools
import sys
from fractions import Fraction as F
from math import comb

CONTROLS = "--controls" in sys.argv


def run(n, mono=lambda cols: len(set(cols)) == 1):
    E = list(itertools.combinations(range(n), 2))
    idx = {e: i for i, e in enumerate(E)}
    quads = [[idx[p] for p in itertools.combinations(q, 2)] for q in itertools.combinations(range(n), 4)]
    total, least = 0, None
    for c in range(1 << len(E)):
        k = sum(1 for qd in quads if mono([(c >> i) & 1 for i in qd]))
        total += k
        least = k if least is None or k < least else least
    return len(E), total, least


def main():
    ok = True
    print("89a row [Erd1962]: C_89 <= 1/32 from the uniformly random 2-colouring")
    m = sum(1 for c in range(64) if len({(c >> i) & 1 for i in range(6)}) == 1)
    print("  1. K_4: %d of the 64 colourings are monochromatic (2/64 = %s)" % (m, F(m, 64)))
    ok &= m == 2
    for n in (4, 5, 6):
        e, total, least = run(n)
        avg = F(total, 2 ** e)
        good = avg == F(comb(n, 4), 32) and least <= avg
        print("  2. n = %d: %d colourings, average monochromatic K_4 = %s = C(%d,4)/32: %s; least %d"
              % (n, 2 ** e, avg, n, avg == F(comb(n, 4), 32), least))
        ok &= good
    g = F(1, 2) ** 6 + F(1, 2) ** 6
    print("  3. one class, red density 1/2: %s" % g)
    ok &= g == F(1, 32)
    print("RESULT [Erd1962]: some colouring of K_n has at most C(n,4)/32 monochromatic K_4 for every n, so "
          "C_89 <= 1/32 = 0.03125: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        e, total, least = run(5, mono=lambda cols: max(cols.count(0), cols.count(1)) >= 5)
        avg = F(total, 2 ** e)
        print("  five of six edges agreeing counted as monochromatic (n = 5): average %s, not C(5,4)/32 = %s: %s"
              % (avg, F(5, 32), "refused, as expected" if avg != F(5, 32) else "UNEXPECTED"))
        c_ok &= avg != F(5, 32)
        claim = F(1, 32) <= F(1, 33)
        print("  a claimed 1/33 from the random colouring: %s" % ("holds (UNEXPECTED)" if claim else "refused, as expected"))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
