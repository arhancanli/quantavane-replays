#!/usr/bin/env python3
"""Exact replay of the 9a lower-bound row [BMRRST1971] (teorth/optimizationproblems, main 87b4ebb):
Theta(C_7) >= 343^(1/5) ~ 3.2141, from an independent set of size 343 = 7^3 in the strong power C_7^(x5).

The paper of Baumert, McEliece, Rodemich, Rumsey, Stanley and Taylor (1971) was not accessible to us and the page prints
no set, so the set is ours, found by exhaustive search over systematic generator matrices: the linear code
S = {m G : m in Z_7^3} with G = [I_3 | P], P = ((1,2),(1,4),(3,1)) (the first P in lexicographic order that works).
Two words of C_7^(x5) are adjacent or equal iff they differ by -1, 0 or +1 (mod 7) in every coordinate. Checked
exactly, by two routes:
  1. every nonzero codeword (342 of them) has a coordinate outside {0, 1, 6}, so no two distinct words of S are
     adjacent; S has 343 distinct words;
  2. all C(343, 2) = 58653 pairs of words compared directly.
Then Theta(C_7) >= alpha(C_7^(x5))^(1/5) >= 343^(1/5) = 7^(3/5) = 3.21410... (the bound Theta(G) >= alpha(G^n)^(1/n)
is Shannon's, read). The printed 3.2141 is the exact value 3.2140958... rounded up, with an approximately-equal sign
next to the exact 343^(1/5).
Red controls (--controls): P with (3,1) replaced by (1,1) has a nonzero codeword in {0, 1, 6}^5 and is refused; the
word (1,1,0,0,0) added to S is adjacent to 0 and refused; a claim of 344 words from this code is refused.
usage: python3 -I bmrrst_check.py [--controls]      (Python 3.9+, standard library only)
"""
import itertools
import sys
from decimal import Decimal, getcontext

CONTROLS = "--controls" in sys.argv
NEAR = {0, 1, 6}


def code(P):
    words = []
    for m in itertools.product(range(7), repeat=3):
        words.append(tuple(list(m) + [sum(m[i] * P[i][j] for i in range(3)) % 7 for j in range(2)]))
    return words


def route1(P):
    bad = [w for w in code(P) if any(w) and all(x in NEAR for x in w)]
    return not bad, bad[:1]


def adjacent_or_equal(x, y):
    return all((a - b) % 7 in NEAR for a, b in zip(x, y))


def route2(words):
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            if adjacent_or_equal(words[i], words[j]):
                return False, (words[i], words[j])
    return True, None


def main():
    ok = True
    P = ((1, 2), (1, 4), (3, 1))
    W = code(P)
    print("9a row [BMRRST1971]: Theta(C_7) >= 343^(1/5) from a 343-word independent set of C_7^(x5)")
    r1, _ = route1(P)
    print("  generator [I_3 | P], P = %s: %d words, %d distinct" % (P, len(W), len(set(W))))
    print("  route 1: every nonzero codeword has a coordinate outside {0, +-1}: %s" % r1)
    r2, _ = route2(W)
    print("  route 2: all %d pairs of words non-adjacent in the strong product: %s" % (len(W) * (len(W) - 1) // 2, r2))
    ok &= r1 and r2 and len(set(W)) == 343
    getcontext().prec = 40
    v = (Decimal(3) * Decimal(7).ln() / 5).exp()
    rounded = v.quantize(Decimal("0.0001"))
    print("  343^(1/5) = 7^(3/5) = %s...; the printed 3.2141 is this value rounded to four decimals (up), printed with "
          "an approximately-equal sign next to the exact 343^(1/5)" % str(v)[:22])
    ok &= str(rounded) == "3.2141"
    print("RESULT [BMRRST1971]: alpha(C_7^(x5)) >= 343, so Theta(C_7) >= 343^(1/5) = 3.2140958...: %s"
          % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g, bad = route1(((1, 2), (1, 4), (1, 1)))
        print("  P with (1,1) for (3,1): %s %s" % ("passes (UNEXPECTED)" if g else "refused, as expected; codeword", bad))
        c_ok &= not g
        g, pair = route2(W + [(1, 1, 0, 0, 0)])
        print("  S plus (1,1,0,0,0): %s %s" % ("independent (UNEXPECTED)" if g else "refused, as expected; pair", pair))
        c_ok &= not g
        claim = len(set(W)) >= 344
        print("  a claim of 344 words from this code: %s" % ("holds (UNEXPECTED)" if claim else "refused, as expected"))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
