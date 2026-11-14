#!/usr/bin/env python3
"""90a, program B: the same covering system by inclusion-exclusion, in exact rationals.

A finite union of residue classes is periodic, so it is all of Z exactly when its density is 1. The density of an
intersection of classes is 0 when two of them are incompatible (a_i != a_j modulo gcd(n_i, n_j)) and 1/lcm otherwise
(Chinese remainder theorem), so

    density of the union = sum over nonempty subsets S of (-1)^(|S| + 1) [S compatible] / lcm(S).

The program evaluates this sum with fractions.Fraction over all 2^k - 1 subsets and accepts only the value 1. No
residue is enumerated, so the route is independent of program A. Moduli must again be distinct and greater than 1.

Standard library only (python3 -I).  usage: c90a_density.py "a1:n1 a2:n2 ..."
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd


def compatible(sub):
    for (a, n), (b, m) in combinations(sub, 2):
        if (a - b) % gcd(n, m):
            return False
    return True


def lcm_of(sub):
    out = 1
    for _, n in sub:
        out = out // gcd(out, n) * n
    return out


def main(argv):
    system = []
    for part in argv[1].split():
        a, n = part.split(":")
        system.append((int(a), int(n)))
    mods = sorted(n for _, n in system)
    fine = mods[0] > 1 and all(mods[i] != mods[i + 1] for i in range(len(mods) - 1))
    print("B: sorted moduli %s, distinct and greater than 1: %s" % (mods, fine))
    total = Fraction(0)
    for size in range(1, len(system) + 1):
        for sub in combinations(system, size):
            if compatible(sub):
                total += Fraction((-1) ** (size + 1), lcm_of(sub))
    print("B: density of the union = %s" % total)
    print("B: least modulus %d" % mods[0])
    print("B: VERDICT %s" % ("PASS" if fine and total == 1 else "REFUSED"))


if __name__ == "__main__":
    main(sys.argv)
