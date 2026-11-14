#!/usr/bin/env python3
"""90a, program A: the page's elementary distinct covering system, checked residue by residue.

A system of classes a_i mod n_i covers Z exactly when it covers one full period, the residues 0 .. L - 1 with
L = lcm(n_i). The program checks that the moduli are integers greater than 1 and pairwise distinct, reports the least
modulus, and lists every residue of the period that no class contains.

Standard library only (python3 -I).  usage: c90a_cover.py "a1:n1 a2:n2 ..."
"""
import sys
from math import gcd


def main(argv):
    classes = [tuple(int(x) for x in tok.split(":")) for tok in argv[1].split()]
    moduli = [n for _, n in classes]
    ok = True
    proper = all(n > 1 for n in moduli)
    distinct = len(set(moduli)) == len(moduli)
    print("A: %d classes, moduli %s, all > 1: %s, pairwise distinct: %s" % (len(classes), moduli, proper, distinct))
    ok = ok and proper and distinct
    period = 1
    for n in moduli:
        period = period * n // gcd(period, n)
    missed = [r for r in range(period) if not any((r - a) % n == 0 for a, n in classes)]
    print("A: period %d, residues not covered: %s" % (period, missed if missed else "none"))
    ok = ok and not missed
    print("A: least modulus %d" % min(moduli))
    print("A: VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main(sys.argv)
