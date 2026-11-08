#!/usr/bin/env python3
"""Independent replay of the 3b lower-bound row [MI2026], C_3b >= 1.77898884 (PR #92), from the 13-point law alone
(c3b_pr92/certificate_3b_13pt.json, extracted by itself from Zenodo 10.5281/zenodo.20794135; the archive's programs
are neither read nor run).

The entropy formulation on the 3b page: C_3b is the least constant with H(X - Y) <= C_3b max(H(X), H(Y), H(X + Y))
for every pair of discrete random variables, so the law gives C_3b >= that ratio. Exact pushforwards; each entropy in
90-digit decimal with correctly rounded logarithms, widened by 1e-80 (as in our c3c_check.py).
  python3 -I mi2026_check.py [CERTS] [--controls]
Python 3.9+, standard library only."""
import json
import os
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction as Fr

getcontext().prec = 90
WIDEN = D("1e-80")
HERE = os.path.dirname(os.path.abspath(__file__))
CLAIM = D("1.77898884")
FORMS = (("X", lambda x, y: x), ("Y", lambda x, y: y), ("X+Y", lambda x, y: x + y), ("X-Y", lambda x, y: x - y))


def entropy(law, f):
    dist = {}
    for x, y, p in law:
        dist[f(x, y)] = dist.get(f(x, y), 0) + p
    s = D(0)
    for p in dist.values():
        s -= (D(p.numerator) / D(p.denominator)) * (D(p.numerator).ln() - D(p.denominator).ln())
    return s - WIDEN, s + WIDEN


def ratio(law):
    H = {n: entropy(law, f) for n, f in FORMS}
    mx_lo = max(H[n][0] for n in ("X", "Y", "X+Y"))
    mx_hi = max(H[n][1] for n in ("X", "Y", "X+Y"))
    return H, H["X-Y"][0] / mx_hi, H["X-Y"][1] / mx_lo


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    certs = args[0] if args else os.path.join(HERE, "certs")
    t0 = time.time()
    d = json.load(open(os.path.join(certs, "mi2026", "certificate_3b_13pt.json"), encoding="utf-8"))
    law = [(w["x"], w["y"], Fr(int(w["num"]), int(w["den"]))) for w in d["weights"]]
    print("3b row [MI2026] (PR #92): C_3b >= 1.77898884, the 13-point law")
    ok = len(law) == 13 and len({(x, y) for x, y, _ in law}) == 13 and all(p > 0 for *_, p in law) and \
        sum(p for *_, p in law) == 1
    print(f"  13 distinct points, positive weights summing to 1 exactly: {ok}")
    sym = {(x, y): p for x, y, p in law} == {(y, x): p for x, y, p in law}
    print(f"  the law is symmetric under (x, y) -> (y, x), as the certificate states: {sym}")
    ok &= sym
    H, lo, hi = ratio(law)
    for n, _ in FORMS:
        print(f"  H({n}) in [{H[n][0]:.30f}, {H[n][1]:.30f}]")
    print(f"  ratio in [{lo:.40f}, {hi:.40f}]")
    good = lo >= CLAIM and str(lo).startswith("1.778988841420693549")
    print(f"  the row's 1.77898884 <= ratio, and the certificate's prefix 1.778988841420693549 matches: {good}")
    ok &= good
    if "--controls" in sys.argv:
        print("Red controls (each must fall below the row's value):")
        heavy = max(range(13), key=lambda i: law[i][2])
        light = min(range(13), key=lambda i: law[i][2])
        moved = [list(t) for t in law]
        m = moved[heavy][2] / 100
        moved[heavy][2] -= m
        moved[light][2] += m
        h2 = ratio([tuple(t) for t in moved])[2]
        good = h2 < CLAIM
        print(f"  1% of the heaviest point's mass moved to the lightest: ratio <= {h2:.12f}: "
              f"{'below, as expected' if good else 'NOT BELOW'}")
        ok &= good
        uni = [(x, y, Fr(1, 13)) for x, y, _ in law]
        h3 = ratio(uni)[2]
        good = h3 < CLAIM
        print(f"  the uniform law on the same 13 points: ratio <= {h3:.12f}: {'below, as expected' if good else 'NOT BELOW'}")
        ok &= good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
