#!/usr/bin/env python3
"""Independent replay of the four entropy lower bounds for C_3c (the 3c page: [G2026], [MI2026], [L2026], [J2026]),
from the laws alone. No program of any certificate was read or run; the laws are read as data.

The entropy formulation on the 3c page: C_3c is the least constant with
    H(X - Y) <= C_3c max(H(X), H(Y), H(X + Y), H(X + 2Y))
for every pair of discrete random variables, so every explicit law of (X, Y) gives C_3c >= that ratio. For a law we
form the five pushforwards exactly (rationals), and enclose H(F) = -sum p ln p rigorously: decimal at 90 digits,
each logarithm correctly rounded (ln p = ln a - ln b for p = a/b), the sum widened by 10^-80 (far above the rounding
of a few hundred operations at 90 digits). [J2026] is a product of independent base-55 digits, n1 from law 1 and n2
from law 2: when every one of the five forms of every support point is at most 27 in absolute value, each form of
(X, Y) = sum_i 55^i (x_i, y_i) has the balanced base-55 digits of its digit-forms, so it determines them, and
H(F) = n1 H_1(F) + n2 H_2(F) for all five forms; we check that condition, and x - y injective on each support (the
row's statement), and evaluate the ratio of the two sums.

  python3 -I c3c_check.py ROW [CERTS] [--controls]       rows: g2026, mi2026, l2026, j2026
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

# [G2026]: the table of PR #70 (x, y, rounded weight); the PR normalises the weights, which sum to
# 0.99999999999852250171
G2026 = [(1, 7, "1.30506564e-12"), (3, 5, "9.25745153e-7"), (3, 6, "1.52498825e-6"), (3, 7, "6.57289640e-8"),
         (5, 4, "0.001244422280"), (5, 5, "0.006113255152"), (5, 6, "0.000037354259"), (7, 2, "0.016223522296"),
         (7, 3, "0.169721258229"), (7, 4, "0.008420798270"), (7, 5, "0.011728626004"), (9, 1, "0.248375293086"),
         (9, 2, "0.046862624648"), (9, 3, "0.243496670072"), (11, 0, "0.010100758400"), (11, 1, "0.188770387301"),
         (11, 2, "0.018228788494"), (13, -1, "0.006653795672"), (13, 0, "0.002961930500"),
         (13, 1, "0.020758390703"), (15, -2, "0.000003025952"), (15, -1, "0.000275897370"),
         (15, 0, "0.000020680540"), (17, -3, "6.99564607e-12"), (17, -2, "1.51075546e-9"),
         (17, -1, "2.78909933e-9")]

ROWS = {
    "g2026": dict(claim="1.67473389", where="3c row [G2026] (Griego, PR #70): 26 points, the table in the PR"),
    "mi2026": dict(claim="1.6747338950208249", file=("mi2026", "certificate_3c_v3_95pt.json"),
                   where="3c row [MI2026] (Mosaic Intelligence, PR #93): 95 points, Zenodo 10.5281/zenodo.20794135"),
    "l2026": dict(claim="1.6747338950414058", file=("l2026", "certificate_3c_147pt.json"),
                  long="1.674733895041405870063135756722213999136383713818148696811828",
                  where="3c row [L2026] (Lin, PR #185): 147 points at denominator 10^320, gist revision 62621d2"),
    "j2026": dict(claim="1.674894211307116", file=("j2026", "certificate.json"),
                  where="3c row [J2026] (Jaswal, PR #221): two laws on 51 and 47 points, base 55, n = 10^20 digits"),
}
FORMS = (("X", lambda x, y: x), ("Y", lambda x, y: y), ("X+Y", lambda x, y: x + y), ("X+2Y", lambda x, y: x + 2 * y),
         ("X-Y", lambda x, y: x - y))


def load_law(key, certs):
    """[(x, y, Fraction p)] with sum p = 1"""
    if key == "g2026":
        w = [(x, y, Fr(s)) for x, y, s in G2026]
        tot = sum(p for _, _, p in w)
        assert tot == Fr("0.99999999999852250171"), tot
        return [(x, y, p / tot) for x, y, p in w]
    d = json.load(open(os.path.join(certs, *ROWS[key]["file"]), encoding="utf-8"))
    if key == "mi2026":
        law = [(x, y, Fr(nd[0], nd[1])) for x, y, nd in d["distribution_x_y_numden"]]
    elif key == "l2026":
        assert d["denominator"] == "10^320"
        law = [(x, y, Fr(n, 10 ** 320)) for x, y, n in d["distribution_x_y_num"]]
    else:
        raise ValueError(key)
    return law


def check_law(law):
    pts = [(x, y) for x, y, _ in law]
    return (len(set(pts)) == len(pts) and all(p > 0 for _, _, p in law) and sum(p for _, _, p in law) == 1)


def push(law, f):
    out = {}
    for x, y, p in law:
        k = f(x, y)
        out[k] = out.get(k, 0) + p
    return out


def entropy(dist):
    """[lo, hi] for -sum p ln p (natural logarithm)"""
    s = D(0)
    for p in dist.values():
        if p == 0:
            continue
        lp = D(p.numerator).ln() - D(p.denominator).ln()
        s -= (D(p.numerator) / D(p.denominator)) * lp
    return s - WIDEN, s + WIDEN


def entropies(law):
    return {name: entropy(push(law, f)) for name, f in FORMS}


def ratio(H):
    """[lo, hi] of H(X-Y) / max of the other four"""
    mx_lo = max(H[n][0] for n in ("X", "Y", "X+Y", "X+2Y"))
    mx_hi = max(H[n][1] for n in ("X", "Y", "X+Y", "X+2Y"))
    return H["X-Y"][0] / mx_hi, H["X-Y"][1] / mx_lo


def floor_dec(x, n):
    q = D(10) ** -n
    return x.quantize(q, rounding="ROUND_FLOOR")


def single(key, certs, controls, log):
    row = ROWS[key]
    law = load_law(key, certs)
    good = check_law(law)
    log(f"  {len(law)} support points, distinct, every weight > 0, total exactly 1: {good}")
    ok = good
    inj = len({x - y for x, y, _ in law}) == len(law)
    log(f"  x - y injective on the support (so H(X - Y) = H(X, Y)): {inj}")
    H = entropies(law)
    for name, _ in FORMS:
        log(f"  H({name}) in [{floor_dec(H[name][0], 30)}, {floor_dec(H[name][1], 30)}]")
    lo, hi = ratio(H)
    log(f"  ratio in [{floor_dec(lo, 40)}, {floor_dec(hi, 40)}]")
    claim = D(row["claim"])
    good = lo >= claim
    log(f"  the row's {row['claim']} <= ratio: {good} (margin {float(lo - claim):.3e})")
    ok &= good
    if "long" in row:
        good = lo >= D(row["long"])
        log(f"  the 60-digit certified value {row['long']} <= ratio: {good}")
        ok &= good
    if controls:
        log("Red controls (each must fall below the row's value):")
        heavy = max(range(len(law)), key=lambda i: law[i][2])
        light = min(range(len(law)), key=lambda i: law[i][2])
        moved = [list(t) for t in law]
        m = moved[heavy][2] / 100
        moved[heavy][2] -= m
        moved[light][2] += m
        lo2, _ = ratio(entropies([tuple(t) for t in moved]))
        hi2 = ratio(entropies([tuple(t) for t in moved]))[1]
        good = hi2 < claim
        log(f"  1% of the heaviest point's mass moved to the lightest: ratio <= {floor_dec(hi2, 12)}: "
            f"{'below, as expected' if good else 'NOT BELOW'}")
        ok &= good
        shifted = [(x, y + 1 if i == heavy else y, p) for i, (x, y, p) in enumerate(law)]
        if len({(x, y) for x, y, _ in shifted}) == len(shifted):
            hi3 = ratio(entropies(shifted))[1]
            good = hi3 < claim
            log(f"  the heaviest point moved by (0, 1): ratio <= {floor_dec(hi3, 12)}: "
                f"{'below, as expected' if good else 'NOT BELOW'}")
            ok &= good
    return ok


def two_law(certs, controls, log):
    row = ROWS["j2026"]
    d = json.load(open(os.path.join(certs, *row["file"]), encoding="utf-8"))
    base, counts = d["base"], [int(c) for c in d["counts"]]
    claim = Fr(*d["claimed_lower_bound"])
    ok = base == 55 and sum(counts) == 10 ** 20 and claim == Fr(row["claim"])
    log(f"  base {base}, counts {counts[0]} + {counts[1]} = 10^20, claimed {claim} = {row['claim']}: {ok}")
    laws = []
    for L in d["laws"]:
        law = [(x, y, Fr(n, den)) for (x, y), (n, den) in zip(L["points"], L["weights"])]
        laws.append(law)
    Hs = []
    for i, law in enumerate(laws):
        good = check_law(law)
        inj = len({x - y for x, y, _ in law}) == len(law)
        small = max(abs(f(x, y)) for x, y, _ in law for _, f in FORMS)
        log(f"  law {i + 1}: {len(law)} points, distinct, weights > 0 summing to 1: {good}; x - y injective: {inj}; "
            f"largest |form| = {small} <= 27: {small <= 27}")
        ok &= good and inj and small <= 27
        Hs.append(entropies(law))
    n1, n2 = counts
    tot = n1 + n2
    H = {name: ((n1 * Hs[0][name][0] + n2 * Hs[1][name][0]) / tot, (n1 * Hs[0][name][1] + n2 * Hs[1][name][1]) / tot)
         for name, _ in FORMS}
    for name, _ in FORMS:
        log(f"  per digit, H({name}) in [{floor_dec(H[name][0], 25)}, {floor_dec(H[name][1], 25)}]")
    lo, hi = ratio(H)
    log(f"  ratio in [{floor_dec(lo, 30)}, {floor_dec(hi, 30)}]")
    good = lo >= D(claim.numerator) / D(claim.denominator)
    log(f"  the row's {row['claim']} <= ratio: {good} (margin {float(lo - D(row['claim'])):.3e})")
    ok &= good
    if controls:
        log("Red controls (each must fall below the row's value):")
        c = D(row["claim"])
        for what, (a, b) in (("law 1 alone", (tot, 0)), ("law 2 alone", (0, tot)), ("the counts swapped", (n2, n1))):
            Hc = {name: ((a * Hs[0][name][0] + b * Hs[1][name][0]) / tot, (a * Hs[0][name][1] + b * Hs[1][name][1]) / tot)
                  for name, _ in FORMS}
            h = ratio(Hc)[1]
            good = h < c
            log(f"  {what}: ratio <= {floor_dec(h, 12)}: {'below, as expected' if good else 'NOT BELOW'}")
            ok &= good
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args or args[0] not in ROWS:
        print(__doc__)
        return 2
    key = args[0]
    certs = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    t0 = time.time()
    print(ROWS[key]["where"] + f": C_3c >= {ROWS[key]['claim']}")
    if key == "j2026":
        ok = two_law(certs, "--controls" in sys.argv, print)
    else:
        ok = single(key, certs, "--controls" in sys.argv, print)
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
