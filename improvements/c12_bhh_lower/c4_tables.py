#!/usr/bin/env python3
"""Rigorous tables for the C_4 programs (upper bounds of union areas of unit discs).

  A3[i][j][k], 1 <= i <= j <= k <= 2N: an UPPER bound of the area of the union of three unit discs whose centres have
  pairwise distances i/N, j/N, k/N (sides capped at 2, where the pair no longer overlaps; when k > i + j the triple
  is expanded to the triangle (i, k - i, k), which dominates every side, so by the planar Kneser-Poulsen theorem the
  area can only grow).
  C1[k], 0 <= k <= 2N: an UPPER bound of c(k/N) = pi - lens(k/N), the area of a unit disc outside another one at
  distance k/N (c(d) = pi for d >= 2).

Two independent versions:
  --prog A : integer interval arithmetic, lower_triples_a.union3 / lower_pairs_a (Green's theorem over arcs);
  --prog B : 40-digit decimals, lower_triples_b.union_ie (inclusion-exclusion) + 1e-30.
Each value is written as a decimal string rounded UP at 20 significant digits.

    python3 -I c4_tables.py --prog A --N 64 --out tables_A_64.txt
"""
import math
import os
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
MST = HERE  # the MST lower-bound programs sit next to this file
sys.path.insert(0, MST)


def up20(x):
    """Decimal/Fraction -> decimal string rounded UP at 20 significant digits."""
    getcontext().prec = 60
    d = Decimal(x.numerator) / Decimal(x.denominator) if isinstance(x, Fraction) else Decimal(x)
    getcontext().prec = 20
    getcontext().rounding = ROUND_CEILING
    r = +d
    getcontext().rounding = "ROUND_HALF_EVEN"
    getcontext().prec = 60
    return str(r)


def down20(x):
    getcontext().prec = 60
    d = Decimal(x.numerator) / Decimal(x.denominator) if isinstance(x, Fraction) else Decimal(x)
    getcontext().prec = 20
    getcontext().rounding = "ROUND_FLOOR"
    r = +d
    getcontext().rounding = "ROUND_HALF_EVEN"
    getcontext().prec = 60
    return str(r)


def pairs_low(prog, N, out):
    """L k: a LOWER bound of the union area of two unit discs at distance k/N (2 pi - lens, increasing in d)."""
    M = 2 * N
    lines = [f"N {N}"]
    if prog == "A":
        import lower_triples_a as TA
        from lower_pairs_a import I, PI, sqrt
        for k in range(M + 1):
            if k == 0:
                v = TA.qlo(PI)
            else:
                d = Fraction(k, N)
                lens = 2 * TA.acos_q(d / 2) - I.frac(d / 2) * sqrt(I.frac(4 - d * d))
                v = TA.qlo(PI * 2 - lens)
            lines.append(f"L {k} {down20(v)}")
    else:
        import lower_triples_b as TB
        import lower_pairs_b as LB
        for k in range(M + 1):
            getcontext().prec = TB.PREC_T
            v = (LB.PI if k == 0 else 2 * LB.PI - TB.lens(TB.dec(Fraction(k, N)))) - TB.EPS
            getcontext().prec = 60
            lines.append(f"L {k} {down20(v)}")
    open(out, "w").write("\n".join(lines) + "\n")
    print(f"prog {prog}: wrote {M + 1} lower pair entries to {out}")


def main():
    prog = sys.argv[sys.argv.index("--prog") + 1]
    if "--pairs-low" in sys.argv:
        pairs_low(prog, int(sys.argv[sys.argv.index("--N") + 1]), sys.argv[sys.argv.index("--out") + 1])
        return
    N = int(sys.argv[sys.argv.index("--N") + 1])
    out = sys.argv[sys.argv.index("--out") + 1]
    M = 2 * N
    lines = [f"N {N}"]
    if prog == "A":
        import lower_triples_a as TA
        from lower_pairs_a import I, PI, sqrt
        def c1(k):
            if k >= M:
                return TA.qhi(PI)
            d = Fraction(k, N)
            if k == 0:
                return Fraction(0)
            # c(d) = pi - lens(d), lens(d) = 2 acos(d/2) - (d/2) sqrt(4 - d^2); upper bound of c = pi.hi - lens.lo
            lens = 2 * TA.acos_q(d / 2) - I.frac(d / 2) * sqrt(I.frac(4 - d * d))
            return TA.qhi(PI - lens)
        def a3(i, j, k):
            a, b, c = Fraction(i, N), Fraction(j, N), Fraction(k, N)
            a, b, c = sorted([a, b, c])
            if c > a + b:
                b = c - a
            return TA.qhi(TA.union3(a, b, c))     # union3(a, b, c): discs at (0,0), (a,0), third at distances b, c
    else:
        import lower_triples_b as TB
        import lower_pairs_b as LB
        getcontext().prec = 60
        def c1(k):
            if k >= M:
                return LB.PI + TB.EPS
            if k == 0:
                return Decimal(0)
            getcontext().prec = TB.PREC_T
            v = LB.PI - TB.lens(TB.dec(Fraction(k, N))) + TB.EPS
            getcontext().prec = 60
            return v
        def a3(i, j, k):
            a, b, c = Fraction(i, N), Fraction(j, N), Fraction(k, N)
            a, b, c = sorted([a, b, c])
            if c > a + b:
                b = c - a
            getcontext().prec = TB.PREC_T
            v = TB.union_ie(a, b, c) + TB.EPS
            getcontext().prec = 60
            return v
    for k in range(M + 1):
        lines.append(f"C {k} {up20(c1(k))}")
    cnt = 0
    for i in range(1, M + 1):
        for j in range(i, M + 1):
            for k in range(j, M + 1):
                lines.append(f"A {i} {j} {k} {up20(a3(i, j, k))}")
                cnt += 1
        if i % 16 == 0:
            print(f"prog {prog}: i = {i}/{M}, {cnt} entries", flush=True)
    open(out, "w").write("\n".join(lines) + "\n")
    print(f"prog {prog}: wrote {cnt} triple entries and {M + 1} pair entries to {out}")


if __name__ == "__main__":
    main()
