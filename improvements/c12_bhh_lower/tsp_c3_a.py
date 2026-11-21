#!/usr/bin/env python3
"""Certified lower bound for Dong and Cao's C_3 (arXiv 2610.02809, Theorem 3.11), program A.

Dong-Cao: beta_TSP >= 5/8 + sum_{k=3}^K C_k,
  C_k = Gamma(k - 1/2)/(2 k!) int_{(R^2)^{k-1}} g(G_1(y)) A(y)^{-(k-1/2)} dy_2 ... dy_k   (y_1 = 0),
with g(G) = q(G) - D(G) for connected G (0 otherwise) and A(y) the area of the union of the unit discs at y_1..y_k.
For k = 3 a connected 3-vertex graph is a path (q = 1, D = 1, g = 0) or a triangle (q = 1, D = 0, g = 1), so
  C_3 = (1/6) int 1[|y2|, |y3|, |y2 - y3| <= 1] G(A(0, y2, y3)) dy2 dy3,   G(A) = int_0^inf r^4 e^{-A r^2} dr
      = 3 sqrt(pi) / (8 A^{5/2}),
and in the coordinates a = |y2|, b = |y3|, theta = the angle between them (dy2 dy3 = 4 pi a b da db dtheta over
theta in [0, pi]):
  C_3 = 4 pi int_0^1 int_0^1 int_0^pi w G(A(a, b, c)) a b dtheta db da,  w = 1/6 if c <= 1 else 0,
  c^2 = a^2 + b^2 - 2 a b cos(theta).
This is exactly the triple term of lower_triples_a.py (our MST lower-bound program, in this folder) with the cut a = 0 (so G3(0, A) = G(A))
and with path cells weighted 0 instead of 1/2.  A cell gets weight 1/6 only if its LARGEST c is <= 1 (every
configuration in it is a triangle), otherwise 0: a lower bound.  Monotonicity of A in the three distances: planar
Kneser-Poulsen (Bezdek and Connelly 2002); the union area by Green's theorem in integer interval arithmetic, from
lower_triples_a.py (imported, unchanged).

Usage: python3 -I tsp_c3_a.py [--na 128] [--nt 256] [--nc 128] [--claim Q] [--controls]
"""
import math
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
MST = HERE
sys.path.insert(0, MST)  # the MST lower-bound programs sit next to this file
import lower_triples_a as TA  # noqa: E402
from lower_pairs_a import I, PI, sqrt  # noqa: E402

# Dong-Cao's certified values (arXiv 2610.02809v1, Section 3.3.1), used as cited, not recomputed here
DC_C3 = Fraction("0.01113528859")
DC_C4 = Fraction("0.005040573276")
DC_C5 = Fraction("0.001015487669")


def G0(A):
    """G(A) = int_0^inf r^4 e^{-A r^2} dr = 3 sqrt(pi) / (8 A^{5/2}) for an interval A > 0 (interval result).
    (lower_triples_a.G3 with cut 0 is not used: its erf/exp series at argument 0 does not terminate quickly.)"""
    return 3 * sqrt(PI) / (8 * A * A * sqrt(A))


def c3_term(na, nt, nc, log=True, wtri=Fraction(1, 6), wmixed=Fraction(0)):
    acut = Fraction(0)
    coses = [TA.cos_int(PI * Fraction(m, nt)) for m in range(nt + 1)]
    table = {}
    total = Fraction(0)
    dth = PI * Fraction(1, nt)
    four_pi = PI * 4
    t0 = time.time()
    for i in range(na):
        a0, a1 = Fraction(i, na), Fraction(i + 1, na)
        for j in range(na):
            b0, b1 = Fraction(j, na), Fraction(j + 1, na)
            wab = (a1 * a1 - a0 * a0) / 2 * (b1 * b1 - b0 * b0) / 2
            corners = [(a0, b0), (a0, b1), (a1, b0), (a1, b1)]
            acc = I(0)
            for m in range(nt):
                k1lo = TA.qlo(coses[m + 1])          # cos(theta) >= k1lo on the cell
                c2max = max(p * p + q * q - 2 * k1lo * p * q for p, q in corners)
                w = wtri if c2max <= 1 else wmixed
                if w == 0:
                    continue
                X = c2max * nc * nc
                k = math.isqrt(X.numerator // X.denominator)
                while k * k < X:
                    k += 1
                k = max(k, 1)
                key = (min(i + 1, j + 1), max(i + 1, j + 1), k)
                if key not in table:
                    A_ = TA.union3(*TA.realise(Fraction(key[0], na), Fraction(key[1], na), Fraction(k, nc)))
                    table[key] = TA.qlo(G0(I(A_.hi)))   # G at the upper end of A: a lower bound
                acc = acc + I.frac(table[key] * w)
            total += TA.qlo(four_pi * dth * acc * wab)
        if log and (i % 16 == 15 or i == na - 1):
            print(f"  rows a < {float(a1):.4f}: C_3 >= {float(total):.12f} ({len(table)} union areas, "
                  f"{time.time() - t0:.0f} s)", flush=True)
    return total, len(table)


def main(argv):
    na, nt, nc = 128, 256, 128
    claim = None
    controls = False
    i = 0
    while i < len(argv):
        k = argv[i]
        if k == "--controls":
            controls = True; i += 1; continue
        v = argv[i + 1]; i += 2
        if k == "--na": na = int(v)
        elif k == "--nt": nt = int(v)
        elif k == "--nc": nc = int(v)
        elif k == "--claim": claim = Fraction(v)
        else: raise SystemExit("unknown argument " + k)
    # self-test: G(A) with acut = 0 equals 3 sqrt(pi)/(8 A^{5/2}) (A = 4: 3 sqrt(pi)/256)
    g = G0(I.frac(4))
    ref = 3.0 * math.sqrt(math.pi) / 256
    assert abs(float(TA.qlo(g)) - ref) < 1e-15 and TA.qlo(g) <= Fraction(ref) * (1 + Fraction(1, 10 ** 12))
    print("self-test: G(4) = 3 sqrt(pi)/256: ok")
    C3, nu = c3_term(na, nt, nc)
    beta = Fraction(5, 8) + C3 + DC_C4 + DC_C5
    print(f"C_3 >= {float(C3):.12f}  (na {na}, nt {nt}, nc {nc}; {nu} union areas)  [Dong-Cao: {float(DC_C3)}]")
    print(f"exact rational lower bound for C_3: {C3.numerator}/{C3.denominator}" if C3.denominator < 10 ** 400 else
          f"C_3 lower bound as a decimal floor: {math.floor(C3 * 10 ** 15)}e-15")
    print(f"with Dong-Cao's C_4 >= {float(DC_C4)} and C_5 >= {float(DC_C5)}: beta_TSP >= 5/8 + C_3 + C_4 + C_5 >= {float(beta):.12f}")
    print(f"Dong-Cao's own total: {float(Fraction(5, 8) + DC_C3 + DC_C4 + DC_C5):.12f}")
    ok = True
    if claim is not None:
        ok = beta >= claim
        print(f"claim beta_TSP >= {claim}: {'HOLDS' if ok else 'NOT IMPLIED'}")
    if controls:
        bad = 0
        if claim is not None:
            hi = Fraction(int(beta * 10 ** 7) + 1, 10 ** 7)
            r = not (beta >= hi); bad += not r
            print(f"red control 1 (claim {float(hi)} one unit above the 7th decimal): {'refused' if r else 'ACCEPTED (bad)'}")
        Cc, _ = c3_term(8, 16, 8, log=False)
        r = Cc < C3; bad += not r
        print(f"red control 2 (coarse grid 8/16/8 gives C_3 >= {float(Cc):.9f} < fine): {'refused' if r else 'ACCEPTED (bad)'}")
        Cm, _ = c3_term(na // 4, nt // 4, nc // 4, log=False, wmixed=Fraction(1, 6))
        Cn, _ = c3_term(na // 4, nt // 4, nc // 4, log=False)
        r = Cm > Cn; bad += not r
        print(f"red control 3 (mixed triangle/path cells wrongly counted as triangles raises the value: "
              f"{float(Cm):.9f} > {float(Cn):.9f}): {'detected' if r else 'NOT DETECTED (bad)'}")
        if bad:
            print("RED CONTROLS FAILED"); return 2
        print("red controls: all behaved as required")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
