#!/usr/bin/env python3
"""Certified lower bound for Dong and Cao's C_3 (arXiv 2610.02809, Theorem 3.11), program B.

Shares no code with tsp_c3_a.py: it uses lower_triples_b.py (our MST lower-bound program, in this folder) (Decimal arithmetic, union area by
inclusion-exclusion with the triple intersection from arcs), which shares no code with lower_triples_a.py.
Same mathematics (see tsp_c3_a.py): C_3 = 4 pi int_0^1 int_0^1 int_0^pi w G(A(a,b,c)) a b dtheta db da with
G(A) = 3 sqrt(pi)/(8 A^{5/2}) (= G3 with cut 0), w = 1/6 on cells whose largest c is <= 1, 0 otherwise.
Rounding: as lower_triples_b.py (40 digits; union area raised by 1e-30, each G value and each cell lowered by 1e-30).

Usage: python3 -I tsp_c3_b.py [--na 120] [--nt 300] [--nc 160] [--claim Q] [--controls]
"""
import math
import os
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
MST = HERE
sys.path.insert(0, MST)  # the MST lower-bound programs sit next to this file
import lower_triples_b as TB  # noqa: E402
import lower_pairs_b as LB  # noqa: E402

DC_C3 = D("0.01113528859")
DC_C4 = D("0.005040573276")
DC_C5 = D("0.001015487669")


def c3_term(na, nt, nc, log=True, wmixed=Fraction(0)):
    cb = [TB.cos_bounds(m, nt) for m in range(nt + 1)]
    getcontext().prec = TB.PREC_T
    tab = {}
    total = D(0)
    pi = LB.PI
    t0 = time.time()
    for i in range(na):
        a0, a1 = Fraction(i, na), Fraction(i + 1, na)
        for j in range(na):
            b0, b1 = Fraction(j, na), Fraction(j + 1, na)
            box = [(a0, b0), (a0, b1), (a1, b0), (a1, b1)]
            cell_sum = D(0)
            for m in range(nt):
                lo1 = cb[m + 1][0]                      # cos(theta) >= lo1 on the cell
                c2hi = max(p * p + q * q - 2 * lo1 * p * q for p, q in box)
                w = Fraction(1, 6) if c2hi <= 1 else wmixed
                if w == 0:
                    continue
                Xn = c2hi * nc * nc
                k = math.isqrt(Xn.numerator // Xn.denominator)
                while Fraction(k * k) < Xn:
                    k += 1
                k = max(k, 1)
                aa, bb = (i + 1, j + 1) if i <= j else (j + 1, i + 1)
                key = (aa, bb, k)
                g = tab.get(key)
                if g is None:
                    A_, B_, C_ = Fraction(aa, na), Fraction(bb, na), Fraction(k, nc)
                    C_ = max(C_, abs(A_ - B_))
                    if C_ > A_ + B_:
                        B_ = C_ - A_
                    U = TB.union_ie(A_, B_, C_) + TB.EPS
                    g = TB.G3(Fraction(0), U) - TB.EPS
                    tab[key] = g
                cell_sum += g * TB.dec(w)
            wab = TB.dec((a1 * a1 - a0 * a0) / 2 * (b1 * b1 - b0 * b0) / 2)
            total += 4 * pi * (pi / nt) * wab * cell_sum - TB.EPS
        if log and (i % 20 == 19 or i == na - 1):
            print(f"  a < {float(a1):.4f}: C_3 >= {total:.12f} ({len(tab)} union areas, {time.time() - t0:.0f} s)",
                  flush=True)
    getcontext().prec = 60
    return total, len(tab)


def main(argv):
    na, nt, nc = 120, 300, 160
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
        elif k == "--claim": claim = D(v)
        else: raise SystemExit("unknown argument " + k)
    getcontext().prec = 60
    g = TB.G3(Fraction(0), D(4))
    assert abs(float(g) - 3 * math.sqrt(math.pi) / 256) < 1e-15
    u = TB.union_ie(Fraction(1), Fraction(1), Fraction(1))
    L1 = 2 * math.pi / 3 - math.sqrt(3) / 2
    assert abs(float(u) - (3 * math.pi - 3 * L1 + (math.pi - math.sqrt(3)) / 2)) < 1e-12
    print("self-tests: G(4) = 3 sqrt(pi)/256, union of three unit discs (equilateral side 1): ok")
    C3, nu = c3_term(na, nt, nc)
    beta = D(5) / D(8) + C3 + DC_C4 + DC_C5
    print(f"C_3 >= {C3:.15f}  (na {na}, nt {nt}, nc {nc}; {nu} union areas)  [Dong-Cao: {DC_C3}]")
    print(f"with Dong-Cao's C_4 >= {DC_C4} and C_5 >= {DC_C5}: beta_TSP >= {beta:.12f}")
    ok = True
    if claim is not None:
        ok = beta >= claim
        print(f"claim beta_TSP >= {claim}: {'HOLDS' if ok else 'NOT IMPLIED'}")
    if controls:
        bad = 0
        if claim is not None:
            hi = (beta * 10 ** 7).to_integral_value(rounding="ROUND_FLOOR") / D(10 ** 7) + D("1e-7")
            r = not (beta >= hi); bad += not r
            print(f"red control 1 (claim {hi}): {'refused' if r else 'ACCEPTED (bad)'}")
        Cc, _ = c3_term(8, 16, 8, log=False)
        r = Cc < C3; bad += not r
        print(f"red control 2 (coarse grid gives {Cc:.9f} < fine): {'refused' if r else 'ACCEPTED (bad)'}")
        Cm, _ = c3_term(max(na // 4, 4), max(nt // 4, 4), max(nc // 4, 4), log=False, wmixed=Fraction(1, 6))
        Cn, _ = c3_term(max(na // 4, 4), max(nt // 4, 4), max(nc // 4, 4), log=False)
        r = Cm > Cn; bad += not r
        print(f"red control 3 (mixed cells counted as triangles: {Cm:.9f} > {Cn:.9f}): {'detected' if r else 'NOT DETECTED (bad)'}")
        if bad:
            print("RED CONTROLS FAILED"); return 2
        print("red controls: all behaved as required")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
