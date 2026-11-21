#!/usr/bin/env python3
"""MST constant in the plane: certified lower bound with the isolated-pair term (program B).

Written separately from lower_pairs_a.py: it shares no code with it, uses the decimal module at 60 significant
digits, parametrises the pair by the half-angle theta of the lens instead of the distance, uses the positive series
for erf, the BBP series for pi, and Taylor series for sin and cos.

Mathematics (the same bound as program A):
  beta >= int_0^a chi(r) dr + int_a^inf e^{-pi r^2} dr + int_a^inf n2(r) dr,
  n2(r) = pi int_0^r t e^{-|B(0,r) u B(t e1, r)|} dt   (isolated pairs per point).
Two unit discs at distance u = 2 cos(theta), theta in [pi/3, pi/2], have union area B = 2 pi - (2 theta - sin 2theta)
and u du = -2 sin(2 theta) d theta, so
  int_a^inf n2 = pi int_{pi/3}^{pi/2} 2 sin(2 theta) G(B(theta)) d theta,   G(B) = int_a^inf r^2 e^{-B r^2} dr.
B decreases in theta, G decreases in B, so G(B(theta)) increases in theta; on a cell [t0, t1] the integral is at
least G(B(t0)) (cos 2t0 - cos 2t1).
int_0^a chi = erf(a sqrt(pi)/2)/2 + (a/2) e^{-pi a^2/4};  int_a^inf e^{-pi r^2} = erfc(a sqrt(pi))/2;
G(B) = a e^{-B a^2}/(2B) + sqrt(pi) erfc(a sqrt(B))/(4 B^{3/2}).

Rounding: 60 significant digits, round-half-even; exp and sqrt of the decimal module are correctly rounded. Every
series is summed until its tail bound is below 1e-58 and the tail bound is added in the safe direction. With at most a
few thousand operations per cell on numbers below 100, the accumulated rounding error per cell is below 1e-52; each
cell's contribution is lowered by 1e-45 and each closed-form term by 1e-45, far more than that.

Usage: python3 -I lower_pairs_b.py [--a 2679/2500] [--cells 150000] [--claim 0.632993] [--controls]
"""
import sys
from decimal import Decimal as D, getcontext

getcontext().prec = 60
MARGIN = D("1e-45")
TINY = D("1e-58")


def pi_bbp():
    s = D(0)
    k = 0
    while True:
        t = (D(4) / (8 * k + 1) - D(2) / (8 * k + 4) - D(1) / (8 * k + 5) - D(1) / (8 * k + 6)) / (D(16) ** k)
        s += t
        if t < TINY:
            break
        k += 1
    return s  # tail after the last term is below the last term / 15; absorbed by MARGIN


PI = pi_bbp()
SQPI = PI.sqrt()


def erf(z):
    """erf(z) = 2/sqrt(pi) e^{-z^2} sum_{n>=0} 2^n z^{2n+1}/(1*3*...*(2n+1)), z >= 0; positive terms with ratio
    2 z^2/(2n+3), so once that ratio is below 1/2 the tail is at most the last term."""
    z2 = z * z
    t = z
    s = t
    n = 0
    while True:
        n += 1
        t = t * 2 * z2 / (2 * n + 1)
        s += t
        if 2 * z2 / (2 * n + 3) < D("0.5") and t < TINY:
            break
    return 2 / SQPI * (-z2).exp() * s


def sin_cos(x):
    """sin and cos for 0 <= x <= 4 by Taylor series (alternating, terms decreasing once n > x)."""
    s = D(0); c = D(0)
    t = D(1)  # x^n/n!
    n = 0
    while True:
        if n % 4 == 0: c += t
        elif n % 4 == 1: s += t
        elif n % 4 == 2: c -= t
        else: s -= t
        n += 1
        t = t * x / n
        if n > x + 2 and t < TINY:
            break
    return s, c


def G(a, B):
    sB = B.sqrt()
    return a * (-(B * a * a)).exp() / (2 * B) + SQPI * (1 - erf(a * sB)) / (4 * B * sB)


def bound(a, cells):
    a = D(a.split("/")[0]) / D(a.split("/")[1]) if "/" in a else D(a)
    head = erf(a * SQPI / 2) / 2 + a / 2 * (-(PI * a * a / 4)).exp() - MARGIN
    t1 = (1 - erf(a * SQPI)) / 2 - MARGIN
    total = D(0)
    h = PI / 6 / cells
    th0 = PI / 3
    s0, c0 = sin_cos(2 * th0)
    for i in range(cells):
        th1 = PI / 3 + h * (i + 1)
        s1, c1 = sin_cos(2 * th1)
        B0 = 2 * PI - (2 * th0 - s0)  # union area at the cell's left end (the largest on the cell)
        total += PI * G(a, B0) * (c0 - c1) - MARGIN
        th0, s0, c0 = th1, s1, c1
    return head + t1 + total, head, t1, total


def main(argv):
    a = "2679/2500"; cells = 150000; claim = D("0.632993"); controls = False
    i = 0
    while i < len(argv):
        if argv[i] == "--a": a = argv[i + 1]; i += 2
        elif argv[i] == "--cells": cells = int(argv[i + 1]); i += 2
        elif argv[i] == "--claim": claim = D(argv[i + 1]); i += 2
        elif argv[i] == "--controls": controls = True; i += 1
        else: raise SystemExit("unknown argument " + argv[i])
    # self-tests
    assert abs(PI - D("3.14159265358979323846264338327950288419716939937510")) < D("1e-50")
    assert abs(erf(D(1)) - D("0.842700792949714869341220635082609259296066997966302908459937")) < D("1e-50")
    s, c = sin_cos(PI / 3)
    assert abs(s * s + c * c - 1) < D("1e-50") and abs(c - D("0.5")) < D("1e-50")
    B1 = 2 * PI - (2 * PI / 3 - sin_cos(2 * PI / 3)[0])  # distance 1: theta = pi/3
    assert abs(B1 - (4 * PI / 3 + D(3).sqrt() / 2)) < D("1e-50")
    print("self-tests: pi, erf(1), sin and cos at pi/3, union area at distance 1: ok")
    lo, head, t1, t2 = bound(a, cells)
    print(f"a = {a}, cells = {cells}")
    print(f"int_0^a chi        >= {head:.15f}")
    print(f"isolated points    >= {t1:.15f}")
    print(f"isolated pairs     >= {t2:.15f}")
    print(f"certified: beta_MST >= {lo:.12f}")
    ok = lo >= claim
    print(f"claim {claim}: {'HOLDS' if ok else 'NOT IMPLIED'}")
    if controls:
        bad = 0
        c1 = claim + D("0.000001")
        r = not (lo >= c1); bad += not r
        print(f"red control 1 (claim {c1}): {'refused' if r else 'ACCEPTED (bad)'}")
        lo2, *_ = bound(a, 30)
        r = not (lo2 >= claim); bad += not r
        print(f"red control 2 (30 cells gives {lo2:.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        r = not (head + t1 >= claim); bad += not r
        print(f"red control 3 (no pair term, {head + t1:.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        lo4, *_ = bound("6/5", 2000)
        r = not (lo4 >= claim); bad += not r
        print(f"red control 4 (a = 6/5 gives {lo4:.9f}): {'refused' if r else 'ACCEPTED (bad)'}")
        if bad:
            print("RED CONTROLS FAILED"); return 2
        print("red controls: 4 of 4 refused")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
