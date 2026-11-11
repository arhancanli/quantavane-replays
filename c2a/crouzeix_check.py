#!/usr/bin/env python3
"""Exact replay of the 2a lower-bound row [C2007] (teorth/optimizationproblems, main 87b4ebb): C_2 >= 2, "achieved by
p(z) = z and A = (0 2; 0 0), for which W(A) is the unit disk".  The matrix is read from the page (pinned).

For p(z) = z the ratio is ||A|| / max_{z in W(A)} |z|.  Checked exactly (integers and rationals):
  1. ||A||^2 is the largest eigenvalue of A*A: its characteristic polynomial and roots (A*A is diagonal here);
  2. W(A) is the closed unit disk.  W(A) is compact and convex (Toeplitz-Hausdorff), so it is the intersection of the
     half-planes Re(e^{-i phi} z) <= h(phi), with h(phi) the largest eigenvalue of the Hermitian part
     H_phi = (e^{-i phi} A + e^{i phi} A*)/2.  For this A, H_phi = (0 e^{-i phi}; e^{i phi} 0): trace 0, determinant
     -|e^{-i phi}|^2 = -1, so h(phi) = 1 for every phi (checked symbolically: the entries are e^{-i phi} times A_12/2
     and its conjugate, with |A_12/2| = 1).  Hence W(A) is the closed unit disk and max |z| = 1.  Attainment is also
     checked directly: v = (1, 1)/sqrt 2 gives v*Av = 1; and |v*Av| <= 1 at 2,000 exact rational unit vectors
     (v = ((1 - s^2)/(1 + s^2), 2s/(1 + s^2)) with a rational phase on the second coordinate);
  3. the ratio ||A|| / 1 = 2.
Red controls (--controls): A = (0 2; 0 1) gives ||A|| = sqrt 5 and W(A) an ellipse with foci 0, 1 and minor axis 2, so
the ratio is sqrt 5 / ((1 + sqrt 5)/2) = 1.382 < 2 (the claim 2 fails); A = (0 3; 0 0) has W(A) of radius 3/2, not 1
(the unit-disk claim fails); a claimed ratio 2 + 10^-12 for the page's A fails.
usage: python3 -I crouzeix_check.py CERTS_ROOT [--controls]
"""
import math
import os
import re
import sys
from fractions import Fraction

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def read_matrix(path):
    text = open(path, encoding="utf-8").read()
    row = [l for l in text.splitlines() if "[C2007]" in l and "p(z)=z" in l]
    need(len(row) == 1, "the [C2007] row was not found")
    m = re.search(r"A=\\begin\{pmatrix\}\s*(-?\d+)\s*&\s*(-?\d+)\s*\\\\\\\\\s*(-?\d+)\s*&\s*(-?\d+)\s*\\end\{pmatrix\}", row[0])
    need(m is not None, "the matrix was not found in the row")
    a, b, c, d = (int(x) for x in m.groups())
    return [[a, b], [c, d]]


def norm_sq(A):
    """Largest eigenvalue of A^T A for a real 2x2 A, as (trace, det) and the exact root when the discriminant is a
    perfect square; returns a Fraction or None."""
    (a, b), (c, d) = A
    g11, g12, g22 = a * a + c * c, a * b + c * d, b * b + d * d
    tr, det = g11 + g22, g11 * g22 - g12 * g12
    disc = tr * tr - 4 * det
    r = math.isqrt(disc)
    if r * r != disc:
        return None, (tr, det)
    return Fraction(tr + r, 2), (tr, det)


def support_is_one(A):
    """For A = (0 b; 0 0): H_phi = (0 e^{-i phi} b/2; e^{i phi} conj(b)/2 0), eigenvalues +-|b|/2 for every phi."""
    (a, b), (c, d) = A
    need(a == c == d == 0, "the symbolic support function is written for (0 b; 0 0) only")
    return Fraction(abs(b), 2)


def main():
    ok = True
    A = read_matrix(os.path.join(ROOT, "page", "2a.md"))
    n2, (tr, det) = norm_sq(A)
    radius = support_is_one(A)
    # attainment and a sample of exact unit vectors (real and with a rational phase e^{i t}: cos t, sin t rational)
    sample_ok = True
    worst = Fraction(0)
    for s in range(-40, 41):
        for u in range(-12, 13):
            x1 = Fraction(1 - s * s, 1 + s * s) if s else Fraction(1)
            x2 = Fraction(2 * s, 1 + s * s)
            ct, st = Fraction(1 - u * u, 1 + u * u), Fraction(2 * u, 1 + u * u)       # phase of the second coordinate
            # v = (x1, x2 e^{it}): v* A v = conj(v1) * 2 * v2 = 2 x1 x2 e^{it}; |.|^2 = 4 x1^2 x2^2
            val2 = 4 * x1 * x1 * x2 * x2 * (ct * ct + st * st)
            worst = max(worst, val2)
            sample_ok &= val2 <= 1
    attain = Fraction(1, 2) * A[0][1]                       # v = (1, 1)/sqrt 2: v*Av = A_12/2
    ratio_sq = n2 / (radius * radius)
    good = (A == [[0, 2], [0, 0]] and n2 == 4 and radius == 1 and attain == 1 and sample_ok and ratio_sq == 4)
    print(f"A from the page: {A}; A^T A has trace {tr} and determinant {det}, largest eigenvalue {n2}, so ||A|| = 2")
    print(f"  Hermitian part H_phi = (0, {Fraction(A[0][1], 2)} e^(-i phi); {Fraction(A[0][1], 2)} e^(i phi), 0): "
          f"eigenvalues +-{radius} for every phi, so W(A) is the closed disk of radius {radius}")
    print(f"  v = (1, 1)/sqrt 2 gives v*Av = {attain}; over {81 * 25} exact unit vectors max |v*Av|^2 = {worst} <= 1: "
          f"{sample_ok} (the maximum 1 needs |v_1| = |v_2|, irrational here)")
    print(f"  ratio ||p(A)|| / max_W |p| = 2/{radius} = {2 / radius}")
    print(f"RESULT [C2007]: p(z) = z and A = (0 2; 0 0) give the ratio 2, so C_2 >= 2: {'HOLDS' if good else 'FAILS'}")
    ok &= good
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        B = [[0, 2], [0, 1]]
        nb, (trb, detb) = norm_sq(B)
        # W(B): ellipse with foci 0 and 1, minor axis sqrt(tr(B*B) - |0|^2 - |1|^2) = 2; farthest point from 0 is at
        # distance (|f1 - f2| + major axis)/2 with major axis sqrt(1 + 4) = sqrt 5 (elliptical range theorem)
        far = (1 + math.sqrt(5)) / 2
        rb = math.sqrt((trb + math.sqrt(trb * trb - 4 * detb)) / 2) / far
        print(f"  A = (0 2; 0 1): ||A|| = {math.sqrt((trb + math.sqrt(trb * trb - 4 * detb)) / 2):.6f}, "
              f"max_W |z| = {far:.6f}, ratio {rb:.6f}")
        res.append(("(0 2; 0 1) reaching 2", rb < 2))
        r3 = support_is_one([[0, 3], [0, 0]])
        res.append(("(0 3; 0 0) with W the unit disk", r3 != 1))
        print(f"  A = (0 3; 0 0): W(A) has radius {r3}")
        res.append(("a claimed ratio 2 + 10^-12", Fraction(2) < 2 + Fraction(1, 10 ** 12) and ratio_sq != (2 + Fraction(1, 10 ** 12)) ** 2))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
