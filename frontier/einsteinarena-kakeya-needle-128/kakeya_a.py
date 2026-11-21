#!/usr/bin/env python3
"""Program A: exact area of the union of the 128 triangles of EinsteinArena's 'Discretized Kakeya Needle (n = 128)'.
T_j = conv{(x_j, 0), (x_j + 1/128, 0), (x_j + j/128, 1)}, j = 1..128. Method: horizontal slices. At height y the
section of T_j is [x_j + j y/128, x_j + 1/128 + (j-1) y/128]; the union length U(y) is continuous and affine between
consecutive heights where two of the 256 endpoint lines cross, so the area is the trapezoid sum over those heights.
All coordinates are scaled to integers (common denominator of the offsets, times 128); U is evaluated at each
breakpoint y = p/q as an exact integer over q; the sum is an exact Fraction. Standard library only.

Usage: python3 -I kakeya_a.py best_26.json SOLUTION_ID [--control shift|drop]
"""
import json
import sys
from fractions import Fraction
from math import gcd


def parse(s):
    s = s.strip()
    neg = s.startswith("-")
    s = s.lstrip("+-")
    ip, _, fp = s.partition(".")
    assert (ip + fp).isdigit() and len(fp) <= 25
    v = Fraction(int(ip + fp), 10 ** len(fp))
    return -v if neg else v


def area(xs):
    n = len(xs)
    D = 1
    for x in xs:
        D = D * x.denominator // gcd(D, x.denominator)
    X = [int(x * D) for x in xs]            # x_j = X_j / D
    # endpoint lines scaled by 128 D:  left_j(y) = 128 X_j + j D y,  right_j(y) = 128 X_j + D + (j-1) D y
    L = [(128 * X[j - 1], j * D) for j in range(1, n + 1)]
    R = [(128 * X[j - 1] + D, (j - 1) * D) for j in range(1, n + 1)]
    lines = L + R
    ys = {Fraction(0), Fraction(1)}
    for i in range(len(lines)):
        a0, a1 = lines[i]
        for k in range(i + 1, len(lines)):
            b0, b1 = lines[k]
            if a1 != b1:
                y = Fraction(b0 - a0, a1 - b1)
                if 0 < y < 1:
                    ys.add(y)
    ys = sorted(ys)

    def U(y):  # union length at height y, times 128 D, as a Fraction
        p, q = y.numerator, y.denominator
        iv = sorted((L[j][0] * q + L[j][1] * p, R[j][0] * q + R[j][1] * p) for j in range(n))
        tot, cl, cr = 0, iv[0][0], iv[0][1]
        for l, r in iv[1:]:
            if l > cr:
                tot += cr - cl
                cl, cr = l, r
            elif r > cr:
                cr = r
        tot += cr - cl
        return Fraction(tot, q)

    total = Fraction(0)
    prev_y, prev_u = ys[0], U(ys[0])
    for y in ys[1:]:
        u = U(y)
        total += (prev_u + u) * (y - prev_y) / 2
        prev_y, prev_u = y, u
    return total / (128 * D), len(ys)


def main():
    path, sid = sys.argv[1], sys.argv[2]
    control = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else ""
    sol = next(s for s in json.load(open(path)) if str(s["id"]) == sid)
    data = sol["data"] if isinstance(sol["data"], dict) else json.loads(sol["data"])
    xs = [parse(s) for s in data["offsets"]]
    assert len(xs) == 128
    if control == "shift":     # red control: move the 64th triangle by 1/1000 (the area must change)
        xs[63] += Fraction(1, 1000)
    if control == "translate":  # invariance control: translate all triangles (the area must not change)
        xs = [x + Fraction(7, 3) for x in xs]
    A, nb = area(xs)
    print("solution %s by %s, stated score %s, created %s" % (sol["id"], sol["agentName"], sol["score"], sol["createdAt"]))
    print("control: %s" % (control or "none"))
    print("breakpoint heights (with 0 and 1): %d" % nb)
    print("exact area = %d/%d" % (A.numerator, A.denominator) if A.denominator < 10 ** 400 else "exact area has a %d-digit denominator" % len(str(A.denominator)))
    s = str(A.numerator * 10 ** 30 // A.denominator)
    print("area = 0.%s... (30 digits, truncated)" % s.rjust(30, "0"))
    f = float(A)
    print("float(area) = %r; equals the stated score: %s" % (f, repr(f) == sol["score"] or f == float(sol["score"])))
    print("RESULT area_fraction_sha=%s" % __import__("hashlib").sha256(("%d/%d" % (A.numerator, A.denominator)).encode()).hexdigest()[:16])


if __name__ == "__main__":
    main()
