#!/usr/bin/env python3
"""C_1c, PR #354 (update of 13 November 2026): the piecewise-polynomial certificate, program A.

The certificate (a JSON file of exact rationals) gives

    f = (1/s) 1_[-s/2, s/2) + h,

with h = 0 outside the listed pieces and h(x) = sum_k c_k (x - lo)^k on each piece [lo, hi). It states

    min over 0 <= t <= 1 of g(t) >= T - margin,   g(t) = int f(x) f(x + t) dx,   int f = 1 + M,

and hence C_1c >= (T - margin) / (1 + M)^2 = bound.

This program, written from the definition on the 1c page and the format description alone (no program of the
certificate's repository was read or run), checks all of it in exact rational arithmetic (fractions.Fraction):

  * the structure: contiguous pieces, the spike disjoint from them;
  * f >= 0: every polynomial piece is non-negative on its interval (Bernstein coefficients, with subdivision);
  * int f = 1 + M;
  * g as an exact polynomial on every cell of the lag mesh, built SYMBOLICALLY: for each ordered pair of pieces
    (i, j) the bivariate antiderivative A(y, tau) = int_0^y P_i(u) P_j(u + tau) du is expanded in y and tau, and the
    limits of integration (0 or -tau below, L_i or L_j - tau above) are substituted on each range of tau;
  * consistency of that g: continuity at every mesh point, g(0) = int f^2, 2 int_0^inf g = (int f)^2, g = 0 at the
    end of its support;
  * g >= T - margin on every cell inside [0, 1] (Bernstein coefficients of g - (T - margin), with subdivision);
  * bound = (T - margin) / (1 + M)^2 exactly, and its decimal expansion.

It prints a SHA-256 digest of the exact piecewise form of g (cell by cell, coefficients in powers of t - u), which
program B, written separately and by a different method, must reproduce.

Standard library only; run with python3 -I.
usage: c1c354pw_a.py CERT.json [--bump PIECE COEFF NUM/DEN] [--claim NUM/DEN]
  --bump adds a rational to one coefficient of one piece of h (red control);
  --claim replaces the stated lower bound T - margin on g by another value (red control).
"""
import hashlib
import json
import sys
import time
from fractions import Fraction as F
from math import comb

ZERO = F(0)


def load(path):
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)

    def rat(x):
        if not isinstance(x, str):
            raise SystemExit("A: a number in the certificate is not a string")
        return F(x)

    cert = {k: rat(raw[k]) for k in ("t1", "s", "T", "margin", "M", "bound")}
    cert["pieces"] = [(rat(lo), rat(hi), [rat(c) for c in cs]) for lo, hi, cs in raw["pieces"]]
    return cert


def taylor_shift(p, c):
    """coefficients of p(x + c)"""
    q = list(p)
    n = len(q)
    for i in range(n - 1):
        for j in range(n - 2, i - 1, -1):
            q[j] += c * q[j + 1]
    return q


def peval(p, x):
    acc = ZERO
    for c in reversed(p):
        acc = acc * x + c
    return acc


def pint(p, a, b):
    """exact integral of p over [a, b]"""
    anti = [ZERO] + [c / (k + 1) for k, c in enumerate(p)]
    return peval(anti, b) - peval(anti, a)


def bernstein_min(local, length, floor):
    """Decide whether sum_k local[k] w^k >= floor for all 0 <= w <= length.

    Returns (ok, certified lower bound of p - floor on the interval, leaves, depth). The Bernstein coefficients of a
    polynomial on an interval bound it from below, and the first and last coefficients are its end values.
    """
    n = len(local) - 1
    scaled = [c * length ** k for k, c in enumerate(local)]
    scaled[0] -= floor
    bern = []
    for i in range(n + 1):
        acc = ZERO
        for k in range(i + 1):
            if scaled[k]:
                acc += F(comb(i, k), comb(n, k)) * scaled[k]
        bern.append(acc)
    stack = [(bern, 0)]
    low = None
    leaves = 0
    deepest = 0
    while stack:
        b, depth = stack.pop()
        if b[0] < 0 or b[-1] < 0:
            return False, min(b[0], b[-1]), leaves, depth
        m = min(b)
        if m >= 0:
            leaves += 1
            deepest = max(deepest, depth)
            low = m if low is None else min(low, m)
            continue
        if depth >= 400:
            return False, m, leaves, depth
        left, right, cur = [b[0]], [b[-1]], b
        while len(cur) > 1:
            cur = [(cur[i] + cur[i + 1]) / 2 for i in range(len(cur) - 1)]
            left.append(cur[0])
            right.append(cur[-1])
        right.reverse()
        stack.append((left, depth + 1))
        stack.append((right, depth + 1))
    return True, low, leaves, deepest


def antiderivative(ci, cj):
    """A[k][n] = coefficient of y^k tau^n in int_0^y P_i(u) P_j(u + tau) du."""
    di, dj = len(ci) - 1, len(cj) - 1
    prod = [[ZERO] * (dj + 1) for _ in range(di + dj + 1)]
    for m in range(dj + 1):
        for n in range(dj + 1 - m):
            e = cj[m + n] * comb(m + n, m)       # P_j(u + tau) = sum e u^m tau^n
            if not e:
                continue
            for a in range(di + 1):
                if ci[a]:
                    prod[a + m][n] += ci[a] * e
    out = [[ZERO] * (dj + 1)]
    for k, row in enumerate(prod):
        out.append([c / (k + 1) for c in row])
    return out


def at_const(A, c):
    out = [ZERO] * len(A[0])
    ck = F(1)
    for row in A:
        for n, a in enumerate(row):
            if a:
                out[n] += a * ck
        ck *= c
    return out


def at_minus_tau(A):
    out = [ZERO] * (len(A) + len(A[0]) - 1)
    for k, row in enumerate(A):
        sign = -1 if k % 2 else 1
        for n, a in enumerate(row):
            if a:
                out[k + n] += sign * a
    return out


def at_c_minus_tau(A, c):
    out = [ZERO] * (len(A) + len(A[0]) - 1)
    power = [F(1)]                               # (c - tau)^k
    for k, row in enumerate(A):
        if k:
            nxt = [ZERO] * (len(power) + 1)
            for r, co in enumerate(power):
                nxt[r] += c * co
                nxt[r + 1] -= co
            power = nxt
        for n, a in enumerate(row):
            if a:
                for r, co in enumerate(power):
                    if co:
                        out[n + r] += a * co
    return out


def psub(p, q):
    n = max(len(p), len(q))
    return [(p[k] if k < len(p) else ZERO) - (q[k] if k < len(q) else ZERO) for k in range(n)]


def main(argv):
    t0 = time.time()
    cert = load(argv[1])
    floor_override = None
    if "--bump" in argv:
        j = argv.index("--bump")
        pi, ci, delta = int(argv[j + 1]), int(argv[j + 2]), F(argv[j + 3])
        cert["pieces"][pi][2][ci] += delta
        print("A: CONTROL: piece %d, coefficient %d changed by %s" % (pi, ci, delta))
    if "--claim" in argv:
        floor_override = F(argv[argv.index("--claim") + 1])
        print("A: CONTROL: the lower bound on g replaced by %s" % float(floor_override))
    s, t1, T, margin, M, bound = (cert[k] for k in ("s", "t1", "T", "margin", "M", "bound"))
    hp = cert["pieces"]
    ok = True

    def check(label, cond):
        nonlocal ok
        print("A: %-74s %s" % (label, "PASS" if cond else "FAIL"))
        ok = ok and bool(cond)

    # structure
    print("A: s = 2^-%d, t1 = %s, %d pieces, degrees %s" % (
        s.denominator.bit_length() - 1, t1, len(hp), [len(c) - 1 for _, _, c in hp]))
    check("s is a positive power of 1/2 and the spike [-s/2, s/2) lies left of t1",
          s.numerator == 1 and s.denominator & (s.denominator - 1) == 0 and 0 < s / 2 <= t1)
    check("pieces contiguous from t1 to 1 + s/2, each of positive length",
          hp[0][0] == t1 and hp[-1][1] == 1 + s / 2 and all(lo < hi for lo, hi, _ in hp)
          and all(hp[k][1] == hp[k + 1][0] for k in range(len(hp) - 1)))
    check("the last two pieces are the constant T on [1 - t1, 1) and [1, 1 + s/2)",
          hp[-2][2] == [T] and hp[-1][2] == [T] and hp[-2][0] == 1 - t1 and hp[-2][1] == 1)
    pieces = [(-s / 2, s / 2, [1 / s])] + hp

    # f >= 0
    worst = None
    nonneg = True
    for lo, hi, cs in pieces:
        good, low, _, _ = bernstein_min(cs, hi - lo, ZERO)
        nonneg = nonneg and good
        if good:
            worst = low if worst is None else min(worst, low)
    check("f >= 0 on every piece (Bernstein)", nonneg)
    if nonneg:
        print("A:   smallest Bernstein coefficient over the pieces: %s" % float(worst))

    # integrals of f and f^2
    integral = sum(pint(cs, ZERO, hi - lo) for lo, hi, cs in pieces)
    check("int f = 1 + M exactly", integral == 1 + M)
    print("A:   int f = %.30f" % float(integral))
    square = ZERO
    for lo, hi, cs in pieces:
        sq = [ZERO] * (2 * len(cs) - 1)
        for a, ca in enumerate(cs):
            for b, cb in enumerate(cs):
                sq[a + b] += ca * cb
        square += pint(sq, ZERO, hi - lo)

    # g, pair by pair
    ranges = []                                   # (t_from, t_to, polynomial in t)
    for lo_i, hi_i, ci in pieces:
        Li = hi_i - lo_i
        for lo_j, hi_j, cj in pieces:
            if hi_j <= lo_i:                      # no overlap for any t >= 0
                continue
            Lj = hi_j - lo_j
            d = lo_j - lo_i
            A = antiderivative(ci, cj)
            cuts = sorted({-Li, ZERO, Lj - Li, Lj})
            cache = {}
            for ta, tb in zip(cuts, cuts[1:]):
                if d + tb <= 0:
                    continue
                mid = (ta + tb) / 2
                upper = "L" if mid < Lj - Li else "R"
                lower = "N" if mid < 0 else "Z"
                if upper not in cache:
                    cache[upper] = at_const(A, Li) if upper == "L" else at_c_minus_tau(A, Lj)
                poly = cache[upper]
                if lower == "N":
                    if "N" not in cache:
                        cache["N"] = at_minus_tau(A)
                    poly = psub(poly, cache["N"])
                ranges.append((max(d + ta, ZERO), d + tb, taylor_shift(poly, -d)))
    mesh = sorted({r[0] for r in ranges} | {r[1] for r in ranges})
    cells = []
    for u, v in zip(mesh, mesh[1:]):
        total = []
        for a, b, p in ranges:
            if a <= u and v <= b:
                if len(p) > len(total):
                    total += [ZERO] * (len(p) - len(total))
                for k, c in enumerate(p):
                    total[k] += c
        local = taylor_shift(total, u)
        while len(local) > 1 and local[-1] == 0:
            local.pop()
        cells.append((u, v, local))
    print("A: lag mesh: %d points from %s to 1 + s; %d cells, largest degree %d" % (
        len(mesh), mesh[0], len(cells), max(len(c[2]) - 1 for c in cells)))
    check("the mesh starts at 0, contains 1 and ends at 1 + s", mesh[0] == 0 and 1 in mesh and mesh[-1] == 1 + s)

    # consistency of g
    check("g is continuous at every interior mesh point",
          all(peval(cells[k][2], cells[k][1] - cells[k][0]) == cells[k + 1][2][0] for k in range(len(cells) - 1)))
    check("g(0) = int f^2", cells[0][2][0] == square)
    check("g = 0 at the end of the support", peval(cells[-1][2], cells[-1][1] - cells[-1][0]) == 0)
    area = sum(pint(local, ZERO, v - u) for u, v, local in cells)
    check("2 int_0^inf g = (int f)^2", 2 * area == integral * integral)

    # g >= T - margin on [0, 1]
    floor = T - margin if floor_override is None else floor_override
    all_good = True
    low_all = None
    leaves_all = 0
    depth_all = 0
    smallest_end = None
    for u, v, local in cells:
        if v > 1:
            continue
        good, low, leaves, depth = bernstein_min(local, v - u, floor)
        if not good:
            all_good = False
            print("A:   g < the stated bound somewhere in the cell [%s, %s] (value - bound = %s)" % (
                float(u), float(v), float(low)))
            continue
        low_all = low if low_all is None else min(low_all, low)
        leaves_all += leaves
        depth_all = max(depth_all, depth)
        for val in (local[0], peval(local, v - u)):
            smallest_end = val if smallest_end is None else min(smallest_end, val)
    check("g >= the stated bound on every cell of [0, 1] (Bernstein, %d leaves, depth <= %d)" % (leaves_all, depth_all),
          all_good)
    if all_good:
        print("A:   certified: min over [0, 1] of g, minus the stated bound, >= %s" % low_all)
        print("A:   smallest g at a mesh point, minus T: %s" % (smallest_end - T))
        exact_T = [(float(u), float(v)) for u, v, local in cells if v <= 1 and local == [T]]
        print("A:   cells on which g = T identically: %d, covering [%s, %s]" % (
            len(exact_T), exact_T[0][0] if exact_T else "-", exact_T[-1][1] if exact_T else "-"))

    # the bound
    value = (T - margin) / (1 + M) ** 2
    check("bound = (T - margin) / (1 + M)^2 exactly", value == bound)
    digits = 45
    shown = bound.numerator * 10 ** digits // bound.denominator
    print("A:   bound = 0.%0*d... (denominator of %d bits)" % (digits, shown, bound.denominator.bit_length()))
    check("bound >= 0.410378195 (the printed row)", bound >= F(410378195, 10 ** 9))

    h = hashlib.sha256()
    for u, v, local in cells:
        h.update(("%s %s %s\n" % (u, v, " ".join(str(c) for c in local))).encode())
    print("A: digest of g (cells, coefficients in powers of t - u): %s" % h.hexdigest())
    print("A: run time %.0f s" % (time.time() - t0))
    print("A: VERDICT %s" % ("PASS" if ok else "REFUSED"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
