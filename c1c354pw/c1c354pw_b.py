#!/usr/bin/env python3
"""C_1c, PR #354 (update of 13 November 2026): the piecewise-polynomial certificate, program B.

A second, separately written exact check, by a different method from program A at every step that matters.

  g, by jumps.  f is a compactly supported piecewise polynomial of degree at most n = 20, so
      f(x) = sum over breakpoints b, orders r of  D[b][r] (x - b)_+^r,
  where D[b][r] is the jump at b of the r-th Taylor coefficient (right minus left). Differentiating n + 1 times
  turns f into a finite sum of derivatives of point masses; doing the same for x -> f(-x), convolving, and
  integrating 2n + 2 times from -infinity gives the closed form

      g(t) = int f(x) f(x + t) dx
           = sum over b, b', r, r' of  (-1)^(r+1) D[b][r] D[b'][r'] r! r'! / (r + r' + 1)!  (t - (b' - b))_+^(r+r'+1).

  No integral of a product of pieces is taken (program A builds g from such integrals). The closed form is tested
  inside the program: the sum over ALL lags must vanish identically to the right of the support, and g is compared
  at eleven rational lags with a direct integral of f(x) f(x + t) over the common refinement of the two partitions.

  Positivity, by centred Taylor forms.  p >= floor on an interval is certified by
      p(c) - floor - sum_{k >= 1} |p^(k)(c) / k!| rad^k >= 0      (c the midpoint, rad the half-length),
  with bisection where that fails (program A uses Bernstein coefficients). A piece of f that vanishes at its left
  end is first divided by the corresponding power of (x - lo).

Everything is exact (fractions.Fraction). Standard library only; run with python3 -I. No program of the
certificate's repository was read or run.

usage: c1c354pw_b.py CERT.json [--bump PIECE COEFF NUM/DEN] [--claim NUM/DEN]
"""
import hashlib
import json
import re
import sys
import time
from fractions import Fraction
from math import factorial

NUMBER = re.compile(r"^-?[0-9]+(/[0-9]+)?$")


def number(text):
    if type(text) is not str or not NUMBER.match(text):
        raise SystemExit("B: not an exact rational: %r" % (text,))
    return Fraction(text)


def shift(coeffs, c):
    """Taylor coefficients of x -> p(x + c), by repeated synthetic division (Horner's scheme)."""
    work = coeffs[:]
    top = len(work) - 1
    for start in range(top):
        k = top - 1
        while k >= start:
            work[k] = work[k] + c * work[k + 1]
            k -= 1
    return work


def value(coeffs, x):
    total = Fraction(0)
    xp = Fraction(1)
    for c in coeffs:
        if c:
            total += c * xp
        xp *= x
    return total


def trim(coeffs):
    out = coeffs[:]
    while len(out) > 1 and out[-1] == 0:
        del out[-1]
    return out


def centred_ok(coeffs, length, floor, max_depth=200):
    """Certify sum coeffs[k] w^k >= floor on [0, length] by centred Taylor forms. Returns (ok, intervals, depth)."""
    todo = [(Fraction(0), Fraction(length), 0)]
    used = 0
    deepest = 0
    while todo:
        a, b, depth = todo.pop()
        mid = (a + b) / 2
        rad = (b - a) / 2
        local = shift(coeffs, mid)
        slack = local[0] - floor
        if slack < 0:
            return False, used, depth           # an actual value below the floor
        rp = Fraction(1)
        for c in local[1:]:
            rp *= rad
            slack -= abs(c) * rp
        if slack >= 0:
            used += 1
            deepest = max(deepest, depth)
            continue
        if depth >= max_depth:
            return False, used, depth
        todo.append((a, mid, depth + 1))
        todo.append((mid, b, depth + 1))
    return True, used, deepest


class Step:
    """f as a list of segments [x0, x1) with Taylor coefficients at x0; zero outside."""

    def __init__(self, segments):
        self.segments = segments
        self.points = [seg[0] for seg in segments] + [segments[-1][1]]

    def taylor_at(self, x):
        """Taylor coefficients of the piece containing x (taken on its right), about x."""
        for x0, x1, cs in self.segments:
            if x0 <= x < x1:
                return shift(cs, x - x0)
        return [Fraction(0)]

    def jumps(self):
        """D[b] = right Taylor coefficients at b minus left Taylor coefficients at b, for every breakpoint b."""
        table = []
        previous = None
        for x0, x1, cs in self.segments:
            left = shift(previous[2], previous[1] - previous[0]) if previous else []
            size = max(len(left), len(cs))
            table.append((x0, [(cs[r] if r < len(cs) else 0) - (left[r] if r < len(left) else 0)
                               for r in range(size)]))
            previous = (x0, x1, cs)
        last = shift(previous[2], previous[1] - previous[0])
        table.append((previous[1], [-c for c in last]))
        return table

    def correlate_at(self, t):
        """int f(x) f(x + t) dx by direct integration over the common refinement (used only as a spot check)."""
        cuts = sorted(set(self.points) | set(p - t for p in self.points))
        total = Fraction(0)
        for a, b in zip(cuts, cuts[1:]):
            p = self.taylor_at(a)
            q = self.taylor_at(a + t)
            if p == [0] or q == [0]:
                continue
            prod = [Fraction(0)] * (len(p) + len(q) - 1)
            for i, pi in enumerate(p):
                if pi:
                    for j, qj in enumerate(q):
                        prod[i + j] += pi * qj
            width = b - a
            wp = Fraction(1)
            for k, c in enumerate(prod):
                wp *= width
                total += c * wp / (k + 1)
        return total


def main(argv):
    began = time.time()
    with open(argv[1], "rb") as fh:
        doc = json.loads(fh.read().decode("utf-8"))
    s, t1, T, margin, M, bound = (number(doc[k]) for k in ("s", "t1", "T", "margin", "M", "bound"))
    h = [[number(lo), number(hi), [number(c) for c in cs]] for lo, hi, cs in doc["pieces"]]
    floor = T - margin
    if "--bump" in argv:
        at = argv.index("--bump")
        h[int(argv[at + 1])][2][int(argv[at + 2])] += Fraction(argv[at + 3])
        print("B: CONTROL: piece %s, coefficient %s changed by %s" % (argv[at + 1], argv[at + 2], argv[at + 3]))
    if "--claim" in argv:
        floor = Fraction(argv[argv.index("--claim") + 1])
        print("B: CONTROL: the lower bound on g replaced by %s" % float(floor))
    verdict = [True]

    def report(text, good):
        print("B: %-76s %s" % (text, "PASS" if good else "FAIL"))
        if not good:
            verdict[0] = False

    # the function
    half = s / 2
    report("s = 2^-k, and the spike ends at or before t1", s > 0 and (1 / s).denominator == 1
           and bin(int(1 / s)).count("1") == 1 and half <= t1)
    chain = h[0][0] == t1 and h[-1][1] == 1 + half
    for k in range(len(h)):
        chain = chain and h[k][0] < h[k][1] and (k == 0 or h[k - 1][1] == h[k][0])
    report("the pieces of h tile [t1, 1 + s/2)", chain)
    report("h = T on [1 - t1, 1) and on [1, 1 + s/2)", h[-2][:2] == [1 - t1, Fraction(1)] and h[-2][2] == [T]
           and h[-1][2] == [T])
    segments = [(-half, half, [1 / s]), (half, t1, [Fraction(0)])] + [tuple(p) for p in h]
    f = Step(segments)

    # f >= 0
    nonneg = True
    for x0, x1, cs in segments:
        body = cs[:]
        while len(body) > 1 and body[0] == 0:    # divide by (x - x0) while the piece vanishes there
            body = body[1:]
        if body == [0]:
            continue
        good, _, _ = centred_ok(body, x1 - x0, Fraction(0))
        nonneg = nonneg and good
    report("f >= 0 on every piece (centred forms, after removing zeros at the left end)", nonneg)

    # the integral of f
    mass = Fraction(0)
    for x0, x1, cs in segments:
        width = x1 - x0
        for k, c in enumerate(cs):
            mass += c * width ** (k + 1) / (k + 1)
    report("int f = 1 + M", mass == 1 + M)

    # g by the jump formula
    jumps = f.jumps()
    kernels = {}                                   # lag -> coefficients in powers of (t - lag)
    for b, Db in jumps:
        for b2, Db2 in jumps:
            lag = b2 - b
            acc = kernels.setdefault(lag, [Fraction(0)] * 44)
            for r, x in enumerate(Db):
                if not x:
                    continue
                sign = 1 if r % 2 else -1
                fr = factorial(r)
                for r2, y in enumerate(Db2):
                    if y:
                        acc[r + r2 + 1] += sign * x * y * Fraction(fr * factorial(r2), factorial(r + r2 + 1))
    lags = sorted(kernels)
    running = [Fraction(0)] * 44                   # sum of the kernels switched on so far, in powers of t
    cells = []
    for idx, lag in enumerate(lags):
        moved = shift(kernels[lag], -lag)
        running = [running[k] + moved[k] for k in range(44)]
        if lag >= 0 and idx + 1 < len(lags):
            cells.append((lag, lags[idx + 1], trim(shift(running, lag))))
    report("the closed form vanishes identically beyond the last lag", all(c == 0 for c in running))
    report("lags run from -(1 + s) to 1 + s and include 0 and 1", lags[0] == -(1 + s) and lags[-1] == 1 + s
           and 0 in kernels and 1 in kernels)
    print("B: %d lags, %d cells on [0, 1 + s], largest degree %d" % (
        len(lags), len(cells), max(len(c[2]) - 1 for c in cells)))

    # spot checks of the closed form against direct integration
    probes = [Fraction(0), s / 4, Fraction(3, 1000), Fraction(1, 7), t1, Fraction(1, 3), Fraction(1, 2),
              Fraction(7, 10), 1 - t1, Fraction(9, 10), Fraction(1)]
    agree = True
    for t in probes:
        for u, v, local in cells:
            if u <= t <= v:
                agree = agree and value(local, t - u) == f.correlate_at(t)
                break
        else:
            agree = False
    report("closed form = direct integral of f(x) f(x + t) at %d rational lags" % len(probes), agree)

    # positivity on [0, 1]
    fine = True
    used_all = 0
    deepest_all = 0
    for u, v, local in cells:
        if v > 1:
            continue
        good, used, deepest = centred_ok(local, v - u, floor)
        used_all += used
        deepest_all = max(deepest_all, deepest)
        if not good:
            fine = False
            print("B:   not certified, or below the bound, in the cell [%s, %s]" % (float(u), float(v)))
    report("g >= the stated bound on every cell of [0, 1] (%d intervals, depth <= %d)" % (used_all, deepest_all),
           fine)
    flat = [c for c in cells if c[1] <= 1 and c[2] == [T]]
    print("B:   cells with g = T identically: %d; g(1) - T = %s" % (
        len(flat), [value(l, v - u) for u, v, l in cells if v == 1][0] - T))

    # the bound
    report("bound = (T - margin) / (1 + M)^2", bound * (1 + M) ** 2 == T - margin)
    report("bound >= 0.410378195", bound * 10 ** 9 >= 410378195)
    scaled = bound.numerator * 10 ** 45 // bound.denominator
    print("B:   bound = 0.%045d..." % scaled)

    sha = hashlib.sha256()
    for u, v, local in cells:
        sha.update(("%s %s %s\n" % (u, v, " ".join(str(c) for c in local))).encode())
    print("B: digest of g (cells, coefficients in powers of t - u): %s" % sha.hexdigest())
    print("B: run time %.0f s" % (time.time() - began))
    print("B: VERDICT %s" % ("PASS" if verdict[0] else "REFUSED"))


if __name__ == "__main__":
    main(sys.argv)
