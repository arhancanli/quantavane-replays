#!/usr/bin/env python3
"""Program B: the covering certificate of the 13a upper bound [K2026] (PR #314), written separately from
c13acov_a.py (no shared code), from the same sources: the paper (Sections 4.1-4.4) and the definitions the author
stated on PR #314 on 14 November 2026. Certificate files read as data; no program of the source repository read or
run.

Differences from program A: directions are ordered by the half-angle tangent tau = y / (1 + x) (an exact, monotone
key of the angle), not by half-planes; the tree's splits are checked by marking a grid of (p, q) points; each row is
stored once as integer vectors over its own denominator, and each leaf is summed in integers over the least common
multiple of the denominators of the rows it uses; every family is evaluated from cumulative sums (differences of
prefix sums), not by running sums; B >= 1 is an integer comparison. It also checks every file's SHA-256 against the
manifest of the certificate directory.

Usage: python3 -I c13acov_b.py CERTS_ROOT [--control NAME]
"""
import gzip
import hashlib
import json
import os
import sys
import time
from fractions import Fraction
from math import gcd

NLAB = 63


def lcm(a, b):
    return a // gcd(a, b) * b


def unit_from_slope(t):
    """n(t) = ((1 - t^2) / (1 + t^2), 2t / (1 + t^2)) for a rational t."""
    t = Fraction(t)
    d = 1 + t * t
    return (Fraction(1) - t * t) / d, 2 * t / d


def tau(v):
    """Half-angle tangent of a unit vector (x, y): angle(v) in (-pi, pi) <-> tau in (-inf, inf); None at (-1, 0)."""
    x, y = v
    return None if x == -1 else y / (1 + x)


def angle_key(v):
    """Sort key for the angle in [0, 2 pi): angles in [0, pi) first (tau >= 0, increasing), then pi (tau = None),
    then (pi, 2 pi), i.e. tau < 0, increasing."""
    t = tau(v)
    if t is None:
        return (1, Fraction(0))
    return (0, t) if t >= 0 else (2, t)


def apply(M, v):
    (a, b), (c, d) = M
    return (a * v[0] + b * v[1], c * v[0] + d * v[1])


def main():
    root = sys.argv[1]
    control = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
    base = os.path.join(root, "k2026up")
    t0 = time.time()
    verdict = [True]

    def report(desc, good, detail=""):
        verdict[0] = verdict[0] and bool(good)
        print("B: check %s: %s%s" % (desc, "PASS" if good else "FAIL", (" (%s)" % detail) if detail else ""),
              flush=True)
        if not good and control:  # a red control stops at its first refusal
            print("B: VERDICT REFUSED")
            sys.exit(0)

    if control:
        print("B: CONTROL", control)
    manifest = json.load(open(os.path.join(base, "manifest.json"), encoding="utf-8"))
    bad = []
    for name, digest in manifest["sha256"].items():
        with open(os.path.join(base, name), "rb") as fh:
            if hashlib.sha256(fh.read()).hexdigest() != digest:
                bad.append(name)
    report("the SHA-256 of all 52 data files equals the certificate's manifest.json", not bad
           and len(manifest["sha256"]) == 52, "mismatch: %s" % bad if bad else "")
    header = json.load(open(os.path.join(base, "header.json"), encoding="utf-8"))
    sig = [Fraction(s) for s in header["t"]]
    hh = [Fraction(s) for s in header["h"]]
    kR, kL = header["kR"], header["kL"]
    specs = header["specs"]
    nrm = [unit_from_slope(s) for s in sig]
    report("normals n(sigma_i) are unit vectors and the side data are the paper's (4.1)",
           all(x * x + y * y == 1 for x, y in nrm) and header["t"] == ["-631/58", "-1", "7/16", "37/56"]
           and header["h"] == ["0", "0", "130783/250000", "446199/1000000"])
    # the area of K as (1/2) sum_i h_i l_i, l_i the length of side i (side i runs from the corner with side i - 1 to
    # the corner with side i + 1; its direction (-y_i, x_i) is a rational unit vector, so l_i is rational)
    def corner(i, j):
        A_ = [[nrm[i][0], nrm[i][1]], [nrm[j][0], nrm[j][1]]]
        dt = A_[0][0] * A_[1][1] - A_[0][1] * A_[1][0]
        return ((hh[i] * A_[1][1] - A_[0][1] * hh[j]) / dt, (A_[0][0] * hh[j] - hh[i] * A_[1][0]) / dt)
    lens = []
    for i in range(4):
        P, R = corner((i - 1) % 4, i), corner(i, (i + 1) % 4)
        tx, ty = -nrm[i][1], nrm[i][0]
        lens.append((R[0] - P[0]) * tx + (R[1] - P[1]) * ty)
    area = sum(hh[i] * lens[i] for i in range(4)) / 2
    exact = Fraction(7170601298323360123534337, 29109471743520000000000000)
    report("area of K by (1/2) sum h_i l_i (every side length positive) = the paper's fraction; 0.24634 and "
           "0.2463322373 are upper bounds", all(x > 0 for x in lens) and area == exact
           and area <= Fraction(24634, 100000) and area < Fraction(2463322373, 10 ** 10), "%.16f..." % float(area))
    # coefficients a_j0, a_j1 of a_j0 n_0 + a_j1 n_1 = -n_j, by Cramer's rule
    (x0, y0), (x1, y1) = nrm[0], nrm[1]
    dd = x0 * y1 - x1 * y0
    coefA = {}
    for j in (2, 3):
        xj, yj = nrm[j]
        coefA[j] = ((-xj * y1 + x1 * yj) / dd, (-x0 * yj + xj * y0) / dd)
    stated = {2: (Fraction(5541045, 8027539), Fraction(24427652, 40137695)),
              3: (Fraction(47299645, 118570699), Fraction(502233812, 592853495))}
    report("a_j0, a_j1 positive, solve the normal relations, and equal the author's four fractions",
           coefA == stated and all(c > 0 for j in coefA for c in coefA[j])
           and all(coefA[j][0] * nrm[0][k] + coefA[j][1] * nrm[1][k] + nrm[j][k] == 0 for j in (2, 3) for k in (0, 1)))

    def matrix(sp):
        if sp[0] == "flush":
            x, y = nrm[sp[1]]
            m = ((-y, x), (-x, -y)) if sp[2] == 0 else ((y, -x), (x, y))
            flip = sp[3]
        else:
            c, s = unit_from_slope(sp[1])
            m = ((c, -s), (s, c))
            flip = sp[2]
        return ((-m[0][0], -m[0][1]), m[1]) if flip == 1 else m

    mats = [matrix(sp) for sp in specs]
    report("each matrix has orthonormal columns", all(
        (m[0][0] ** 2 + m[1][0] ** 2, m[0][1] ** 2 + m[1][1] ** 2, m[0][0] * m[0][1] + m[1][0] * m[1][1]) == (1, 1, 0)
        for m in mats))
    found = []
    for sp, m in zip(specs, mats):
        if sp[0] in ("flush", "rot"):
            for v in nrm:
                w = apply(m, v)
                if w not in found:
                    found.append(w)
    ups = [w for w in found if w[0] == 0]
    rights = sorted([w for w in found if w[0] > 0], key=lambda w: w[1])
    lefts = sorted([w for w in found if w[0] < 0], key=lambda w: w[1], reverse=True)
    report("contact directions: (0, +-1) plus %d right and %d left, as kR = kL = 30" % (len(rights), len(lefts)),
           sorted(ups) == [(0, -1), (0, 1)] and len(rights) == kR == 30 and len(lefts) == kL == 30
           and len(set(w[1] for w in rights)) == 30 and len(set(w[1] for w in lefts)) == 30)
    dirs = [(Fraction(0), Fraction(-1)), (Fraction(0), Fraction(-1)), (Fraction(0), Fraction(1))] + rights + lefts
    if control == "dir":
        dirs[3 + 5] = (dirs[3 + 5][0], -dirs[3 + 5][1])
    which = {}
    for lab in range(NLAB - 1, 0, -1):  # (0, -1) ends at label 1 (P+), never 0
        which[dirs[lab]] = lab
    if control == "pminus":  # against the stated rule: (0, -1) from a spec given to P- instead of P+
        which[dirs[0]] = 0
    circle = sorted({dirs[lab] for lab in range(1, NLAB)}, key=angle_key)
    keys = [angle_key(w) for w in circle]
    expansions = []
    nonneg = True
    fr_ok = True
    for sp, m in zip(specs, mats):
        per = []
        for v in nrm:
            w = apply(m, v)
            if w in which:
                per.append(((which[w], Fraction(1)),))
                continue
            if sp[0] != "vrot":
                fr_ok = False
                per.append(())
                continue
            k = angle_key(w)
            nxt = 0
            while nxt < len(keys) and not keys[nxt] > k:
                nxt += 1
            if nxt == len(keys):
                nxt = 0
            w2, w1 = circle[nxt], circle[nxt - 1]
            det = w1[0] * w2[1] - w2[0] * w1[1]
            k1 = (w[0] * w2[1] - w2[0] * w[1]) / det
            k2 = (w1[0] * w[1] - w[0] * w1[1]) / det
            nonneg = nonneg and k1 >= 0 and k2 >= 0 and (k1 * w1[0] + k2 * w2[0], k1 * w1[1] + k2 * w2[1]) == w
            per.append(((which[w1], k1), (which[w2], k2)))
        expansions.append(per)
    report("flush and rot directions M n_i are all label directions", fr_ok)
    report("vrot directions: exact nonnegative two-term expansions in angularly adjacent label directions", nonneg)

    # rows as integer vectors over a row denominator: rowvec[k] = (den, {label: (X, Y)}, hnum, hden)
    def as_int(coef):
        den = 1
        for x, y in coef.values():
            den = lcm(den, lcm(x.denominator, y.denominator))
        return den, {lab: (int(x * den), int(y * den)) for lab, (x, y) in coef.items() if x or y}

    support_cache = {}

    def row(k):
        if k in support_cache:
            return support_cache[k]
        if k < 3969:
            lam, eta = k // NLAB, k % NLAB
            if lam == eta:
                return None
            ux, uy = dirs[lam]
            r = (as_int({lam: (ux, uy), eta: (-ux, -uy)}) if lam != eta else None, None)
        elif k == 3969:
            r = (as_int({1: (Fraction(1), Fraction(0)), 0: (Fraction(-1), Fraction(0))}), None)
        else:
            s, a = divmod(k - 3970, 2)
            if s >= len(specs):
                return None
            j = 2 + a
            coef = {}
            for i, wt in ((0, coefA[j][0]), (1, coefA[j][1]), (j, Fraction(1))):
                for lab, kap in expansions[s][i]:
                    cx, cy = coef.get(lab, (Fraction(0), Fraction(0)))
                    coef[lab] = (cx + wt * kap * dirs[lab][0], cy + wt * kap * dirs[lab][1])
            r = (as_int(coef), Fraction(446198, 1000000) if control == "h3" and j == 3 else hh[j])
        support_cache[k] = r
        return r

    # the tree, recursively, with a grid of (p, q) points for every split
    tree = json.loads(gzip.open(os.path.join(base, "tree.json.gz")).read())
    leaves = {}
    count = {"nodes": 0, "depth": 0, "ok": True}

    def visit(node, box, path, depth):
        count["nodes"] += 1
        count["depth"] = max(count["depth"], depth)
        kind = node[0]
        if kind == "l":
            if node[1] in leaves:
                count["ok"] = False
            leaves[node[1]] = (box, path)
        elif kind == "p":
            s = node[1]
            if any(x[0] == s for x in path):
                count["ok"] = False
            visit(node[2], box, path + ((s, 0),), depth + 1)
            visit(node[3], box, path + ((s, 1),), depth + 1)
        elif kind == "s":
            plo, phi, qlo, qhi = box
            grid = {}
            for cell, child in node[1]:
                a0, a1, b0, b1 = cell
                if control == "gap" and box == (0, kR, 0, kL) and cell == node[1][0][0]:
                    a1 -= 1
                for p in range(a0, a1 + 1):
                    for q in range(b0, b1 + 1):
                        grid[(p, q)] = grid.get((p, q), 0) + 1
                visit(child, (a0, a1, b0, b1), path, depth + 1)
            inside = all(plo <= p <= phi and qlo <= q <= qhi for p, q in grid)
            once = all(c == 1 for c in grid.values())
            full = len(grid) == (phi - plo + 1) * (qhi - qlo + 1)
            if not (inside and once and full):
                count["ok"] = False
        else:
            count["ok"] = False

    visit(tree, (0, kR, 0, kL), (), 0)
    report("tree: splits partition their rectangles (grid marking), specs not repeated on a path, leaves distinct",
           count["ok"])
    report("tree: %d nodes, %d leaves, depth %d; the paper's table says 24,898, 12,663 and 26"
           % (count["nodes"], len(leaves), count["depth"]),
           count["nodes"] == 24898 and sorted(leaves) == list(range(12663)) and count["depth"] == 26)

    # the leaves
    LR = [2 + i for i in range(1, kR + 1)]
    LL = [2 + kR + j for j in range(1, kL + 1)]
    worst_num, worst_den = 0, 1
    least_B = None
    failures = []
    total_rows = 0
    lines = 0
    for name in manifest["files"]:
        for raw in gzip.open(os.path.join(base, name)):
            lines += 1
            leaf, D, entries = json.loads(raw)
            if control == "mult" and leaf == 0:
                # every multiplier of leaf 0 doubled: B = 2 and the sum is still 0, but family norms grow fourfold
                entries = [[k, 2 * c] for k, c in entries]
            box, path = leaves.pop(leaf)
            allowed = set(path)
            good = D > 0
            G = 1
            parts = []
            for k, c in entries:
                total_rows += 1
                good = good and c > 0
                r = row(k) if k >= 0 else None
                if r is None or r[0] is None:
                    good = False
                    continue
                if k >= 3970 and divmod(k - 3970, 2) not in allowed:
                    good = False
                    continue
                parts.append((c, r))
                G = lcm(G, r[0][0])
            vx = [0] * NLAB
            vy = [0] * NLAB
            bnum = Fraction(0)
            for c, ((den, vec), hj) in parts:
                m = c * (G // den)
                for lab, (X, Y) in vec.items():
                    vx[lab] += m * X
                    vy[lab] += m * Y
                if hj is not None:
                    bnum += c * hj
            Bval = bnum / D
            least_B = Bval if least_B is None else min(least_B, Bval)
            good = good and bnum >= D and sum(vx) == 0 and sum(vy) == 0
            # true f = (vx, vy) / (D G); the bound |sum|^2 <= 1 becomes <= (D G)^2
            cap = (D * G) ** 2
            cx = [0] * (kL + 1)
            cy = [0] * (kL + 1)
            for j in range(1, kL + 1):
                cx[j] = cx[j - 1] + vx[LL[j - 1]]
                cy[j] = cy[j - 1] + vy[LL[j - 1]]
            rx = [0] * (kR + 1)
            ry = [0] * (kR + 1)
            for i in range(1, kR + 1):
                rx[i] = rx[i - 1] + vx[LR[i - 1]]
                ry[i] = ry[i - 1] + vy[LR[i - 1]]
            plo, phi, qlo, qhi = box
            sums = []
            for a in range(qlo + 1, kL + 1):          # {L_a, ..., L_b}
                for b in range(a, kL + 1):
                    sums.append((cx[b] - cx[a - 1], cy[b] - cy[a - 1]))
            for j in range(1, qhi + 2):                # {P-, L_j, ..., L_kL}
                sums.append((vx[0] + cx[kL] - cx[j - 1], vy[0] + cy[kL] - cy[j - 1]))
            for s in range(plo, kR + 1):               # {P+, R_1, ..., R_s}
                sums.append((vx[1] + rx[s], vy[1] + ry[s]))
            for a in range(0, phi + 1):                # {R_{a+1}, ..., R_b}
                for b in range(a, phi + 1):
                    sums.append((rx[b] - rx[a], ry[b] - ry[a]))
            top = max(x * x + y * y for x, y in sums)
            if top * worst_den > worst_num * cap:
                worst_num, worst_den = top, cap
            good = good and top <= cap
            if not good:
                failures.append(leaf)
                if control:
                    report("leaf %d: (4.7) and the row rules" % leaf, False)
    report("one multiplier line per leaf (%d lines, %d rows; the manifest says 6,778,020 rows)" % (lines, total_rows),
           not leaves and lines == 12663 and total_rows == manifest["normalized_rows"])
    worst = Fraction(worst_num, worst_den)
    report("all 12,663 leaves pass (4.7): B >= 1, sum f = 0, every family sum of squared norm <= 1, multipliers > 0, "
           "escape rows from the path only", not failures,
           "largest squared norm %.12f, least B %.12f" % (float(worst), float(least_B)) if not failures
           else "failing leaves %s" % failures[:10])
    print("B: largest |sum_I f|^2 over all leaves and families: %s" % worst)
    print("B: least B over all leaves: %s" % least_B)
    print("B: run time %.1f s" % (time.time() - t0))
    print("B: VERDICT %s" % ("PASS" if verdict[0] else "REFUSED"))


if __name__ == "__main__":
    main()
