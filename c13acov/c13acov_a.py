#!/usr/bin/env python3
"""Program A: the covering certificate of the 13a upper bound [K2026] (Keller, arXiv 2609.21968v1, Section 4.4),
C_13a <= area(K) = 7170601298323360123534337/29109471743520000000000000 (PR #314).

Written from the paper (Sections 4.1-4.4, read as text) and from the definitions the author stated on PR #314 on
14 November 2026 (the contact directions, the matrix of each spec kind, the rows and the expansion rule). The
certificate files (ethan-keller/moser-worm-improved-bounds @ 9a33941, certificates/upper/) are read as data; no
program of that repository was read or run. Exact rational arithmetic (Python fractions); per leaf the vectors
f_lambda are brought to one common denominator and every prefix check is an integer comparison.

What is checked:
  1. the side data: n_i = n(sigma_i) unit vectors, a_j0 n_0 + a_j1 n_1 + n_j = 0 with a_j0, a_j1 > 0 (solved here;
     equal to the four fractions the author states);
  2. every spec matrix is orthogonal; the 25 flush/rot specs give, after removing duplicates, exactly the vertical
     directions (0, +-1) and 30 + 30 rational unit directions u^R (x > 0, sorted by increasing y) and u^L (x < 0,
     sorted by decreasing y), strictly ordered;
  3. every M_s n_i (1463 specs, i = 0..3) is a label direction (flush/rot) or, for vrot, a label direction or a
     nonnegative combination of the two angularly adjacent label directions, exactly;
  4. the tree: root [0,30] x [0,30]; every split's cells lie in the parent, are pairwise disjoint and cover it;
     no spec branched twice on a path; every leaf index 0..12662 exactly once;
  5. every leaf: D > 0 and every n_k > 0; rows valid (k < 3969 with lambda != eta, k = 3969, escape rows only for
     the (s, a) chosen on the leaf's path); B >= 1; sum of f_lambda = 0; |sum_{lambda in I} f_lambda|^2 <= 1 for
     every I in the four families F(R) of the paper's table (k_R = k_L = 30).

Usage: python3 -I c13acov_a.py CERTS_ROOT [--control NAME]
"""
import gzip
import json
import os
import sys
import time
from fractions import Fraction as Q
from math import lcm

STATED_A = {(2, 0): Q(5541045, 8027539), (2, 1): Q(24427652, 40137695),
            (3, 0): Q(47299645, 118570699), (3, 1): Q(502233812, 592853495)}


def normal(s):
    return ((1 - s * s) / (1 + s * s), 2 * s / (1 + s * s))


def mat_apply(M, v):
    return (M[0][0] * v[0] + M[0][1] * v[1], M[1][0] * v[0] + M[1][1] * v[1])


def spec_matrix(spec, n):
    kind = spec[0]
    if kind == "flush":
        i, a, b = spec[1], spec[2], spec[3]
        x, y = n[i]
        M = [[-y, x], [-x, -y]] if a == 0 else [[y, -x], [x, y]]
        mirror = b
    elif kind in ("rot", "vrot"):
        c, s = normal(Q(spec[1]))
        M = [[c, -s], [s, c]]
        mirror = spec[2]
    else:
        raise ValueError("unknown spec kind %r" % kind)
    assert mirror in (0, 1)
    if mirror:
        M = [[-M[0][0], -M[0][1]], M[1]]
    return M


def main():
    root = sys.argv[1]
    control = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
    d = os.path.join(root, "k2026up")
    t0 = time.time()
    ok = True

    def check(desc, good, detail=""):
        nonlocal ok
        ok = ok and bool(good)
        print("A: check %s: %s%s" % (desc, "PASS" if good else "FAIL", (" (%s)" % detail) if detail else ""),
              flush=True)
        if not good and control:  # a red control stops at its first refusal, as the paper's checker rejects
            print("A: VERDICT REFUSED")
            sys.exit(0)

    if control:
        print("A: CONTROL", control)
    with open(os.path.join(d, "header.json"), "rb") as fh:
        hdr = json.loads(fh.read().decode("utf-8"))
    sigma = [Q(x) for x in hdr["t"]]
    h = [Q(x) for x in hdr["h"]]
    kR, kL = hdr["kR"], hdr["kL"]
    specs = hdr["specs"]
    n = [normal(s) for s in sigma]
    check("side data as in the paper: sigma = (-631/58, -1, 7/16, 37/56), h = (0, 0, 130783/250000, "
          "446199/1000000), kR = kL = 30", sigma == [Q(-631, 58), Q(-1), Q(7, 16), Q(37, 56)]
          and h == [0, 0, Q(130783, 250000), Q(446199, 1000000)] and kR == kL == 30 and hdr["leaves"] == 12663)
    # 0. the quadrilateral K and its area (shoelace over the four vertices of consecutive side lines)
    def meet(i, j):
        (a, b), (c, e) = n[i], n[j]
        dt = a * e - b * c
        return ((h[i] * e - h[j] * b) / dt, (a * h[j] - c * h[i]) / dt)
    V = [meet(0, 1), meet(1, 2), meet(2, 3), meet(3, 0)]
    area = sum(V[i][0] * V[(i + 1) % 4][1] - V[(i + 1) % 4][0] * V[i][1] for i in range(4)) / 2
    inside = all(n[i][0] * v[0] + n[i][1] * v[1] <= h[i] for v in V for i in range(4))
    diag = (V[1][0] - V[3][0]) ** 2 + (V[1][1] - V[3][1]) ** 2
    AREA = Q(7170601298323360123534337, 29109471743520000000000000)
    check("K: the four vertices satisfy all four side inequalities; shoelace area = 7170601298323360123534337/"
          "29109471743520000000000000; |V2 - V4|^2 > 1", inside and area == AREA and diag > 1,
          "area %s = %.16f..., |V2 - V4|^2 = %.6f" % (area, float(area), float(diag)))
    check("the printed upper bounds 0.24634 (the README cell of PR #314) and 0.2463322373 hold, and 0.2463322372 "
          "is the truncation", area <= Q(24634, 100000) and Q(2463322372, 10 ** 10) <= area < Q(2463322373, 10 ** 10))
    # 1. normal relations
    A = {}
    det = n[0][0] * n[1][1] - n[0][1] * n[1][0]
    for j in (2, 3):
        # solve a0 n0 + a1 n1 = -n_j
        rx, ry = -n[j][0], -n[j][1]
        a0 = (rx * n[1][1] - ry * n[1][0]) / det
        a1 = (n[0][0] * ry - n[0][1] * rx) / det
        A[(j, 0)], A[(j, 1)] = a0, a1
    check("a_j0 n_0 + a_j1 n_1 + n_j = 0 has positive solutions equal to the four fractions stated on PR #314",
          A == STATED_A and all(v > 0 for v in A.values()), ", ".join("a%d%d = %s" % (j, i, A[(j, i)])
                                                                          for (j, i) in sorted(A)))
    # 2. spec matrices and the contact directions
    mats = [spec_matrix(sp, n) for sp in specs]
    kinds = [sp[0] for sp in specs]
    check("1463 specs: 16 flush, 9 rot, 1438 vrot", len(specs) == 1463 and kinds.count("flush") == 16
          and kinds.count("rot") == 9 and kinds.count("vrot") == 1438)
    orth = all(M[0][0] * M[0][0] + M[1][0] * M[1][0] == 1 and M[0][1] * M[0][1] + M[1][1] * M[1][1] == 1
               and M[0][0] * M[0][1] + M[1][0] * M[1][1] == 0 for M in mats)
    check("every spec matrix is orthogonal (exactly)", orth)
    dirs = set()
    for M, kind in zip(mats, kinds):
        if kind != "vrot":
            for i in range(4):
                dirs.add(mat_apply(M, n[i]))
    vert = sorted(v for v in dirs if v[0] == 0)
    right = sorted((v for v in dirs if v[0] > 0), key=lambda v: v[1])
    left = sorted((v for v in dirs if v[0] < 0), key=lambda v: -v[1])
    if control == "swapRL":
        right[0], right[1] = right[1], right[0]
    check("the 25 flush/rot specs give the vertical directions (0, -1), (0, 1) and 30 + 30 others",
          vert == [(0, -1), (0, 1)] and len(right) == kR and len(left) == kL, "%d distinct" % len(dirs))
    check("u^R strictly upward in the open right half-plane, u^L strictly downward in the open left half-plane, "
          "all unit vectors", all(right[i][1] < right[i + 1][1] for i in range(kR - 1))
          and all(left[i][1] > left[i + 1][1] for i in range(kL - 1))
          and all(v[0] * v[0] + v[1] * v[1] == 1 for v in right + left))
    U = [(Q(0), Q(-1)), (Q(0), Q(-1)), (Q(0), Q(1))] + right + left  # labels: P-, P+, P0, R1..R30, L1..L30
    NL = len(U)
    label_of = {(Q(0), Q(-1)): 1, (Q(0), Q(1)): 2}
    for i, v in enumerate(right):
        label_of[v] = 3 + i
    for j, v in enumerate(left):
        label_of[v] = 3 + kR + j

    # 3. expansions of M_s n_i into label directions
    def quadrant_key(v):
        # angle order on [0, 2 pi): (half, then cross products); returns a sortable key with exact comparisons
        x, y = v
        half = 0 if (y > 0 or (y == 0 and x > 0)) else 1
        # within a half-plane, angle increases with -x/|v| ... use the cotangent: for half 0, angle in [0, pi):
        # decreasing x / y; compare by atan2-free key (-x, y) on the unit circle: x = cos, monotone decreasing in
        # angle on [0, pi] and on [pi, 2 pi] increasing
        return (half, -x if half == 0 else x)

    ring = sorted(set([(Q(0), Q(-1)), (Q(0), Q(1))] + right + left), key=quadrant_key)
    assert len(ring) == 62
    expand_ok = True
    nonneg_ok = True
    vrot_two = 0
    EXP = []  # EXP[s][i] = list of (label, kappa)
    for s, (M, kind) in enumerate(zip(mats, kinds)):
        row = []
        for i in range(4):
            v = mat_apply(M, n[i])
            if v in label_of:
                row.append([(label_of[v], Q(1))])
                continue
            if kind != "vrot":
                expand_ok = False
                row.append([])
                continue
            kv = quadrant_key(v)
            pos = next((p for p, w in enumerate(ring) if quadrant_key(w) > kv), 0)
            w2 = ring[pos]
            w1 = ring[pos - 1]
            dt = w1[0] * w2[1] - w1[1] * w2[0]
            k1 = (v[0] * w2[1] - v[1] * w2[0]) / dt
            k2 = (w1[0] * v[1] - w1[1] * v[0]) / dt
            if not (k1 >= 0 and k2 >= 0 and k1 * w1[0] + k2 * w2[0] == v[0] and k1 * w1[1] + k2 * w2[1] == v[1]):
                nonneg_ok = False
            vrot_two += 1
            row.append([(label_of[w1], k1), (label_of[w2], k2)])
        EXP.append(row)
    check("every M_s n_i of a flush/rot spec is a label direction", expand_ok)
    check("every other M_s n_i of a vrot spec is kappa_1 w_1 + kappa_2 w_2 with adjacent label directions and "
          "kappa >= 0, exactly", nonneg_ok, "%d two-term expansions" % vrot_two)

    # escape row coefficient vectors: row 3970 + 2s + a -> {label: vector}, right side h_j
    def escape_row(s, a):
        j = 2 + a
        coef = {}
        for i, w in ((0, A[(j, 0)]), (1, A[(j, 1)]), (j, Q(1))):
            for lab, kap in EXP[s][i]:
                c = coef.get(lab, (Q(0), Q(0)))
                coef[lab] = (c[0] + w * kap * U[lab][0], c[1] + w * kap * U[lab][1])
        return coef, (Q(130782, 250000) if control == "h2" and j == 2 else h[j])

    ESC = {}

    # 4. the tree
    with gzip.open(os.path.join(d, "tree.json.gz"), "rb") as fh:
        tree = json.loads(fh.read().decode("utf-8"))
    leaf_info = {}
    stats = {"nodes": 0, "depth": 0, "split_ok": True, "branch_ok": True, "dup_leaf": False}

    def cells_ok(parent, cells):
        p0, p1, q0, q1 = parent
        area = 0
        for c in cells:
            a0, a1, b0, b1 = c
            if not (p0 <= a0 <= a1 <= p1 and q0 <= b0 <= b1 <= q1):
                return False
            area += (a1 - a0 + 1) * (b1 - b0 + 1)
        for x in range(len(cells)):
            for y in range(x + 1, len(cells)):
                a, b = cells[x], cells[y]
                if a[0] <= b[1] and b[0] <= a[1] and a[2] <= b[3] and b[2] <= a[3]:
                    return False
        return area == (p1 - p0 + 1) * (q1 - q0 + 1)

    stack = [(tree, (0, kR, 0, kL), frozenset(), frozenset(), 0)]  # depth: edges from the root
    while stack:
        node, rect, used, chosen, depth = stack.pop()
        stats["nodes"] += 1
        stats["depth"] = max(stats["depth"], depth)
        tag = node[0]
        if tag == "l":
            idx = node[1]
            if idx in leaf_info:
                stats["dup_leaf"] = True
            leaf_info[idx] = (rect, chosen)
        elif tag == "p":
            s = node[1]
            if s in used or not (0 <= s < len(specs)) or len(node) != 4:
                stats["branch_ok"] = False
            stack.append((node[2], rect, used | {s}, chosen | {(s, 0)}, depth + 1))
            stack.append((node[3], rect, used | {s}, chosen | {(s, 1)}, depth + 1))
        elif tag == "s":
            cells = [tuple(c[0]) for c in node[1]]
            if control == "overlap" and rect == (0, kR, 0, kL):
                cells = [(0, 23, 0, 30)] + cells[1:]
            if not cells_ok(rect, cells):
                stats["split_ok"] = False
            for c, child in zip(cells, [c[1] for c in node[1]]):
                stack.append((child, c, used, chosen, depth + 1))
        else:
            stats["branch_ok"] = False
    check("tree: every split's cells lie in the parent, are pairwise disjoint and cover it", stats["split_ok"])
    check("tree: no spec branched twice on a path", stats["branch_ok"])
    check("tree: leaves 0..12662, each exactly once; 24,898 nodes, depth 26 (the paper's table)",
          not stats["dup_leaf"] and sorted(leaf_info) == list(range(12663)) and stats["nodes"] == 24898
          and stats["depth"] == 26, "%d nodes, depth %d" % (stats["nodes"], stats["depth"]))

    # 5. the leaves
    PM, PP, P0 = 0, 1, 2
    R = [None] + [2 + i for i in range(1, kR + 1)]          # R[i] = label of R_i
    L = [None] + [2 + kR + j for j in range(1, kL + 1)]     # L[j] = label of L_j
    files = sorted(f for f in os.listdir(d) if f.startswith("multipliers-"))
    seen = set()
    nrows = 0
    worst = Q(0)
    minB = None
    bad_leaves = []
    checks_done = 0
    for fn in files:
        with gzip.open(os.path.join(d, fn), "rb") as fh:
            for line in fh:
                idx, D, mult = json.loads(line)
                if control == "mult" and idx == 0:
                    # leaf 0 has B = 1 exactly from its one escape row (4004); that multiplier one unit lower
                    mult = [list(x) for x in mult]
                    for x in mult:
                        if x[0] >= 3970:
                            x[1] -= 1
                if control == "esc" and idx == 0:
                    mult = [list(x) for x in mult] + [[3970 + 2 * 1462 + 1, 1]]
                seen.add(idx)
                rect, chosen = leaf_info[idx]
                good = D > 0 and all(c > 0 for _, c in mult)
                f = [[Q(0), Q(0)] for _ in range(NL)]
                Bsum = Q(0)
                for k, c in mult:
                    nrows += 1
                    mu = Q(c, D)
                    if 0 <= k < 3969:
                        lam, eta = divmod(k, NL)
                        if lam == eta:
                            good = False
                            continue
                        u = U[lam]
                        f[lam][0] += mu * u[0]
                        f[lam][1] += mu * u[1]
                        f[eta][0] -= mu * u[0]
                        f[eta][1] -= mu * u[1]
                    elif k == 3969:
                        f[PP][0] += mu
                        f[PM][0] -= mu
                    elif 3970 <= k < 3970 + 2 * len(specs):
                        s, a = divmod(k - 3970, 2)
                        if (s, a) not in chosen:
                            good = False
                            continue
                        if (s, a) not in ESC:
                            ESC[(s, a)] = escape_row(s, a)
                        coef, hj = ESC[(s, a)]
                        for lab, vec in coef.items():
                            f[lab][0] += mu * vec[0]
                            f[lab][1] += mu * vec[1]
                        Bsum += mu * hj
                    else:
                        good = False
                minB = Bsum if minB is None else min(minB, Bsum)
                good = good and Bsum >= 1
                good = good and sum(v[0] for v in f) == 0 and sum(v[1] for v in f) == 0
                den = 1
                for v in f:
                    den = lcm(den, v[0].denominator, v[1].denominator)
                F = [(int(v[0] * den), int(v[1] * den)) for v in f]
                lim = den * den
                p_lo, p_hi, q_lo, q_hi = rect
                big = 0
                # family 1: {L_a..L_b}, q_lo+1 <= a <= b <= kL
                for a in range(q_lo + 1, kL + 1):
                    x = y = 0
                    for b in range(a, kL + 1):
                        x += F[L[b]][0]
                        y += F[L[b]][1]
                        big = max(big, x * x + y * y)
                        checks_done += 1
                # family 2: {P-, L_j..L_kL}, 1 <= j <= q_hi + 1
                x, y = F[PM]
                tail = [(0, 0)] * (kL + 2)
                tx = ty = 0
                for j in range(kL, 0, -1):
                    tx += F[L[j]][0]
                    ty += F[L[j]][1]
                    tail[j] = (tx, ty)
                for j in range(1, q_hi + 2):
                    sx, sy = x + tail[j][0], y + tail[j][1]
                    big = max(big, sx * sx + sy * sy)
                    checks_done += 1
                # family 3: {P+, R_1..R_s}, p_lo <= s <= kR
                x, y = F[PP]
                for s in range(0, kR + 1):
                    if s >= 1:
                        x += F[R[s]][0]
                        y += F[R[s]][1]
                    if s >= p_lo:
                        big = max(big, x * x + y * y)
                        checks_done += 1
                # family 4: {R_{a+1}..R_b}, 0 <= a <= b <= p_hi (a = b: empty)
                for a in range(0, p_hi + 1):
                    x = y = 0
                    for b in range(a + 1, p_hi + 1):
                        x += F[R[b]][0]
                        y += F[R[b]][1]
                        big = max(big, x * x + y * y)
                        checks_done += 1
                ratio = Q(big, lim)
                worst = max(worst, ratio)
                good = good and big <= lim
                if not good:
                    bad_leaves.append(idx)
                    if control:
                        check("leaf %d: (4.7) and the row rules" % idx, False)
    check("multiplier files: every leaf 0..12662 has exactly one line", seen == set(range(12663)),
          "%d rows in all (the manifest states 6,778,020)" % nrows)
    check("every leaf: D > 0, every multiplier > 0, rows valid, escape rows only from its path, B >= 1, "
          "sum f = 0, and |sum_I f|^2 <= 1 for every I in F(R)", not bad_leaves and nrows == 6778020,
          "%d prefix checks; largest |sum_I f|^2 = %.12f; least B = %.12f" % (checks_done, float(worst), float(minB))
          if not bad_leaves else "failing leaves: %s" % bad_leaves[:10])
    print("A: largest |sum_I f|^2 over all leaves and families: %s" % worst)
    print("A: least B over all leaves: %s" % minB)
    print("A: run time %.1f s" % (time.time() - t0))
    print("A: VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
