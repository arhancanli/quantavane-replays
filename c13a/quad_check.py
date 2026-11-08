"""Replay of the 13a upper-bound row [K2026] proposed in PR #314: the quadrilateral K of Keller (arXiv 2609.21968v1,
Section 4, equation (4.1) and the table after it) and its area 7170601298323360123534337/29109471743520000000000000.

Exact rational arithmetic (Python fractions) on the side data read from the paper's LaTeX source: the outward normals
n(sigma) = ((1 - sigma^2)/(1 + sigma^2), 2 sigma/(1 + sigma^2)), the vertices V1..V4 as intersections of consecutive
side lines (lines 0, 1 first), convexity and orientation, the shoelace area, |V2 - V4|^2 > 1 (a unit segment inside
K), and the normal relations a_j0 n_0 + a_j1 n_1 + n_j = 0 (j = 2, 3) with positive coefficients, as the paper states.
Not checked: that K covers every unit arc (the paper's 12,663-leaf certificate; its encoding of contact directions and
orientations is defined only in the author's code, which we do not read).

Our own code, standard library only. Usage: python3 -I quad_check.py CERTS_ROOT [--controls]
"""
import hashlib
import os
import re
import sys
from fractions import Fraction as Q


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def read_data(root):
    with open(os.path.join(root, "k2026/sections/upper.tex"), encoding="utf-8") as fh:
        up = fh.read()
    with open(os.path.join(root, "k2026/main.tex"), encoding="utf-8") as fh:
        mn = fh.read()
    sig = re.search(r"\\sigma_i&([^\\]*)\\\\", up).group(1).split("&")
    hh = re.search(r"\n h_i&([^\n]*)\n", up).group(1).split("&")
    sigma = [Q(x.strip()) for x in sig]
    h = [Q(x.strip()) for x in hh]
    frac = re.search(r"\\newcommand\{\\uboundvalue\}\{\\frac\{(\d+)\}\{(\d+)\}\}", mn)
    dec = re.search(r"\\newcommand\{\\ubounddecimal\}\{([0-9.]+)\\ldots\}", mn).group(1)
    return sigma, h, Q(int(frac.group(1)), int(frac.group(2))), dec, up


def normal(s):
    return ((1 - s * s) / (1 + s * s), 2 * s / (1 + s * s))


def meet(n1, h1, n2, h2):
    det = n1[0] * n2[1] - n1[1] * n2[0]
    if det == 0:
        raise ValueError("parallel side lines")
    return ((h1 * n2[1] - h2 * n1[1]) / det, (n1[0] * h2 - n2[0] * h1) / det)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def analyse(sigma, h, out, tag=""):
    n = [normal(s) for s in sigma]
    unit = all(dot(v, v) == 1 for v in n)
    V = [meet(n[0], h[0], n[1], h[1]), meet(n[1], h[1], n[2], h[2]), meet(n[2], h[2], n[3], h[3]),
         meet(n[3], h[3], n[0], h[0])]
    # each vertex satisfies all four inequalities, with equality exactly on its two lines
    on = [{i for i in range(4) if dot(n[i], v) == h[i]} for v in V]
    inside = all(dot(n[i], v) <= h[i] for v in V for i in range(4))
    exact_lines = on == [{0, 1}, {1, 2}, {2, 3}, {3, 0}]
    cross = []
    for k in range(4):
        a, b, c = V[k], V[(k + 1) % 4], V[(k + 2) % 4]
        cross.append((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]))
    convex = all(x > 0 for x in cross)
    area = sum(V[k][0] * V[(k + 1) % 4][1] - V[(k + 1) % 4][0] * V[k][1] for k in range(4)) / 2
    d24 = (V[1][0] - V[3][0]) ** 2 + (V[1][1] - V[3][1]) ** 2
    d13 = (V[0][0] - V[2][0]) ** 2 + (V[0][1] - V[2][1]) ** 2
    rel = {}
    for j in (2, 3):
        # a0 n0 + a1 n1 = -n_j
        det = n[0][0] * n[1][1] - n[0][1] * n[1][0]
        a0 = (-n[j][0] * n[1][1] + n[j][1] * n[1][0]) / det
        a1 = (-n[0][0] * n[j][1] + n[0][1] * n[j][0]) / det
        ok = a0 * n[0][0] + a1 * n[1][0] + n[j][0] == 0 and a0 * n[0][1] + a1 * n[1][1] + n[j][1] == 0
        rel[j] = (a0, a1, ok)
    out("%snormals are unit vectors: %s" % (tag, unit))
    for k, v in enumerate(V):
        out("%sV%d = (%s, %s)" % (tag, k + 1, v[0], v[1]))
    out("%severy vertex satisfies all four side inequalities: %s; each lies on exactly its two side lines: %s"
        % (tag, inside, exact_lines))
    out("%sturns at the four vertices all positive (convex, counterclockwise): %s" % (tag, convex))
    out("%sarea (shoelace) = %s = %s..." % (tag, area, dec_trunc(area, 25)))
    out("%s|V2 - V4|^2 = %s (> 1: %s); |V1 - V3|^2 = %s" % (tag, d24, d24 > 1, d13))
    for j in (2, 3):
        a0, a1, ok = rel[j]
        out("%snormal relation j = %d: a_j0 = %s, a_j1 = %s, exact: %s, both positive: %s"
            % (tag, j, a0, a1, ok, a0 > 0 and a1 > 0))
    good = unit and inside and exact_lines and convex and d24 > 1 and all(r[2] and r[0] > 0 and r[1] > 0
                                                                         for r in rel.values())
    return good, area


def dec_trunc(q, digits):
    s = "-" if q < 0 else ""
    q = abs(q)
    ip = q.numerator // q.denominator
    frac = q - ip
    ds = []
    for _ in range(digits):
        frac *= 10
        d = frac.numerator // frac.denominator
        ds.append(str(d))
        frac -= d
    return "%s%d.%s" % (s, ip, "".join(ds))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    controls = "--controls" in sys.argv

    def out(s):
        print(s)

    out("13a upper-bound row [K2026] (PR #314): C_13a <= 7170601298323360123534337/29109471743520000000000000*")
    for rel in ("k2026/main.tex", "k2026/sections/upper.tex"):
        out("%s sha256 %s" % (rel, sha256(os.path.join(root, rel))))
    sigma, h, stated, dec, up = read_data(root)
    out("side data read from the table: sigma = %s; h = %s" % ([str(s) for s in sigma], [str(x) for x in h]))
    ours = ([Q(-631, 58), Q(-1), Q(7, 16), Q(37, 56)], [Q(0), Q(0), Q(130783, 250000), Q(446199, 1000000)])
    tr = (sigma, h) == ours
    out("equal to our transcription of the table: %s" % tr)
    out("stated area: %s; stated decimal: %s..." % (stated, dec))
    good, area = analyse(sigma, h, out)
    eq = area == stated
    out("the shoelace area equals the stated fraction exactly: %s" % eq)
    tr_ok = dec_trunc(area, 10) == dec
    out("the stated decimal %s... is the area truncated at 10 decimals: %s" % (dec, tr_ok))
    readme = Q(24633, 100000)
    out("the README cell 0.24633* of PR #314 is %s the area (by %s...)" % (
        "below" if readme < area else "not below", dec_trunc(area - readme, 12)))
    status = tr and good and eq and tr_ok
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        # (i) h_2 changed in its last digit: the area moves off the stated fraction
        _, a2 = analyse(sigma, [h[0], h[1], h[2] + Q(1, 250000), h[3]], lambda s: None)
        c1 = a2 != stated
        out("[ctl] h_2 = 130784/250000: area %s..., no longer the stated fraction: %s" % (dec_trunc(a2, 12), c1))
        ctl &= c1
        # (ii) sigma_3 = 37/57: the vertices move and the area changes
        _, a3 = analyse([sigma[0], sigma[1], sigma[2], Q(37, 57)], h, lambda s: None)
        c2 = a3 != stated
        out("[ctl] sigma_3 = 37/57: area %s..., refused: %s" % (dec_trunc(a3, 12), c2))
        ctl &= c2
        # (iii) a decimal one unit higher in the 10th place is refused as the truncation
        c3 = dec_trunc(area, 10) != "0.2463322373"
        out("[ctl] 0.2463322373 refused as the truncation: %s" % c3)
        ctl &= c3
        # (iv) swapping two side lines breaks the vertex structure
        g4, _ = analyse([sigma[1], sigma[0], sigma[2], sigma[3]], [h[1], h[0], h[2], h[3]], lambda s: None)
        out("[ctl] lines 0 and 1 swapped: refused as the quadrilateral of the paper: %s" % (not g4))
        ctl &= not g4
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT: %s" % ("PASS (the area and the stated structure of K; the cover property is not checked)"
                        if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
