#!/usr/bin/env python3
"""5a, upper bound [COT26] C_5a <= 0.94304, program B (written separately from program A): the same WT-triple
conditions (Definition 2.2 of arXiv 2610.36207v1) checked with a generic exact step-function engine, no index algebra.

Each w_i and t_i is built as a list of pieces (left, right, value) straight from the paper's formulas (Section 6.1):
w_i = sum_j g_j 1_[j delta, (j+1) delta), t_i(x) = sum_j h_(3N-1-j) 1_[j delta, (j+1) delta)(x + upsilon), so the t-piece
for index j is [j delta - upsilon, (j+1) delta - upsilon). Then:
  - norms: ||f||_2^2 = sum (right - left) value^2 over pieces;
  - (w * t)(x) = int w(y) t(x - y) dy by merging the sorted breakpoints of w(y) and of y -> t(x - y) (the piece [a, b) of t
    becomes y in (x - b, x - a]), exact Fractions;
  - (vi) on [0, sigma]: at every knot j delta and also at every midpoint (j + 1/2) delta, where linearity between knots is
    checked (the midpoint value must be the average of its neighbours), so the knot minimum is the minimum;
  - (i) t_i(x) = 1 evaluated at the midpoints of [0, upsilon) on the delta grid and at 0; (ii) int w_i = 1; w, t >= 0.
alpha = E||t||^2 - 2 upsilon, beta = E||w||^2; 2 sqrt(alpha beta), 4 (1 + upsilon sqrt(beta/alpha))^2 and sqrt(alpha beta)
in 50-digit decimals moved outward by 10^-40. Red controls (--control): the claim 0.94303; one atom's probability
moved to another (p_0 -> p_0 + 1/1000, p_1 -> p_1 - 1/1000; (vi) or the bound must fail); window height 100 of atom 0
raised by 1/1000 (int w = 1 fails). (A first version doubled height 5 of atom 2, which is 0, so that control changed
nothing and was not refused; replaced.) Standard library only.
"""
import os
import sys
from fractions import Fraction as F
from decimal import Decimal as D, getcontext

getcontext().prec = 50
PAD = D(10) ** -40
ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else os.path.join(os.path.dirname(os.path.abspath(__file__)), "certs")
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None


def read_rows(path):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append([F(tok) for tok in line.split(",")])
    return out


def conv_at(wp, tp, x):
    """int w(y) t(x - y) dy for step functions given as sorted piece lists"""
    # y -> t(x - y): piece [a, b) of t maps to y in (x - b, x - a]; sort these by left end
    tq = sorted(((x - b, x - a, v) for (a, b, v) in tp), key=lambda r: r[0])
    i = j = 0; acc = F(0)
    while i < len(wp) and j < len(tq):
        a1, b1, v1 = wp[i]; a2, b2, v2 = tq[j]
        lo = max(a1, a2); hi = min(b1, b2)
        if hi > lo:
            acc += (hi - lo) * v1 * v2
        if b1 <= b2: i += 1
        else: j += 1
    return acc


def value_at(pieces, x):
    for a, b, v in pieces:
        if a <= x < b:
            return v
    return F(0)


def run(d, b2, k0max, row):
    G = read_rows(os.path.join(ROOT, d, "G.txt")); H = read_rows(os.path.join(ROOT, d, "H.txt"))
    p = read_rows(os.path.join(ROOT, d, "p.txt"))[0]
    ups = F(open(os.path.join(ROOT, d, "upsilon.txt")).read().strip())
    if CONTROL == "prob" and d == "asymptotic":
        p[0] += F(1, 1000); p[1] -= F(1, 1000)
    if CONTROL == "height" and d == "asymptotic":
        G[0][100] += F(1, 1000)
    K = len(p); N = len(G[0]); dl = ups / N; sig = 2 * ups
    W = [[(j * dl, (j + 1) * dl, g[j]) for j in range(N)] for g in G]
    T = [[(j * dl - ups, (j + 1) * dl - ups, h[3 * N - 1 - j]) for j in range(3 * N)] for h in H]
    shapes = all(len(h) == 3 * N for h in H) and len(G) == len(H) == K and sum(p) == 1 and min(p) >= 0
    c1 = all(value_at(T[i], m) == 1 for i in range(K) for m in [F(0)] + [(k + F(1, 2)) * dl for k in range(N)])
    c2 = all(sum((b - a) * v for a, b, v in W[i]) == 1 for i in range(K))
    strong = all(v >= 0 for P in W + T for _, _, v in P)
    n2 = lambda P: sum((b - a) * v * v for a, b, v in P)
    alpha = sum(p[i] * n2(T[i]) for i in range(K)) - sig
    beta = sum(p[i] * n2(W[i]) for i in range(K))
    knots = []
    for j in range(2 * N + 1):
        knots.append(sum(p[i] * conv_at(W[i], T[i], j * dl) for i in range(K) if p[i] != 0))
    linear = True
    for j in range(2 * N):
        mid = sum(p[i] * conv_at(W[i], T[i], (j + F(1, 2)) * dl) for i in range(K) if p[i] != 0)
        if mid != (knots[j] + knots[j + 1]) / 2:
            linear = False; break
    mn = min(knots); c6 = mn >= 1
    print("== %s: upsilon %s, K %d, N %d; shapes and p: %s; (i): %s; (ii): %s; strong: %s" % (d, ups, K, N, shapes, c1, c2, strong))
    print("   (vi): E(w * t) linear between knots (all %d midpoints): %s; min at knots %s (>= 1: %s)"
          % (2 * N, linear, format(D(mn.numerator) / D(mn.denominator), ".15f"), c6))
    A = D(alpha.numerator) / D(alpha.denominator); B = D(beta.numerator) / D(beta.denominator)
    s = (A * B).sqrt()
    bval = 2 * s + PAD
    k0 = 4 * (1 + ups.numerator / D(ups.denominator) * (B / A).sqrt()) ** 2 + PAD
    print("   alpha %s, beta %s; 2 sqrt(alpha beta) <= %s (<= %s: %s); 4(1 + upsilon sqrt(beta/alpha))^2 <= %s (<= %d: %s)"
          % (format(A, ".20f"), format(B, ".20f"), format(bval, ".12f"), b2, bval <= D(b2), format(k0, ".6f"), k0max, k0 <= k0max))
    ok = shapes and c1 and c2 and strong and linear and c6 and bval <= D(b2) and k0 <= k0max
    if row is not None:
        print("   the row: sqrt(alpha beta) <= %s (<= %s: %s)" % (format(s + PAD, ".12f"), row, s + PAD <= D(row)))
        ok = ok and s + PAD <= D(row)
    print("   key %s %s %s" % (d, format(A, ".25f")[:22], format(B, ".25f")[:22]))
    return ok


def main():
    row = "0.94303" if CONTROL == "claim" else "0.94304"
    ok = run("asymptotic", "1.88607", 483, row)
    if CONTROL is None:
        ok = run("explicit", "1.8898801", 30, None) and ok
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
