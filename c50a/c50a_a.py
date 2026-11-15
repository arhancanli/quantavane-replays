#!/usr/bin/env python3
"""50a, program A: [GP19] product-state approximation ratio 0.498 for Quantum Max Cut (Gharibian and Parekh,
arXiv 1909.08846v1, Theorem 4 and its proof), recomputed from its definition (a value from its definition, not a
certificate).

The proof rounds the SDP vectors to Bloch vectors in R^r (r = 3 for the anti-ferromagnetic Heisenberg model) and
gets, on an edge with SDP inner product t, the ratio R_r(t) = (1 - F_r(t)) / (1 - r t), where (Lemma of Briet,
de Oliveira Filho and Vallentin, the paper's appendix)
    F_r(t) = (2/r) (Gamma((r+1)/2) / Gamma(r/2))^2  t  2F1(1/2, 1/2; r/2 + 1; t^2).
The ratio of the algorithm is the minimum of R_r over t in [-1, 1/r); the paper evaluates it on a grid of step 0.01
and prints 0.878 (r = 1, Goemans-Williamson), 0.649 (r = 2) and 0.498 (r = 3).

Route A: the hypergeometric series in exact rationals, 600 terms and a tail bound
(c_n = ((1/2)_n)^2 / ((r/2+1)_n n!) <= C(2n,n)/4^n * (1/2)_n/(r/2+1)_n, C(2n,n)/4^n <= 1/sqrt(3n)), pi by Machin's
formula; F_r is odd with positive Taylor coefficients, so increasing, and on a cell [a, b] with b < 1/r the ratio is
at least (1 - F_r(b)) / (1 - r a). The cells are bisected until every one clears the target; for t >= 0 one cell
suffices (R_r(t) >= 1 - F_r(1/r)). Every interval end is rounded outward to a multiple of 2^-200.
Usage: c50a_a.py [--r 2|3] [--target T] [--g-unsquared]
"""
import sys
from fractions import Fraction as Q
from math import floor, isqrt

PREC = 200
def dn(q):
    return Q(floor(q * (1 << PREC)), 1 << PREC)
def upq(q):
    return -dn(-q)

# pi by Machin: pi = 16 atan(1/5) - 4 atan(1/239), alternating series with remainder
def atan_inv(m, terms):
    s, lo, hi = Q(0), None, None
    for k in range(terms):
        s += Q((-1) ** k, (2 * k + 1) * m ** (2 * k + 1))
    rem = Q(1, (2 * terms + 1) * m ** (2 * terms + 1))
    return (s - rem, s + rem)
a5, a239 = atan_inv(5, 80), atan_inv(239, 30)
PI = (16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0])

def gfactor(r):
    """(2/r) (Gamma((r+1)/2)/Gamma(r/2))^2 as an interval: rational times pi or 1/pi"""
    # Gamma ratios for small r in closed form: r = 1: Gamma(1)/Gamma(1/2) = 1/sqrt(pi); r = 2: Gamma(3/2)/Gamma(1) =
    # sqrt(pi)/2; r = 3: Gamma(2)/Gamma(3/2) = 2/sqrt(pi); r = 4: Gamma(5/2)/Gamma(2) = 3 sqrt(pi)/4
    if r == 1:
        return (2 / PI[1], 2 / PI[0])                   # 2/pi
    if r == 2:
        return (PI[0] / 4, PI[1] / 4)                   # pi/4
    if r == 3:
        return (Q(8, 3) / PI[1], Q(8, 3) / PI[0])       # 8/(3 pi)
    if r == 4:
        return (Q(9, 32) * PI[0], Q(9, 32) * PI[1])     # (1/2)(9 pi / 16)
    raise SystemExit("r must be 1..4")

CN = {}
def coeffs(r, NT):
    """exact c_n = ((1/2)_n)^2 / ((r/2+1)_n n!), n < NT"""
    if (r, NT) not in CN:
        c, out = Q(1), []
        for n in range(NT):
            out.append(c)
            c = c * Q(2 * n + 1, 2) ** 2 / (Q(2 * n + r + 2, 2) * (n + 1))
        CN[(r, NT)] = out
    return CN[(r, NT)]

def tail_bound(r, NT):
    """sum_(n >= NT) c_n, from (1/2)_n/n! = C(2n,n)/4^n <= 1/sqrt(3n) and
    r = 3: (1/2)_n/(5/2)_n = (3/4)/((n+1/2)(n+3/2)) <= (3/4)/n^2, so c_n <= (1/(2 sqrt 3)) (3/2) n^(-5/2)... summed:
          sum_(n>=NT) c_n <= (3/4)/sqrt(3) * (2/3) (NT-1)^(-3/2) <= (29/100) / ((NT-1) isqrt(NT-1));
    r = 2: (1/2)_n/(2)_n = ((1/2)_n/n!)/(n+1), so c_n <= 1/(3 n^2) and the tail is <= 1/(3 (NT-1))."""
    m = NT - 1
    if r == 3:
        return Q(29, 100) / (m * isqrt(m))
    if r == 2:
        return Q(1, 3 * m)
    raise SystemExit("series route for r = 2, 3 only")

def F_upper(r, t, NT=600):
    """an upper bound for F_r(t), t rational in [-1, 1]"""
    z = t * t
    c = coeffs(r, NT)
    s, zn = Q(0), Q(1)
    for n in range(NT):
        s = dn(s + c[n] * zn)            # a lower bound of the partial sum (zn is exact or rounded down)
        zn = dn(zn * z)
    g = gfactor(r)
    if t >= 0:
        return upq(g[1] * t * (s + tail_bound(r, NT) + Q(NT, 1 << (PREC - 20))))   # rounding slack for the n terms
    return upq(-(g[0] * (-t) * s))       # F odd: F(t) = -F(|t|) <= -g_lo |t| S_lower

args = sys.argv[1:]
r, target, gfree = 3, Q("0.4987"), False
while args:
    k = args.pop(0)
    if k == "--r":
        r = int(args.pop(0))
    elif k == "--target":
        target = Q(args.pop(0))
    elif k == "--g-unsquared":
        gfree = True
    else:
        raise SystemExit("bad argument")
def qsqrt(q, upper):     # rational square root bounds via integer square roots at scale 2^400
    n = isqrt(floor(q * (1 << 800)))
    return Q(n + 1 if upper else n, 1 << 400)
if gfree:            # red control: (Gamma(2)/Gamma(3/2)) not squared, i.e. 8/(3 pi) replaced by 4/(3 sqrt(pi))
    gfactor = lambda rr: (Q(4, 3) / qsqrt(PI[1], True), Q(4, 3) / qsqrt(PI[0], False))
    print("CONTROL variant: (Gamma(2)/Gamma(3/2)) not squared")

print("program A (c50a_a.py): hypergeometric series in exact rationals, pi by Machin")
def dec(q, n=30):
    s = str(abs(q.numerator) * 10 ** n // q.denominator).rjust(n + 1, "0")
    return ("-" if q < 0 else "") + s[:-n] + "." + s[-n:]
print("pi in [%s, %s]" % (dec(PI[0]), dec(PI[1])))
if r == 1:
    # anchor: F_1(t) = (2/pi) arcsin t; Goemans-Williamson; not run through the series (its tail is too slow here)
    raise SystemExit("r = 1 is the anchor in program B")
# t >= 0: R >= 1 - F(1/r)
F_third = F_upper(r, Q(1, r))
low_pos = 1 - F_third
print("t in [0, 1/%d): R >= 1 - F(1/%d) >= %.6f" % (r, r, float(low_pos)))
ok = low_pos >= target
# t in [-1, 0]: cells, bisected
stack = [(Q(-1), Q(0))]
cells, worst, worst_cell = 0, None, None
depth_cap = 2 ** -24
memo = {}
def Fu(t):
    if t not in memo:
        memo[t] = F_upper(r, t)
    return memo[t]
# start from 2000 equal cells
stack = [(Q(-1) + Q(i, 2000), Q(-1) + Q(i + 1, 2000)) for i in range(2000)]
while stack:
    a, b = stack.pop()
    lb = (1 - Fu(b)) / (1 - r * a)
    if lb >= target:
        cells += 1
        if worst is None or lb < worst:
            worst, worst_cell = lb, (a, b)
        continue
    if b - a < depth_cap:
        ok = False
        print("cell [%.9f, %.9f] does not clear the target (lower bound %.9f)" % (float(a), float(b), float(lb)))
        break
    m = (a + b) / 2
    stack += [(a, m), (m, b)]
if not stack:
    print("t in [-1, 0]: %d cells, every lower bound >= %s (least %.9f, on [%.6f, %.6f])"
          % (cells, target, float(worst), float(worst_cell[0]), float(worst_cell[1])))
else:
    print("t in [-1, 0]: stopped at a cell of width 2^-24 that does not clear the target (%d cells cleared before it)" % cells)
ok = ok and stack == []
print("RESULT r = %d: min over [-1, 1/%d) of (1 - F_r(t))/(1 - %d t) >= %s: %s" % (r, r, r, target, ok))
print("VERDICT", "PASS" if ok else "REFUSED")
