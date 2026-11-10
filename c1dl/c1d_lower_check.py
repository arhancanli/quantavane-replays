#!/usr/bin/env python3
"""Replay of the lower bound C_1d >= 1.152* (teorth/optimizationproblems #270, vidigoat, update of 10 Nov 2026, second,
commit 3c99a12; the argument as the c1d_lower README states it at vidigoat/optimization-records @ a5fcdef). Written
from that README's statement of the proof alone: check_c1d_lower.py and sanity_random.py were neither read nor run.

The argument (re-derived by hand; the computation below checks its one finite inequality).
f integrable on [-1/4, 1/4], int f = 1, g = f*f on [-1/2, 1/2], int g = 1, M = ||g||_oo >= 1. F(z) = int f e^{-2 pi i z x}
is entire and g^ = F^2. h = (M - g) 1_[-1/2,1/2] has 0 <= h <= 2M, int h = M - 1 and h^(z) = M S(z) - F(z)^2 with
S(z) = sin(pi z)/(pi z).
 (1) For z = x + iy, |e^{-2 pi i z t}| = e^{2 pi y t}, a MONOTONE weight in t. Over 0 <= h <= 2M with int h = M - 1 the
     integral int h e^{2 pi y t} dt is largest for h = 2M on ONE end interval of length (M - 1)/(2M) (bathtub), which
     gives |h^(z)| <= B_M(|y|) = (M/(pi |y|)) (e^{pi |y|} - e^{pi |y|/M}), B_M(0) = M - 1.
     (The README writes the intermediate weight as e^{2 pi |y| |t|}; for that symmetric weight the maximiser splits
     the mass between both ends and the maximum is larger, B'_M(y) = (2M/(pi y)) (e^{pi y} - e^{pi y (M+1)/(2M)}),
     with which the contour inequality FAILS at M = 1.152: red control 2. The bound B_M itself is right, by (1).)
 (2) B_M(y)/M = phi_M(pi y) with phi_M(u) = int_{1/M}^{1} e^{s u} ds = sum_k u^k (1 - M^{-(k+1)})/(k+1)!, which
     increases in u >= 0 and in M; so the inequality at M0 = 1.152 gives it for every M in [1, M0].
 (3) |S(x + iy)|^2 = (sin^2(pi x) + sinh^2(pi y)) / (pi^2 (x^2 + y^2)).
 (4) Rouche on R = [0.87, 1.23] x [-0.34, 0.34]: if |S(z)| > phi_M0(pi |Im z|) on dR then F^2 and M S have equally
     many zeros in R; M S has exactly one (z = 1, simple), every zero of F^2 has even multiplicity: contradiction,
     so M > 1.152 for every admissible f.

What is computed here: the strict inequality |S(z)|^2 > phi_M0(pi |y|)^2 on all four sides of dR, M0 = 144/125.
Both sides are even in y, so y >= 0 suffices. On a vertical side (x fixed) a cell [c, d] of y is accepted when
    (sin^2(pi x) + sinh^2(pi c)) / (pi^2 (x^2 + d^2)) > phi_M0(pi d)^2,
and on a horizontal side (y = 0.34) a cell [a, b] of x when
    (min_{[a,b]} sin^2(pi x) + sinh^2(0.34 pi)) / (pi^2 (max(a^2, b^2) + 0.34^2)) > phi_M0(0.34 pi)^2,
using |sin(pi x)| = sin(pi |1 - x|) (|1 - x| <= 0.23), increasing in |1 - x|. Two implementations that share no code:
  A. exact rationals (Fraction), every series bounded above and below with an explicit remainder, pi from Machin's
     formula (the formula confirmed by Gaussian integers), adaptive bisection;
  B. Python's decimal at 60 digits with directed rounding (ROUND_FLOOR / ROUND_CEILING contexts), exp from the
     correctly rounded Decimal.exp widened by one unit, pi from Gauss's arctangent formula, a fixed uniform grid.
Red controls (each must be refused): 1. M0 = 1.154 (the corner z = 1.23 + 0.34i then violates the inequality,
shown rigorously); 2. the two-sided remainder bound B'_M of the README's intermediate step, at M0 = 1.152; 3. the
point z = 1 inside R, where S = 0 (the evaluator must refuse it); 4. (floats, a sanity check of the lemma, not part of
the proof) |h^(z)| <= B_M(|y|) for random signed step functions f at random z.

Usage: python3 c1d_lower_check.py [--controls]     (standard library only; under a minute)
"""
import cmath
import decimal
import math
import random
import sys
import time
from fractions import Fraction as F

M0 = F(144, 125)          # 1.152
X0, X1, Y1 = F(87, 100), F(123, 100), F(34, 100)
P = 200                   # outward rounding grid 2^-P for method A
ONE = 1 << P


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ------------------------------------------------------------------ method A: exact rationals, outward rounded
def rd(x):
    return F((x.numerator * ONE) // x.denominator, ONE)


def ru(x):
    return F(-((-x.numerator * ONE) // x.denominator), ONE)


def atan_inv_bounds(n, terms=60):
    """arctan(1/n) between two consecutive partial sums of the alternating series (terms decrease)."""
    s, prev = F(0), F(0)
    for k in range(terms + 1):
        prev = s
        s += F((-1) ** k, (2 * k + 1) * n ** (2 * k + 1))
    return min(s, prev), max(s, prev)


def pi_bounds_machin():
    # Machin: pi/4 = 4 arctan(1/5) - arctan(1/239); confirmed: (5+i)^4 (239-i) is a positive real multiple of 1+i
    z = (5, 1)
    w = (1, 0)
    for _ in range(4):
        w = (w[0] * z[0] - w[1] * z[1], w[0] * z[1] + w[1] * z[0])
    w = (w[0] * 239 + w[1] * 1, -w[0] * 1 + w[1] * 239)
    need(w[0] == w[1] and w[0] > 0, "Machin's formula: (5+i)^4 (239-i) is not a positive multiple of 1+i")
    a_lo, a_hi = atan_inv_bounds(5)
    b_lo, b_hi = atan_inv_bounds(239)
    return rd(16 * a_lo - 4 * b_hi), ru(16 * a_hi - 4 * b_lo)


PI_LO, PI_HI = pi_bounds_machin()


def sin_bounds(t, terms=14):
    """0 <= t <= 1: sin t between consecutive partial sums of its alternating series."""
    need(0 <= t <= 1, "sin_bounds range")
    s, prev, term = F(0), F(0), t
    for k in range(terms + 1):
        prev = s
        s += term if k % 2 == 0 else -term
        term = term * t * t / ((2 * k + 2) * (2 * k + 3))
    return rd(min(s, prev)), ru(max(s, prev))


def sinh_bounds(t, terms=16):
    """0 <= t <= 2: partial sum, plus the tail bound t^(2N+3)/(2N+3)! / (1 - t^2/((2N+4)(2N+5)))."""
    need(0 <= t <= 2, "sinh_bounds range")
    s, term = F(0), t
    for k in range(terms + 1):
        s += term
        term = term * t * t / ((2 * k + 2) * (2 * k + 3))
    n = 2 * terms + 3                     # term is now t^n/n!
    tail = term / (1 - t * t / ((n + 1) * (n + 2)))
    return rd(s), ru(s + tail)


def phi_bounds(u, M, terms=40):
    """phi_M(u) = sum_k u^k (1 - M^-(k+1))/(k+1)!, 0 <= u <= 2: partial sum and tail <= u^(N+1)/(N+2)!/(1-u/(N+3))."""
    need(0 <= u <= 2, "phi_bounds range")
    s, upow, fact, minv = F(0), F(1), 1, 1 / M
    mp = minv
    for k in range(terms + 1):
        fact *= (k + 1)                   # (k+1)!
        s += upow * (1 - mp) / fact
        upow *= u
        mp *= minv
    # upow = u^(N+1); next factorial (N+2)!
    tail = upow / (fact * (terms + 2)) / (1 - u / (terms + 3))
    return rd(s), ru(s + tail)


def phi_bounds_twosided(u, M, terms=40):
    """The two-sided bound of red control 2, divided by M: (2/u)(e^u - e^{u (M+1)/(2M)}) =
    2 sum_k u^k (1 - c^(k+1))/(k+1)!, c = (M+1)/(2M)."""
    lo, hi = phi_bounds(u, 1 / ((M + 1) / (2 * M)), terms)   # phi with 1/M replaced by c
    return 2 * lo, 2 * hi


def lhs_lower(sin_arg_min, y_lo, x_absmax, y_hi):
    """Lower bound of |S|^2 = (sin^2 + sinh^2)/(pi^2 (x^2 + y^2)) over a box, given the smallest |1 - x| (so the
    smallest |sin(pi x)|), the smallest |y| and the largest |x|, |y|."""
    s_lo = sin_bounds(PI_LO * sin_arg_min)[0] if sin_arg_min > 0 else F(0)
    sh_lo = sinh_bounds(PI_LO * y_lo)[0] if y_lo > 0 else F(0)
    num = s_lo * s_lo + sh_lo * sh_lo
    den = PI_HI * PI_HI * (x_absmax * x_absmax + y_hi * y_hi)
    return rd(num / den)


def rhs_upper(y_hi, M, twosided=False):
    u = PI_HI * y_hi
    hi = (phi_bounds_twosided(u, M) if twosided else phi_bounds(u, M))[1]
    return ru(hi * hi)


def check_sides(M, twosided=False, n0=64, maxdepth=30):
    """Adaptive bisection over the three distinct sides (left x = 0.87, right x = 1.23, top/bottom |y| = 0.34).
    Returns (all_ok, boxes, smallest ratio lhs_lower/rhs_upper seen, list of failures at max depth)."""
    stats = {"boxes": 0, "minratio": None, "fail": [], "side": {}}

    def note(lo, hi, side):
        r = lo / hi
        if stats["minratio"] is None or r < stats["minratio"]:
            stats["minratio"] = r
        if side not in stats["side"] or r < stats["side"][side]:
            stats["side"][side] = r

    def vertical(x, c, d, depth):
        stats["boxes"] += 1
        lo = lhs_lower(abs(1 - x), c, x, d)
        hi = rhs_upper(d, M, twosided)
        if lo > hi:
            note(lo, hi, "x = %s" % x)
            return True
        if depth >= maxdepth:
            stats["fail"].append(("x", x, c, d))
            return False
        m = (c + d) / 2
        return vertical(x, c, m, depth + 1) & vertical(x, m, d, depth + 1)

    def horizontal(a, b, depth):
        stats["boxes"] += 1
        dmin = F(0) if a <= 1 <= b else min(abs(1 - a), abs(1 - b))
        lo = lhs_lower(dmin, Y1, max(abs(a), abs(b)), Y1)
        hi = rhs_upper(Y1, M, twosided)
        if lo > hi:
            note(lo, hi, "|y| = %s" % Y1)
            return True
        if depth >= maxdepth:
            stats["fail"].append(("y", Y1, a, b))
            return False
        m = (a + b) / 2
        return horizontal(a, m, depth + 1) & horizontal(m, b, depth + 1)

    ok = True
    for x in (X0, X1):
        for i in range(n0):
            ok &= vertical(x, Y1 * i / n0, Y1 * (i + 1) / n0, 0)
    for i in range(n0):
        ok &= horizontal(X0 + (X1 - X0) * i / n0, X0 + (X1 - X0) * (i + 1) / n0, 0)
    return ok, stats


# ------------------------------------------------------------------ method B: decimal, directed rounding
def method_b(M_str="1.152", n=4000):
    D = decimal.Decimal
    prec = 60
    cf = decimal.Context(prec=prec, rounding=decimal.ROUND_FLOOR)
    cc = decimal.Context(prec=prec, rounding=decimal.ROUND_CEILING)
    ch = decimal.Context(prec=prec, rounding=decimal.ROUND_HALF_EVEN)
    ulp = D(10) ** (-(prec - 5))          # generous: every quantity here is below 10^4

    def atan_inv_fixed(nn, scale, terms=80):
        """arctan(1/nn) * scale, enclosed in integers: each term floored and ceiled, and the two consecutive partial
        sums S_N, S_{N+1} of the alternating series bracket the value."""
        lo_a = hi_a = 0
        sums = []
        for k in range(terms + 1):
            d = (2 * k + 1) * nn ** (2 * k + 1)
            t_lo, t_hi = scale // d, -((-scale) // d)
            if k % 2 == 0:
                lo_a, hi_a = lo_a + t_lo, hi_a + t_hi
            else:
                lo_a, hi_a = lo_a - t_hi, hi_a - t_lo
            sums.append((lo_a, hi_a))
        return min(sums[-1][0], sums[-2][0]), max(sums[-1][1], sums[-2][1])

    # Gauss: pi = 48 arctan(1/18) + 32 arctan(1/57) - 20 arctan(1/239), in integers at the scale 10^75; the formula
    # holds because (18+i)^12 (57+i)^8 (239-i)^5 is a positive real multiple of 1+i
    g = (1, 0)
    for zz, e in (((18, 1), 12), ((57, 1), 8), ((239, -1), 5)):
        for _ in range(e):
            g = (g[0] * zz[0] - g[1] * zz[1], g[0] * zz[1] + g[1] * zz[0])
    need(g[0] == g[1] and g[0] > 0, "Gauss's formula: the Gaussian-integer product is not a positive multiple of 1+i")
    scale = 10 ** 75
    a = atan_inv_fixed(18, scale)
    b = atan_inv_fixed(57, scale)
    c = atan_inv_fixed(239, scale)
    pi_lo = cf.scaleb(D(48 * a[0] + 32 * b[0] - 20 * c[1]), -75)
    pi_hi = cc.scaleb(D(48 * a[1] + 32 * b[1] - 20 * c[0] + 1), -75)
    need(pi_lo < pi_hi and pi_hi - pi_lo < D(10) ** -50, "method B: pi enclosure")

    def exp_lo(x):
        return cf.subtract(ch.exp(x), cc.multiply(ulp, ch.exp(x)))

    def exp_hi(x):
        return cc.add(ch.exp(x), cc.multiply(ulp, ch.exp(x)))

    def sin_lo(t, terms=20):           # 0 <= t <= 1, alternating: S_N, S_{N+1} bracket; lower = min, rounded down
        s_a = D(0)
        term_lo = term_hi = t
        sums_lo = []
        for k in range(terms + 1):
            s_a = cf.add(s_a, term_lo if k % 2 == 0 else -term_hi)
            sums_lo.append(s_a)
            term_lo = cf.divide(cf.multiply(cf.multiply(term_lo, t), t), D((2 * k + 2) * (2 * k + 3)))
            term_hi = cc.divide(cc.multiply(cc.multiply(term_hi, t), t), D((2 * k + 2) * (2 * k + 3)))
        return min(sums_lo[-1], sums_lo[-2])

    def sinh_lo(t):
        return cf.divide(cf.subtract(exp_lo(t), exp_hi(-t)), D(2))

    M = D(M_str)

    def phi_hi(u):                     # (e^u - e^{u/M})/u for u >= 1/100; the series below that
        if u >= D("0.01"):
            return cc.divide(cc.subtract(exp_hi(u), exp_lo(cf.divide(u, M))), u)
        s = D(0)
        upow = D(1)
        fact = 1
        minv_k = cf.divide(D(1), M)          # M^-(k+1) rounded down, so 1 - M^-(k+1) is rounded up
        for k in range(30):
            fact *= (k + 1)
            s = cc.add(s, cc.divide(cc.multiply(upow, cc.subtract(D(1), minv_k)), D(fact)))
            upow = cc.multiply(upow, u)
            minv_k = cf.divide(minv_k, M)
        return cc.add(s, cc.multiply(upow, D(2)))   # tail <= 2 u^30 / 31! < 2 u^30 (u < 1/100)

    def lhs_lo(dmin, ylo, xmax, yhi):
        s = sin_lo(cf.multiply(pi_lo, dmin)) if dmin > 0 else D(0)
        sh = sinh_lo(cf.multiply(pi_lo, ylo)) if ylo > 0 else D(0)
        num = cf.add(cf.multiply(s, s), cf.multiply(sh, sh))
        den = cc.multiply(cc.multiply(pi_hi, pi_hi), cc.add(cc.multiply(xmax, xmax), cc.multiply(yhi, yhi)))
        return cf.divide(num, den)

    def rhs_hi(yhi):
        p = phi_hi(cc.multiply(pi_hi, yhi))
        return cc.multiply(p, p)

    x0, x1, y1 = D("0.87"), D("1.23"), D("0.34")
    worst = None
    fails = 0
    for x in (x0, x1):
        for i in range(n):
            c_ = y1 * i / n
            d_ = y1 * (i + 1) / n
            lo = lhs_lo(abs(1 - x), c_, x, d_)
            hi = rhs_hi(d_)
            r = lo / hi
            worst = r if worst is None or r < worst else worst
            fails += lo <= hi
    for i in range(n):
        a_ = x0 + (x1 - x0) * i / n
        b_ = x0 + (x1 - x0) * (i + 1) / n
        dmin = D(0) if a_ <= 1 <= b_ else min(abs(1 - a_), abs(1 - b_))
        lo = lhs_lo(dmin, y1, max(abs(a_), abs(b_)), y1)
        hi = rhs_hi(y1)
        r = lo / hi
        worst = r if worst is None or r < worst else worst
        fails += lo <= hi
    return fails == 0, 3 * n, worst, (pi_lo, pi_hi)


# ------------------------------------------------------------------ helpers for the report and the controls
def corner_violation(M):
    """Rigorous: at z = 1.23 + 0.34i, upper bound of |S|^2 < lower bound of phi_M(0.34 pi)^2 (the inequality fails)."""
    s_hi = sin_bounds(PI_HI * F(23, 100))[1]
    sh_hi = sinh_bounds(PI_HI * Y1)[1]
    lhs_hi = ru((s_hi * s_hi + sh_hi * sh_hi) / (PI_LO * PI_LO * (X1 * X1 + Y1 * Y1)))
    ph_lo = phi_bounds(PI_LO * Y1, M)[0]
    return lhs_hi < rd(ph_lo * ph_lo), float(lhs_hi), float(ph_lo * ph_lo)


def sanity_lemma(trials=300, seed=1152):
    """Floats only, a sanity check of |h^(z)| <= B_M(|y|): random signed step functions on [-1/4, 1/4]."""
    rng = random.Random(seed)
    worst = 0.0
    for _ in range(trials):
        n = rng.randint(1, 40)
        v = [rng.uniform(-1, 3) for _ in range(n)]
        s = sum(v) / (2 * n)
        if abs(s) < 1e-3:
            continue
        v = [x / s for x in v]                 # int f = 1
        w = 1 / (2 * n)
        conv = [0.0] * (2 * n - 1)
        for i in range(n):
            for j in range(n):
                conv[i + j] += v[i] * v[j]
        M = max(abs(c) * w for c in conv)      # f*f is linear between knots; its knot values are w * conv
        for _ in range(5):
            z = complex(rng.uniform(0.8, 1.3), rng.uniform(-0.4, 0.4))
            Fz = sum(v[j] * (cmath.exp(-2j * math.pi * z * (-0.25 + j * w)) - cmath.exp(-2j * math.pi * z * (-0.25 + (j + 1) * w)))
                     for j in range(n)) / (2j * math.pi * z)
            Sz = cmath.sin(math.pi * z) / (math.pi * z)
            hz = abs(M * Sz - Fz * Fz)
            y = abs(z.imag)
            B = (M / (math.pi * y)) * (math.exp(math.pi * y) - math.exp(math.pi * y / M)) if y > 1e-12 else M - 1
            worst = max(worst, hz / B if B > 0 else 0.0)
    return worst


def main():
    t0 = time.time()
    controls = "--controls" in sys.argv
    print("C_1d >= 1.152* (#270, vidigoat, 10 Nov 2026, second update; c1d_lower README @ a5fcdef): the contour")
    print("inequality |S(z)| > B_M(|Im z|)/M on the boundary of R = [0.87, 1.23] x [-0.34, 0.34], M = 1.152 = 144/125")
    print("  pi (Machin, method A): an enclosure of width %.1e" % float(PI_HI - PI_LO))
    ok_a, st = check_sides(M0)
    print("  method A (exact rationals, adaptive from 64 cells per side): %s; %d boxes"
          % ("HOLDS on all four sides" if ok_a else "FAILS", st["boxes"]))
    ok_f, stf = check_sides(M0, n0=2048, maxdepth=0)
    print("  method A on 2,048 equal cells per side (no bisection): %s; smallest certified ratio of the squares, "
          "|S|^2 over (B/M)^2, per side: %s" % ("HOLDS" if ok_f else "FAILS",
          ", ".join("%s: %.6f" % (k, float(v)) for k, v in sorted(stf["side"].items()))))
    ok_a = ok_a and ok_f
    ok_b, nb, worst_b, pib = method_b()
    print("  method B (decimal, 60 digits, directed rounding, Gauss's pi, uniform grid): %s; %d cells; smallest ratio %.6f"
          % ("HOLDS on all four sides" if ok_b else "FAILS", nb, float(worst_b)))
    # the README's figures, for comparison only: ratio^2 about 1.0107 on three sides, 1.0416 on the left side
    zeros = [k for k in range(-3, 4) if k != 0 and X0 <= k <= X1]
    print("  zeros of S in R: the non-zero integers in [0.87, 1.23] on the real axis: %s (simple)" % zeros)
    ok = ok_a and ok_b and zeros == [1]
    print("  so by Rouche every f with int f = 1 has ||f*f||_oo > 1.152, and C_1d >= 1.152: %s" % ok)
    if controls:
        print("Red controls (each must be refused):")
        okc, stc = check_sides(F(1154, 1000), maxdepth=6)
        viol, l_, r_ = corner_violation(F(1154, 1000))
        c1 = (not okc) and viol
        print("  1. M0 = 1.154: the check fails (%d boxes left at depth 6); at z = 1.23 + 0.34i, |S|^2 <= %.8f < %.8f <= "
              "phi^2 rigorously, so the inequality is false there: refused: %s" % (len(stc["fail"]), l_, r_, c1))
        ok2, st2 = check_sides(M0, twosided=True, maxdepth=6)
        c2 = not ok2
        print("  2. the two-sided remainder bound B'_M of the README's intermediate weight e^{2 pi |y| |t|}, at M = 1.152:"
              " fails on %d boxes at depth 6: refused: %s" % (len(st2["fail"]), c2))
        lo = lhs_lower(F(0), F(0), F(1), F(0))
        hi = rhs_upper(F(0), M0)
        c3 = not (lo > hi)
        print("  3. the point z = 1 (inside R, S(1) = 0): lhs lower bound %s, rhs %.6f: refused: %s" % (lo, float(hi), c3))
        w = sanity_lemma()
        c4 = w <= 1
        print("  4. (floats, sanity only) max |h^(z)| / B_M(|y|) over random signed step functions: %.4f <= 1: %s" % (w, c4))
        ok = ok and c1 and c2 and c3 and c4
    print("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("run time %.0f s" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
