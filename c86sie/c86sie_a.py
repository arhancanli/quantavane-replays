#!/usr/bin/env python3
"""86a, lower bound [Sie45] ~ 1.73361: Siegel's bound e (1 + 1/nu)^(-nu), nu the positive root of
    f(nu) = (1 + nu) log(1 + 1/nu) + log(nu)/(1 + nu) - 1 = 0
(as stated in Orloski, Talebizadeh Sardari and Smith, arXiv 2401.03252v3, Section 1.1 and Corollary 1.4; Siegel's
paper itself was not read). Program A: exact rational intervals.
  - f < 0 on (0, 1/100] (written out in main()) and f > 0 for nu >= 1 (log(1 + x) >= x - x^2/2 gives
    (1 + nu) log(1 + 1/nu) >= 1 + 1/(2 nu) - 1/(2 nu^2) >= 1, and log nu >= 0); on [1/100, 1] the sign of f is
    enclosed on 2000 cells with a bound on |f'| (f' = log(1 + 1/nu) - 1/nu + 1/(nu (1 + nu)^2) - log(nu)/(1 + nu)^2,
    each term bounded by its value at the cell's left end), and the signs must read -1 ... -1, then undecided cells,
    then +1 ... +1, with f' > 0 on the undecided span: so f has exactly one positive root nu_0;
  - nu_0 by bisection to width 2^-110 inside the undecided cells; the bound B(nu) = e (1 + 1/nu)^(-nu)
    = exp(1 - nu log(1 + 1/nu)) is enclosed on the final bracket from both ends plus 10 * width (|B'| < 10 there).
Logarithms by atanh series, exp by Taylor with halving, all in exact rationals with explicit tail bounds.
Red controls (--control): the page figure claimed as 1.73362; the equation with + replaced by - before log(nu)/(1+nu)
(a different root: the bound changes); the root taken as 1/3 (the bound at nu = 1/3 differs in the sixth digit).
"""
import sys
from fractions import Fraction as F

CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
EPS = F(1, 10 ** 45)


def ln_iv(x):
    """log of a positive rational, enclosed: x = 2^k y, y in [1, 2), log y = 2 atanh((y-1)/(y+1))"""
    k = 0
    while x >= 2: x /= 2; k += 1
    while x < 1: x *= 2; k -= 1
    def at2(u):
        s = F(0); p = u; j = 0
        while p > EPS:
            s += p / (2 * j + 1); p *= u * u; j += 1
        return 2 * s, 2 * (s + p / (1 - u * u))
    lo, hi = at2((x - 1) / (x + 1))
    l2lo, l2hi = at2(F(1, 3))
    return (lo + k * (l2lo if k >= 0 else l2hi), hi + k * (l2hi if k >= 0 else l2lo))


def exp_iv(x):
    """exp of a rational |x| <= 4: exp(x/2^m)^(2^m)"""
    m = 0; y = x
    while abs(y) > F(1, 8):
        y /= 2; m += 1
    s = F(0); t = F(1); n = 0
    while abs(t) > EPS / 100:
        s += t; n += 1; t = t * y / n
    lo, hi = s - 2 * abs(t), s + 2 * abs(t)
    for _ in range(m):
        lo, hi = lo * lo, hi * hi
    return lo, hi


def f_iv(nu, sign=1):
    a = ln_iv(1 + 1 / nu)
    b = ln_iv(nu)
    lo = (1 + nu) * a[0] + sign * (b[0] if sign > 0 else b[1]) / (1 + nu) - 1
    hi = (1 + nu) * a[1] + sign * (b[1] if sign > 0 else b[0]) / (1 + nu) - 1
    return lo, hi


def fprime_bound(a, b):
    """a crude bound for |f'| on [a, b] (0 < a < b <= 1): each term bounded by monotone pieces"""
    la = ln_iv(1 + 1 / a)[1]
    lnmin = ln_iv(a)[0]           # log nu >= log a (negative)
    return la + 1 / a + 1 / (a * (1 + a) ** 2) + abs(lnmin) / (1 + a) ** 2


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def main():
    sign = -1 if CONTROL == "sign" else 1
    ok = True
    # on (0, 1/100]: f(nu) = (1 + nu) log(1 + nu) - nu (2 + nu)/(1 + nu) log(1/nu) - 1 <= (1 + nu) nu
    # + 3 nu log(1/nu) - 1 <= 0.0101 + 0.1382 - 1 < 0 (nu log(1/nu) increases on (0, 1/e)); on [1, oo): f > 0 (above)
    small = F(1, 100)
    near0 = F(101, 100) * small + 3 * small * ln_iv(1 / small)[1] - 1 < 0
    n = 2000
    grid = [small + (1 - small) * F(i, n) for i in range(n + 1)]
    signs = []
    for i in range(n):
        a, b = grid[i], grid[i + 1]
        fa = f_iv(a, sign)
        L = fprime_bound(a, b)
        if fa[0] - L * (b - a) > 0: signs.append(1)
        elif fa[1] + L * (b - a) < 0: signs.append(-1)
        else: signs.append(0)
    # pattern: -1 ... -1, then zero or more 0 (the cells that may hold a root), then +1 ... +1 (a single sign change)
    i = 0
    while i < n and signs[i] == -1: i += 1
    j = i
    while j < n and signs[j] == 0: j += 1
    one_change = near0 and i > 0 and j < n and all(x == 1 for x in signs[j:])
    print("f < 0 on (0, 1/100]: %s; signs on %d cells of [1/100, 1]: %d cells -1, then %d undecided, then %d of which +1: %d;"
          " one sign change only: %s; f > 0 on [1, oo) by the series bound"
          % (near0, n, i, j - i, n - j, sum(1 for x in signs[j:] if x == 1), one_change))
    # f increasing on the undecided span [a, b] (so the root there is single): f' >= log(1 + 1/b) - 1/a + 1/(b (1 + b)^2)
    # - log(b)/(1 + b)^2 (each term at its smallest; log b < 0 there)
    a_, b_ = grid[i], grid[j]
    fp_lo = ln_iv(1 + 1 / b_)[0] - 1 / a_ + 1 / (b_ * (1 + b_) ** 2) - ln_iv(b_)[1] / (1 + b_) ** 2
    print("on the undecided span [%s, %s]: f' >= %s > 0: %s" % (dec(a_, 6), dec(b_, 6), dec(fp_lo, 6), fp_lo > 0))
    ok &= one_change and fp_lo > 0
    lo, hi = grid[i], grid[j]
    if CONTROL == "root":
        lo = hi = F(1, 3)
    else:
        for _ in range(110):
            m = (lo + hi) / 2
            fm = f_iv(m, sign)
            if fm[1] < 0: lo = m
            elif fm[0] > 0: hi = m
            else: break
    print("nu_0 in [%s, %s]" % (dec(lo, 32), dec(hi, 32)))
    def Bv(nu):
        l = ln_iv(1 + 1 / nu)
        e1 = exp_iv(1 - nu * l[1]); e2 = exp_iv(1 - nu * l[0])
        return e1[0], e2[1]
    b1, b2 = Bv(lo), Bv(hi)
    Blo = min(b1[0], b2[0]) - 10 * (hi - lo); Bhi = max(b1[1], b2[1]) + 10 * (hi - lo)
    print("e (1 + 1/nu_0)^(-nu_0) in [%s, %s]" % (dec(Blo, 30), dec(Bhi, 30)))
    printed = F("1.73362") if CONTROL == "claim" else F("1.73361")
    u = F(1, 10 ** 5)
    rd = abs(Blo - printed) <= u / 2 and abs(Bhi - printed) <= u / 2
    lit = Blo >= F("1.7336105") and Bhi < F("1.7336106")
    print("page ~ %s: correct rounding %s; the source's 1.7336105... (a truncation): %s" % (printed, rd, lit))
    ok &= rd and lit
    print("key %s" % dec(Blo, 25))
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
