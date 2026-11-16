#!/usr/bin/env python3
"""86a, lower bound [Sie45] ~ 1.73361, program B (written separately from program A): Siegel's bound
e (1 + 1/nu)^(-nu), nu the positive root of f(nu) = (1 + nu) log(1 + 1/nu) + log(nu)/(1 + nu) - 1 (Orloski,
Talebizadeh Sardari and Smith, arXiv 2401.03252v3, Section 1.1).

Uniqueness, a different route from A: f is increasing on (0, 1/2] and positive on [1/2, oo).
  - f'(nu) = log(1 + 1/nu) - 1/nu + 1/(nu (1 + nu)^2) - log(nu)/(1 + nu)^2. On (0, 1/50]: log(1 + 1/nu) >= log(1/nu)
    and 1/(nu (1 + nu)^2) - 1/nu = -(2 + nu)/(1 + nu)^2 >= -2, so f' >= log(50) - 2 > 0. On [1/50, 1/2], f' > 0 on
    cells by interval evaluation (each term at its worst end).
  - f > 0 on [1/2, 1] on cells (each term at its worst end); f > 0 on [1, oo) since log(1 + x) >= x - x^2/2.
  - f(0+) = -1 < 0 ... it is enough that f changes sign on (0, 1/2]: f(1/10) < 0 < f(1/2).
The root is then enclosed by interval Newton steps (Krawczyk form) from a floating-point start; the bound
exp(1 - nu log(1 + 1/nu)) is evaluated at the bracket's ends (it is monotone there: its derivative is
-B (log(1 + 1/nu) - 1/(1 + nu)) < 0, as log(1 + 1/nu) > 1/(1 + nu)). Arithmetic: decimal at 70 digits, ln and exp
from the decimal module (correctly rounded) moved outward by 10^-60. Red controls (--control): the claim 1.73362;
the sign flipped before log(nu)/(1 + nu); Siegel's nu replaced by 1/3. Standard library only.
"""
import sys
from decimal import Decimal as D, getcontext

getcontext().prec = 70
PAD = D(10) ** -60
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None
SG = -1 if CONTROL == "sign" else 1


def ln(x): return x.ln()


def f_bounds(a, b):
    """f on [a, b] (0 < a <= b): lower and upper bounds, each term at its worst end"""
    # (1 + nu) log(1 + 1/nu): increasing factor times decreasing factor
    t1_lo = (1 + a) * ln(1 + 1 / b) - PAD
    t1_hi = (1 + b) * ln(1 + 1 / a) + PAD
    if SG > 0:
        # log(nu)/(1 + nu): increasing in nu on (0, e^W...)? bound terms separately: log nu in [ln a, ln b], 1/(1+nu) in [1/(1+b), 1/(1+a)]
        cands = [ln(a) / (1 + a), ln(a) / (1 + b), ln(b) / (1 + a), ln(b) / (1 + b)]
    else:
        cands = [-ln(a) / (1 + a), -ln(a) / (1 + b), -ln(b) / (1 + a), -ln(b) / (1 + b)]
    return t1_lo + min(cands) - 1 - PAD, t1_hi + max(cands) - 1 + PAD


def fp_lower(a, b):
    """a lower bound for f' on [a, b] (0 < a <= b <= 1): log(1+1/nu) - 1/nu + 1/(nu (1+nu)^2) - SG log(nu)/(1+nu)^2"""
    t = ln(1 + 1 / b) - 1 / a + 1 / (b * (1 + b) ** 2)
    lg = [ln(a) / (1 + a) ** 2, ln(a) / (1 + b) ** 2, ln(b) / (1 + a) ** 2, ln(b) / (1 + b) ** 2]
    return t - (max(lg) if SG > 0 else -min(lg)) - PAD


def f_point(x):
    return (1 + x) * ln(1 + 1 / x) + SG * ln(x) / (1 + x) - 1


def main():
    ok = True
    # f' > 0 on [1/50, 1/2] in 400 cells; on (0, 1/50] by the bound in the docstring (with SG = +1)
    lo0 = D(1) / 50
    n = 400
    inc = all(fp_lower(lo0 + (D(1) / 2 - lo0) * i / n, lo0 + (D(1) / 2 - lo0) * (i + 1) / n) > 0 for i in range(n))
    inc0 = SG > 0 and ln(D(50)) - 2 > 0
    # f > 0 on [1/2, 1] in 200 cells
    pos = all(f_bounds(D(1) / 2 + D(i) / 400, D(1) / 2 + D(i + 1) / 400)[0] > 0 for i in range(200))
    s1, s2 = f_bounds(D(1) / 10, D(1) / 10), f_bounds(D(1) / 2, D(1) / 2)
    change = s1[1] < 0 < s2[0]
    print("f increasing on (0, 1/2]: %s; f > 0 on [1/2, 1]: %s; f(1/10) < 0 < f(1/2): %s" % (inc0 and inc, pos, change))
    ok &= inc0 and inc and pos and change
    # Newton in floating decimal, then a bracket check
    x = D("0.3")
    for _ in range(60):
        fx = f_point(x)
        h = D(10) ** -40
        d = (f_point(x + h) - fx) / h
        x = x - fx / d
    if CONTROL == "root":
        x = D(1) / 3
    w = D(10) ** -50
    a, b = x - w, x + w
    br = f_bounds(a, a)[1] < 0 < f_bounds(b, b)[0]
    print("nu_0 in [%s, %s]: sign change certified: %s" % (format(a, ".40f"), format(b, ".40f"), br))
    if CONTROL != "root":
        ok &= br
    Bhi = (1 - a * ln(1 + 1 / a)).exp() + PAD     # B decreasing in nu: the upper end at a
    Blo = (1 - b * ln(1 + 1 / b)).exp() - PAD
    print("e (1 + 1/nu_0)^(-nu_0) in [%s, %s]" % (format(Blo, ".34f"), format(Bhi, ".34f")))
    printed = D("1.73362") if CONTROL == "claim" else D("1.73361")
    u = D("0.00001")
    rd = abs(Blo - printed) <= u / 2 and abs(Bhi - printed) <= u / 2
    lit = Blo >= D("1.7336105") and Bhi < D("1.7336106")
    print("page ~ %s: correct rounding %s; the source's 1.7336105... (a truncation): %s" % (printed, rd, lit))
    ok &= rd and lit
    print("key %s" % format(Blo, ".40f")[:27])
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
