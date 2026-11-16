#!/usr/bin/env python3
"""Closed forms on pages that had no replay folder (6a, 11b, 35a, 88b, 23b, 46a, 56a, 61a), program B (written
separately from program A).

Decimal intervals at 60 significant digits: square roots and logarithms from the decimal module (each result is
correctly rounded, so one unit in the last place on each side encloses the exact value); log2 x as ln x / ln 2 with
outward rounding; rationals as exact numerator/denominator pairs compared by cross-multiplication. Each printed
decimal of the page is compared with the value through its digits: truncation (the printed digits are the leading
digits of the value), rounding to nearest, rounding up. Standard library only.
"""
import sys
from decimal import Decimal as D, getcontext, Context, ROUND_FLOOR, ROUND_CEILING

PREC = 60
getcontext().prec = PREC
DN = Context(prec=PREC, rounding=ROUND_FLOOR)
UP = Context(prec=PREC, rounding=ROUND_CEILING)


def widen(x):
    """a correctly rounded result x: return [x - ulp, x + ulp]"""
    return (DN.next_minus(x), UP.next_plus(x))


def sqrt_iv(n):
    return widen(D(n).sqrt())


def ln_iv(lo, hi):
    a = widen(lo.ln())[0]
    b = widen(hi.ln())[1]
    return (a, b)


def add(x, y): return (DN.add(x[0], y[0]), UP.add(x[1], y[1]))
def sub(x, y): return (DN.subtract(x[0], y[1]), UP.subtract(x[1], y[0]))


def divpos(x, y):           # x, y positive intervals
    return (DN.divide(x[0], y[1]), UP.divide(x[1], y[0]))


def log2_iv(x):
    return divpos(ln_iv(*x), ln_iv(D(2), D(2)))


class Rat:
    def __init__(self, n, d=1):
        from math import gcd
        g = gcd(n, d)
        self.n, self.d = n // g, d // g


def rel_rat(r, printed):
    """r: Rat. printed: string. Exact comparison through integers."""
    k = len(printed.split(".")[1]) if "." in printed else 0
    p = int(printed.replace(".", ""))            # printed * 10^k
    v10 = r.n * 10 ** k                          # value * 10^k * d
    out = []
    if p * r.d == v10:
        return "exact"
    if p * r.d <= v10 < (p + 1) * r.d:
        out.append("truncation")
    if (p - 1) * r.d < v10 <= p * r.d:
        out.append("rounded up")
    if abs(2 * (v10 - p * r.d)) <= r.d:
        out.append("correct rounding")
    return ", ".join(out) or "NONE"


def rel_iv(x, printed):
    k = len(printed.split(".")[1]) if "." in printed else 0
    d = D(printed)
    u = D(1).scaleb(-k)
    out = []
    if d <= x[0] and x[1] < d + u:
        out.append("truncation")
    if d - u < x[0] and x[1] <= d:
        out.append("rounded up")
    if d - u / 2 <= x[0] and x[1] <= d + u / 2:
        out.append("correct rounding")
    return ", ".join(out) or "NONE"


def table(ctrl):
    s5 = sqrt_iv(5)
    s2 = sqrt_iv(2)
    three, two, one = (D(3), D(3)), (D(2), D(2)), (D(1), D(1))
    psi = sub(three, s5) if ctrl != "6a-plus" else add(three, s5)
    psi = (DN.divide(psi[0], 2), UP.divide(psi[1], 2))
    T = [
        ("6a", "[AHS2022], [CL2022], [P2022], [S2022]", "(3 - sqrt 5)/2", psi, "0.381966", {"truncation"}),
        ("11b", "[KP2020]", "log2(3/2)", log2_iv((D("1.5"), D("1.5"))), "0.58497" if ctrl == "11b-print" else "0.58496",
         {"correct rounding"}),
        ("35a", "[AP24]", "log2(1 + sqrt 2)", log2_iv(add(one, s2)), "1.271", {"truncation"}),
        ("88b", "[Pol14b]", "4 - 28/157", Rat(4 * 157 - 28, 157), "3.821656", {"truncation"}),
        ("23b", "[AJM2025]", "17/92", Rat(17, 92), "0.18478", {"correct rounding"}),
        ("46a", "[Bo1991]", "4 - 2/15", Rat(4 * 13 - 2, 13) if ctrl == "46a-15" else Rat(4 * 15 - 2, 15), "3.86667",
         {"correct rounding"}),
        ("46a", "[Wo1995], [MVV1996]", "4 - 2/11", Rat(42, 11), "3.81818", {"correct rounding"}),
        ("46a", "[TVV1998]", "4 - 2/9", Rat(34, 9), "3.77778", {"correct rounding"}),
        ("46a", "[TV2000]", "4 - 2/7", Rat(26, 7), "3.71429", {"correct rounding"}),
        ("46a", "[Tao2003], [BG2011]", "10/3", Rat(10, 3), "3.33333", {"correct rounding"}),
        ("46a", "[Gut2016]", "13/4", Rat(13, 4), "3.25", {"exact"}),
        ("46a", "[Wan2022]", "3 + 3/13", Rat(42, 13), "3.23077", {"correct rounding"}),
        ("46a", "[WW2022]", "3 + 3/14", Rat(45, 14), "3.21429", {"correct rounding"}),
        ("46a", "[WW2024]", "22/7", Rat(22, 7), "3.142857", {"correct rounding"}),
        ("56a", "[LRS1999]", "1/2 - 1/(2^2 + 1)", Rat(5 - 2, 10), "0.3", {"exact"}),
        ("56a", "[KS2003]", "7/64", Rat(7, 64), "0.109375", {"exact"}),
    ]
    return T


def digits(v, k=30):
    if isinstance(v, Rat):
        q = v.n * 10 ** k // v.d
        s = str(q).rjust(k + 1, "0")
        return s[:-k] + "." + s[-k:]
    return str(v[0].quantize(D(1).scaleb(-k), rounding=ROUND_FLOOR))


def main():
    a = sys.argv[1:]
    ctrl = a[a.index("--control") + 1] if "--control" in a else None
    print("program B (decimal intervals, %d digits)%s" % (PREC, (" CONTROL " + ctrl) if ctrl else ""))
    # self-tests: ln(e) = 1 by the decimal exp, sqrt(5)^2 = 5, log2(2^10) = 10
    e = D(1).exp()
    li = ln_iv(DN.next_minus(e), UP.next_plus(e))
    s = sqrt_iv(5)
    l2 = log2_iv((D(1024), D(1024)))
    ok = li[0] <= 1 <= li[1] and DN.multiply(s[0], s[0]) <= 5 <= UP.multiply(s[1], s[1]) and l2[0] <= 10 <= l2[1]
    if not ok:
        print("SELF-TEST FAILED"); sys.exit(2)
    print("self-tests: ln e = 1, sqrt(5)^2 = 5, log2 1024 = 10: ok")
    nbad = 0
    for page, tag, expr, v, printed, want in table(ctrl):
        r = rel_rat(v, printed) if isinstance(v, Rat) else rel_iv(v, printed)
        good = bool(want & set(r.split(", ")))
        nbad += not good
        print("row %s %s: %s = %s...; page %s: %s  %s" % (page, tag, expr, digits(v), printed, r, "ok" if good else "MISMATCH"))
    # 61a: theta from lambda_1 = s(1 - s), s = 1/2 + theta, i.e. theta = sqrt(1/4 - lambda_1), squared and compared exactly
    for tag, lam_n, lam_d, th_n, th_d in (("[LRS1995]", 3, 16, 1, 4), ("[KRS2003]", 975, 4096, 7, 64)):
        # 1/4 - lam = theta^2  <=>  (lam_d - 4 lam_n) * th_d^2 = 4 lam_d th_n^2
        lhs = (lam_d - 4 * lam_n) * th_d ** 2 if not (ctrl == "61a-sign" and tag == "[KRS2003]") else (lam_d + 4 * lam_n) * th_d ** 2
        good = lhs == 4 * lam_d * th_n ** 2
        nbad += not good
        print("row 61a %s: lambda_1 >= %d/%d gives theta = %d/%d  %s" % (tag, lam_n, lam_d, th_n, th_d, "ok" if good else "MISMATCH"))
    print("VERDICT " + (("REFUSED" if nbad else "NOT REFUSED") if ctrl else ("PASS" if nbad == 0 else "FAILED")))


if __name__ == "__main__":
    main()
