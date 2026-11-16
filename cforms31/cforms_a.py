#!/usr/bin/env python3
"""Closed forms on pages that had no replay folder (6a, 11b, 35a, 88b, 23b, 46a, 56a, 61a), program A.

Each value is enclosed in an interval of integers at scale 2^-200 (fixed point, directed rounding) with our own square
root (isqrt) and our own natural logarithm (the atanh series ln x = 2 sum_k u^(2k+1)/(2k+1), u = (x-1)/(x+1), after
halving by powers of 2, with an explicit geometric tail). Rationals are exact Fractions. The printed decimal of the
page is then classified: the correct rounding of the value, its truncation, its rounding up, or none of these.
Standard library only.
"""
import sys
from fractions import Fraction as Fr
from math import isqrt

P = 200
ONE = 1 << P


class I:
    """closed interval [lo, hi] of integers at scale 2^-P"""
    def __init__(self, lo, hi):
        assert lo <= hi
        self.lo, self.hi = lo, hi

    @staticmethod
    def frac(q):
        q = Fr(q)
        n, d = q.numerator * ONE, q.denominator
        return I(n // d, -((-n) // d))

    def __add__(s, o): return I(s.lo + o.lo, s.hi + o.hi)
    def __sub__(s, o): return I(s.lo - o.hi, s.hi - o.lo)

    def mul(s, o):
        c = [s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi]
        return I(min(c) >> P, -((-max(c)) >> P))

    def div(s, o):
        assert o.lo > 0 or o.hi < 0
        c = []
        for a in (s.lo, s.hi):
            for b in (o.lo, o.hi):
                c.append(Fr(a * ONE, b))
        lo, hi = min(c), max(c)
        return I(lo.numerator // lo.denominator, -((-hi.numerator) // hi.denominator))

    def lo_fr(s): return Fr(s.lo, ONE)
    def hi_fr(s): return Fr(s.hi, ONE)


def isqrt_iv(q):
    """sqrt of a positive rational q, as an interval"""
    q = Fr(q)
    n = q.numerator * ONE * ONE // q.denominator          # floor(q * 2^(2P))
    lo = isqrt(n)
    n2 = -((-q.numerator * ONE * ONE) // q.denominator)  # ceil
    hi = isqrt(n2)
    if hi * hi < n2:
        hi += 1
    return I(lo, hi)


def ln_iv(x):
    """ln of an interval x with x.lo > 0 (monotone): ln(lo) lower, ln(hi) upper"""
    def ln_point(v, up):
        # v > 0 as Fraction; ln v = k ln 2 + ln m, m in [1, 2)
        k = 0
        while v >= 2:
            v /= 2; k += 1
        while v < 1:
            v *= 2; k -= 1
        return k, v
    def ln_m(m):
        # ln m = 2 atanh(u), u = (m-1)/(m+1) in [0, 1/3); terms u^(2j+1)/(2j+1); tail <= u^(2N+1)/((2N+1)(1-u^2))
        u = (m - 1) / (m + 1)
        s = Fr(0); term = u; j = 0
        u2 = u * u
        while True:
            s += term / (2 * j + 1)
            j += 1
            term *= u2
            if term < Fr(1, 1 << (P + 20)):
                break
        tail = term / (2 * j + 1) / (1 - u2)
        return 2 * s, 2 * (s + tail)
    ln2lo, ln2hi = ln_m(Fr(2))   # m = 2 gives u = 1/3: fine for the series too
    out = []
    for v in (x.lo_fr(), x.hi_fr()):
        k, m = ln_point(v, None)
        a, b = ln_m(m)
        if k >= 0:
            out.append((a + k * ln2lo, b + k * ln2hi))
        else:
            out.append((a + k * ln2hi, b + k * ln2lo))
    lo, hi = out[0][0], out[1][1]
    return I(lo.numerator * ONE // lo.denominator, -((-hi.numerator * ONE) // hi.denominator))


def log2_iv(x):
    return ln_iv(x).div(ln_iv(I.frac(2)))


class Exact:
    """an exact rational value (rows that are rational numbers)"""
    def __init__(self, q): self.q = Fr(q)
    def lo_fr(self): return self.q
    def hi_fr(self): return self.q


def classify(iv, printed):
    """printed: decimal string. Returns the relation of the printed decimal to the value."""
    d = Fr(printed)
    digits = len(printed.split(".")[1]) if "." in printed else 0
    ulp = Fr(1, 10 ** digits)
    lo, hi = iv.lo_fr(), iv.hi_fr()
    rel = []
    if lo == hi == d:
        return "exact"
    if d <= lo and hi < d + ulp:
        rel.append("truncation")
    if d - ulp < lo and hi <= d:
        rel.append("rounded up")
    if d - ulp / 2 <= lo and hi <= d + ulp / 2:
        rel.append("correct rounding")
    return ", ".join(rel) if rel else "NONE"


def dstr(iv, k=30):
    q = iv.lo_fr()      # the lower end (equal to the value for exact rows), truncated
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def rows(control=None):
    s5 = isqrt_iv(5)
    s2 = isqrt_iv(2)
    R = []
    # (page, tag, side, expression, interval, printed, accepted relations)
    psi = I.frac(3).__sub__(s5).div(I.frac(2)) if control != "6a-plus" else (I.frac(3) + s5).div(I.frac(2))
    R.append(("6a", "[AHS2022], [CL2022], [P2022], [S2022]", "lower", "(3 - sqrt 5)/2", psi, "0.381966", ("truncation",)))
    p11 = "0.58497" if control == "11b-print" else "0.58496"
    R.append(("11b", "[KP2020]", "upper", "log2(3/2)", log2_iv(I.frac(Fr(3, 2))), p11, ("correct rounding",)))
    R.append(("35a", "[AP24]", "lower", "log2(1 + sqrt 2)", log2_iv(I.frac(1) + s2), "1.271", ("truncation",)))
    R.append(("88b", "[Pol14b]", "upper", "4 - 28/157", Exact(4 - Fr(28, 157)), "3.821656", ("truncation",)))
    R.append(("23b", "[AJM2025]", "lower", "17/92", Exact(Fr(17, 92)), "0.18478", ("correct rounding",)))
    q15 = Fr(2, 13) if control == "46a-15" else Fr(2, 15)
    for tag, expr, val, pr in (("[Bo1991]", "4 - 2/15", 4 - q15, "3.86667"),
                               ("[Wo1995], [MVV1996]", "4 - 2/11", 4 - Fr(2, 11), "3.81818"),
                               ("[TVV1998]", "4 - 2/9", 4 - Fr(2, 9), "3.77778"),
                               ("[TV2000]", "4 - 2/7", 4 - Fr(2, 7), "3.71429"),
                               ("[Tao2003], [BG2011]", "10/3", Fr(10, 3), "3.33333"),
                               ("[Gut2016]", "13/4", Fr(13, 4), "3.25"),
                               ("[Wan2022]", "3 + 3/13", 3 + Fr(3, 13), "3.23077"),
                               ("[WW2022]", "3 + 3/14", 3 + Fr(3, 14), "3.21429"),
                               ("[WW2024]", "22/7", Fr(22, 7), "3.142857")):
        R.append(("46a", tag, "upper", expr, Exact(val), pr, ("correct rounding", "exact")))
    R.append(("56a", "[LRS1999]", "upper", "1/2 - 1/(2^2 + 1)", Exact(Fr(1, 2) - Fr(1, 2 ** 2 + 1)), "0.3", ("exact",)))
    R.append(("56a", "[KS2003]", "upper", "7/64", Exact(Fr(7, 64)), "0.109375", ("exact",)))
    return R


def main():
    a = sys.argv[1:]
    control = a[a.index("--control") + 1] if "--control" in a else None
    print("program A (fixed point at 2^-%d, own sqrt and log)%s" % (P, (" CONTROL " + control) if control else ""))
    # self-tests
    ln2 = ln_iv(I.frac(2))
    assert ln2.hi - ln2.lo < 1 << 12
    t = log2_iv(I.frac(8))
    assert t.lo_fr() <= 3 <= t.hi_fr(), "log2 8"
    s = isqrt_iv(2)
    assert s.mul(s).lo_fr() <= 2 <= s.mul(s).hi_fr()
    e = ln_iv(I.frac(Fr(3, 2))) + ln_iv(I.frac(Fr(2, 3)))
    assert e.lo_fr() <= 0 <= e.hi_fr()
    print("self-tests: ln 2 width, log2 8 = 3, sqrt(2)^2 = 2, ln(3/2) + ln(2/3) = 0: ok")
    bad = 0
    for page, tag, side, expr, iv, printed, ok in rows(control):
        rel = classify(iv, printed)
        good = any(r in rel.split(", ") for r in ok)
        if not good:
            bad += 1
        safe = ""
        if rel not in ("exact", "NONE"):
            below = Fr(printed) <= iv.lo_fr()
            safe = " (printed %s the value: %s for the %s bound)" % ("below" if below else "above",
                    "safe" if (below == (side == "lower")) else "not on the safe side, an approximation", side)
        print("row %s %s: %s = %s...; page %s: %s%s  %s" % (page, tag, expr, dstr(iv), printed, rel, safe,
                                                          "ok" if good else "MISMATCH"))
    # 61a: the two exponents from the spectral bounds lambda_1 >= 1/4 - theta^2
    for tag, lam, theta in (("[LRS1995]", Fr(3, 16), Fr(1, 4)), ("[KRS2003]", Fr(975, 4096), Fr(7, 64))):
        if control == "61a-sign" and tag == "[KRS2003]":
            good = Fr(1, 4) + theta ** 2 == lam
        else:
            good = Fr(1, 4) - theta ** 2 == lam
        if not good:
            bad += 1
        print("row 61a %s: lambda_1 >= %s = 1/4 - (%s)^2  %s" % (tag, lam, theta, "ok" if good else "MISMATCH"))
    if control:
        print("VERDICT " + ("REFUSED" if bad else "NOT REFUSED"))
    else:
        print("VERDICT " + ("PASS" if bad == 0 else "FAILED"))


if __name__ == "__main__":
    main()
