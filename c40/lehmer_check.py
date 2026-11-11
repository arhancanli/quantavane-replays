#!/usr/bin/env python3
"""Exact replay of the upper-bound row of the 40a page (teorth/optimizationproblems, main 87b4ebb):
[Leh1933]/[BDM2007] L <= M(l) = 1.176280..., l(x) = x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1.

M(f) is |leading coefficient| times the product of |roots| outside the unit disc (the page's definition).  Written
from the page alone; the polynomial is read from the page's text as data.  Exact rational arithmetic only:

route 1  l is palindromic, so l(x) = x^5 T(x + 1/x) with T of degree 5 (computed from x^k + x^-k = V_k(y)).  Sturm
         sequences of T (squarefree, checked) count its real roots in (-inf, -2), (-2, 2), (2, inf); all five real,
         with T(+-2) != 0, means every root of l is on the unit circle except lambda > 1 and 1/lambda, where
         lambda + 1/lambda = y0 is the root of T in (2, inf).  So M(l) = lambda.  y0 is enclosed by bisection and
         lambda = (y0 + sqrt(y0^2 - 4))/2 by integer square roots, rounded outward.
route 2  Sturm sequences of l itself: exactly one real root in (1, inf) and exactly one in (0, 1), none <= 0; that
         root is enclosed by bisection directly.
The two enclosures must overlap; the printed 1.176280... must be its leading digits.
usage: python3 -I lehmer_check.py CERTS_ROOT [--controls]
"""
import re
import sys
from fractions import Fraction
from math import isqrt

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv
BITS = 200          # bisection to width 2^-200


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def parse_page(path):
    """Coefficients (index = power) of the polynomial l(x) printed in the 40a page's upper-bound table."""
    text = open(path, encoding="utf-8").read()
    m = re.search(r"\\ell\(x\)=([^$]*)\$", text)
    need(m is not None, "Lehmer's polynomial was not found on the page")
    coeffs = {}
    for term in re.findall(r"[+-]?[^+-]+", m.group(1).replace("{", "").replace("}", "").replace(" ", "")):
        sign = -1 if term[0] == "-" else 1
        body = term.lstrip("+-")
        if "x" in body:
            c, _, rest = body.partition("x")
            p = int(rest[1:]) if rest.startswith("^") else 1
            need(rest == "" or re.fullmatch(r"\^\d+", rest) is not None, f"unexpected term {term!r}")
        else:
            c, p = body, 0
        coeffs[p] = coeffs.get(p, 0) + sign * (int(c) if c else 1)
    deg = max(coeffs)
    return [coeffs.get(i, 0) for i in range(deg + 1)], m.group(1)


def pstrip(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def pev(p, x):
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def pderiv(p):
    return pstrip([i * p[i] for i in range(1, len(p))] or [0])


def pdivmod(a, b):
    a = [Fraction(c) for c in a]
    b = [Fraction(c) for c in b]
    q = [Fraction(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        k = len(a) - len(b)
        f = a[-1] / b[-1]
        q[k] = f
        for i in range(len(b)):
            a[i + k] -= f * b[i]
        a = pstrip(a)
        if len(a) < len(b) or (len(a) == 1 and a[0] == 0):
            break
    return pstrip(q), pstrip(a)


def sturm(p):
    seq = [pstrip([Fraction(c) for c in p]), pderiv([Fraction(c) for c in p])]
    while len(seq[-1]) > 1 or seq[-1][0] != 0:
        _, r = pdivmod(seq[-2], seq[-1])
        if len(r) == 1 and r[0] == 0:
            break
        seq.append([-c for c in r])
    return seq


def sign_changes(seq, x):
    """Sign changes of the Sturm sequence at x (x = 'inf' or '-inf' allowed)."""
    vals = []
    for p in seq:
        if x == "inf":
            v = p[-1]
        elif x == "-inf":
            v = p[-1] * (-1) ** (len(p) - 1)
        else:
            v = pev(p, x)
        if v != 0:
            vals.append(v > 0)
    return sum(1 for a, b in zip(vals, vals[1:]) if a != b)


def count(seq, a, b):
    """Distinct real roots in (a, b], for a squarefree p with p(b) != 0 the open interval."""
    return sign_changes(seq, a) - sign_changes(seq, b)


def squarefree(p):
    seq = sturm(p)
    return len(seq[-1]) == 1     # the last nonzero remainder is a constant: gcd(p, p') = 1


def bisect(p, lo, hi, bits=BITS):
    flo = pev(p, lo)
    need(flo != 0 and pev(p, hi) != 0 and (flo > 0) != (pev(p, hi) > 0), "no sign change on the bracket")
    width = Fraction(1, 1 << bits)
    while hi - lo > width:
        mid = (lo + hi) / 2
        fm = pev(p, mid)
        if fm == 0:
            return mid, mid
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return lo, hi


def to_T(l):
    """T with l(x) = x^m T(x + 1/x), for a palindromic l of even degree 2m."""
    n = len(l) - 1
    need(n % 2 == 0 and all(l[i] == l[n - i] for i in range(n + 1)), "not a palindromic polynomial of even degree")
    m = n // 2
    V = [[2], [0, 1]]                       # V_k(y) = x^k + x^-k
    for k in range(2, m + 1):
        a = [0] + V[k - 1]
        b = V[k - 2] + [0] * (len(a) - len(V[k - 2]))
        V.append([a[i] - b[i] for i in range(len(a))])
    T = [0] * (m + 1)
    T[0] += l[m]
    for k in range(1, m + 1):
        for i, c in enumerate(V[k]):
            T[i] += l[m + k] * c
    return T


def sqrt_bounds(x, scale):
    """floor and ceil of sqrt(x) * scale for a non-negative rational x."""
    num = x.numerator * scale * scale
    den = x.denominator
    q, r = divmod(num, den)
    lo = isqrt(q)
    hi = lo if lo * lo * den == num else lo + 1
    return lo, hi


def measure(l, label, verbose=True):
    """Return (lower, upper) rational enclosure of M(l) by route 1, after checking route 1's hypotheses, or None."""
    T = to_T(l)
    need(squarefree(T), f"{label}: T is not squarefree")
    seq = sturm(T)
    lo, hi = Fraction(-2), Fraction(2)
    need(pev(T, lo) != 0 and pev(T, hi) != 0, f"{label}: T(+-2) = 0 (a root +-1 of l)")
    big = sum(abs(c) for c in T)            # every real root lies in (-big, big)
    below, mid, above = count(seq, -big, lo), count(seq, lo, hi), count(seq, hi, big)
    if verbose:
        print(f"  {label}: T = {T} (constant term first); real roots of T: {below} in (-inf, -2), {mid} in (-2, 2), "
              f"{above} in (2, inf); degree {len(T) - 1}")
    if below + mid + above != len(T) - 1:
        if verbose:
            print(f"  {label}: T has non-real roots, so this route does not apply; refused")
        return None
    need(below == 0, f"{label}: a root of T below -2 (a negative real root of l)")
    if above == 0:
        return Fraction(1), Fraction(1)
    need(above == 1, f"{label}: more than one root of T above 2")
    a, b = bisect(T, Fraction(2), Fraction(big))
    scale = 10 ** 60
    # lambda = (y + sqrt(y^2 - 4)) / 2 is increasing in y > 2
    lo_s = sqrt_bounds(a * a - 4, scale)[0]
    hi_s = sqrt_bounds(b * b - 4, scale)[1]
    lam_lo = (a + Fraction(lo_s, scale)) / 2
    lam_hi = (b + Fraction(hi_s, scale)) / 2
    return lam_lo * abs(l[-1]), lam_hi * abs(l[-1])


def direct(l, label):
    """Route 2: Sturm on l itself, and its root in (1, 2] by bisection."""
    need(squarefree(l), f"{label}: l is not squarefree")
    seq = sturm(l)
    big = sum(abs(c) for c in l)
    neg = count(seq, -big, Fraction(0))
    in01 = count(seq, Fraction(0), Fraction(1)) - (1 if pev(l, Fraction(1)) == 0 else 0)
    above = count(seq, Fraction(1), big)
    print(f"  {label}: real roots of l: {neg} in (-inf, 0], {in01} in (0, 1), {above} in (1, inf)")
    need(neg == 0 and in01 == 1 and above == 1, f"{label}: not the root pattern of a Salem polynomial")
    a, b = bisect(l, Fraction(1), Fraction(2))
    return a, b


def digits(x, n):
    return f"{x.numerator * 10 ** n // x.denominator}"


def main():
    l, text = parse_page(f"{ROOT}/page/40a.md")
    print(f"page: {ROOT}/page/40a.md; l(x) = {text}; coefficients (constant term first) {l}")
    need(l == [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], "the parsed polynomial is not Lehmer's")
    ok = True
    r1 = measure(l, "route 1")
    r2 = direct(l, "route 2")
    need(r1 is not None, "route 1 did not apply")
    lo, hi = max(r1[0], r2[0]), min(r1[1], r2[1])
    need(lo <= hi, "the two routes do not overlap")
    print(f"  route 1: M(l) in [{float(r1[0])!r}, ...], width {float(r1[1] - r1[0]):.1e}; route 2: width "
          f"{float(r2[1] - r2[0]):.1e}; they overlap")
    d_lo, d_hi = digits(lo, 30), digits(hi, 30)
    common = 0
    while common < len(d_lo) and d_lo[common] == d_hi[common]:
        common += 1
    print(f"  M(l) = {d_lo[0]}.{d_lo[1:common]}... (the first {common - 1} decimals certain)")
    printed = Fraction(1176280, 10 ** 6)
    good = printed <= lo < printed + Fraction(1, 10 ** 6)
    ok &= good
    print(f"RESULT [Leh1933]/[BDM2007]: M(l) = {d_lo[0]}.{d_lo[1:21]}..., so L <= M(l); the printed 1.176280... is "
          f"{'its leading digits: HOLDS' if good else 'NOT its leading digits: FAILS'}")
    # anchor: a cyclotomic polynomial has measure 1 (the page's Kronecker row)
    phi7 = [1] * 7
    print("anchor: Phi_7 = 1 + x + ... + x^6 (palindromic, degree 6):")
    m7 = measure(phi7, "Phi_7")
    print(f"  M(Phi_7) = {m7[0]} (every root on the unit circle), as Kronecker's theorem says")
    need(m7 == (1, 1), "anchor failed")
    if CONTROLS:
        print("red controls (each must be refused or give a different value):")
        res = []
        l2 = list(l)
        l2[5] = -2                       # -2 x^5 instead of -x^5 (still palindromic)
        r = measure(l2, "l with -2x^5")
        if r is None:
            res.append(("-2 x^5 for -x^5", True))
        else:
            ld = digits(r[0], 6)
            print(f"  l with -2x^5: M = {ld[0]}.{ld[1:]}..., not 1.176280")
            res.append(("-2 x^5 for -x^5", not (printed <= r[0] < printed + Fraction(1, 10 ** 6))))
        l3 = list(l)
        l3[4] = 0                        # drop -x^4: no longer palindromic
        try:
            measure(l3, "l without -x^4")
            res.append(("l without -x^4", False))
        except AssertionError as e:
            print(f"  l without -x^4: refused ({e})")
            res.append(("l without -x^4", True))
        wrong = Fraction(1176281, 10 ** 6)
        refused = not (wrong <= lo < wrong + Fraction(1, 10 ** 6))
        print(f"  a printed 1.176281...: {'refused' if refused else 'NOT refused'}")
        res.append(("printed 1.176281", refused))
        for label, r in res:
            print(f"  {label}: {'as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
