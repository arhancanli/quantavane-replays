"""24a: the row C24 <= 6.9013* (Akbas and Sra, arXiv 2609.27172v2, Theorem 2.4 (label thm:cosine), added in PR #250),
replayed with our own code (Quantavane replays).  Standard library only; run with python3 -I.

The theorem's constant C_cos is the unique solution of equation (eq:cosine-constant),
    F(C) := (3 sqrt(2 pi) / C) * int_0^3 (1 - r/3)^2 exp(-pi^2 r^2 / (2 C^2)) dr = 1,
and the row is C_cos < 6.9013 (the paper: C_cos in (6.9012574, 6.9012575)).  The equation is read from the paper's
LaTeX source (the arXiv e-print of v2, SHA-256 checked by fetch.sh) and compared with our transcription first.

F is strictly decreasing in C (so the solution is unique and C_cos < c exactly when F(c) < 1): with s = pi/C and
u = s r, dF/ds is a positive multiple of int_0^{3s} w(u) e^{-u^2/2} (1 - u^2) du with w(u) = (1 - u/(3s))^2
non-negative and decreasing; int_0^T e^{-u^2/2}(1 - u^2) du = T e^{-T^2/2} > 0, and (w(u) - w(1))(1 - u^2) >= 0, so
the integral is at least w(1) T e^{-T^2/2} > 0 (and trivially positive when T <= 1).  This is our own argument; the
paper states the monotonicity.

Two evaluations of F, both in fixed-point interval arithmetic (2^-400 units, outward rounding), with pi enclosed by
Gauss's formula pi = 48 atan(1/18) + 32 atan(1/57) - 20 atan(1/239) in exact rationals:
  route A: term by term, int_0^3 (1-r/3)^2 r^(2j) dr = 2 * 3^(2j+1) / ((2j+1)(2j+2)(2j+3)), an alternating series in
           z = 9 pi^2 / (2 C^2) with the remainder bounded by the first omitted term;
  route B: through the error-function form, I = J0 (1 + 1/(18a)) - (1 - E)/(3a) - E/(6a) with a = pi^2/(2C^2),
           E = e^(-9a) (Taylor series, geometric tail) and J0 = int_0^3 e^(-a r^2) dr = 3 E sum_n (18a)^n/(2n+1)!!
           (positive terms, geometric tail).
Both must decide the sign of F(c) - 1 at every c used, and agree.  C_cos is then enclosed by bisection.

What rests on mathematics (read, not computed): the proof of Theorem 2.4 (the cosine signing and its Fourier
representation) and Theorem 1.1's reduction; this replay checks the constant, not the discrepancy argument.

usage: python3 -I komlos_check.py PAPER_DIR [--controls]
"""
import hashlib
import math
import os
import re
import sys
from fractions import Fraction as Q

TRANSCRIPTION = (r"\frac{3\sqrt{2\pi}}{C_{\mathrm{cos}}}\int_0^3(1-r/3)^2"
                 r"\exp\!\left(-\frac{\pi^2r^2}{2C_{\mathrm{cos}}^2}\right)\,dr=1.")
SH = 400


def read_equation(paper_dir):
    data = open(os.path.join(paper_dir, "komlos-paper.tex"), "rb").read()
    text = data.decode("utf-8", "replace")
    m = re.search(r"\\label\{eq:cosine-constant\}\s*(.*?)\\end\{equation\}", text, re.S)
    if not m:
        raise SystemExit("FAIL: equation (eq:cosine-constant) not found")
    return hashlib.sha256(data).hexdigest(), re.sub(r"\s+", "", m.group(1))


# ------------------------------------------------------------------ fixed-point intervals
def fl(q):
    q = Q(q)
    return (q.numerator << SH) // q.denominator


def ce(q):
    q = Q(q)
    return -((-q.numerator << SH) // q.denominator)


def iv(q):
    return (fl(q), ce(q))


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def sub(a, b):
    return (a[0] - b[1], a[1] - b[0])


def mul(a, b):
    p = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return (min(p) >> SH, -((-max(p)) >> SH))


def div(a, b):
    assert b[0] > 0, "division by a non-positive interval"
    qs = [Q(x, y) for x in a for y in b]
    return (fl(min(qs)), ce(max(qs)))


def sqrt_iv(a):
    assert a[0] >= 0
    return (math.isqrt(a[0] << SH), math.isqrt(a[1] << SH) + 1)


def hull(a, b):
    return (min(a[0], b[0]), max(a[1], b[1]))


def to_q(x):
    return Q(x, 1 << SH)


def atan_inv(x, terms=120):
    """atan(1/x) for an integer x >= 2, as a rational enclosure (alternating series, remainder < next term)"""
    s = Q(0)
    for k in range(terms):
        s += Q((-1) ** k, (2 * k + 1) * x ** (2 * k + 1))
    nxt = Q(1, (2 * terms + 1) * x ** (2 * terms + 1))
    return (s - nxt, s + nxt)


def pi_iv():
    a, b, c = atan_inv(18), atan_inv(57), atan_inv(239)
    lo = 48 * a[0] + 32 * b[0] - 20 * c[1]
    hi = 48 * a[1] + 32 * b[1] - 20 * c[0]
    assert Q(314159, 100000) < lo < hi < Q(314160, 100000)
    return (fl(lo), ce(hi))


PI = pi_iv()
SQRT_2PI = sqrt_iv(mul(iv(2), PI))


# ------------------------------------------------------------------ route A
def F_route_a(C, p=2, N=90, half=Q(1, 2)):
    """F(C) by the term-by-term series; p is the power of the weight (1 - r/3)^p (2 in the paper; a control uses 1);
    half is the factor in exp(-half * pi^2 r^2 / C^2) (1/2 in the paper)."""
    Ci = iv(C)
    a = div(mul(iv(half), mul(PI, PI)), mul(Ci, Ci))       # a = half pi^2 / C^2
    # int_0^3 (1-r/3)^p r^(2j) dr = 3^(2j+1) B(2j+1, p+1) = 3^(2j+1) (2j)! p! / (2j+p+1)!
    total = iv(0)
    apow = iv(1)
    last = None
    for j in range(N + 1):
        coef = Q(3 ** (2 * j + 1) * math.factorial(2 * j) * math.factorial(p), math.factorial(2 * j + p + 1) * math.factorial(j))
        term = mul(apow, iv(coef))
        if j == N:
            last = term
            break
        total = add(total, term) if j % 2 == 0 else sub(total, term)
        apow = mul(apow, a)
    # |t_(j+1)| / |t_j| = 9a (2j+1)(2j+2) / ((j+1)(2j+p+2)(2j+p+3)) <= 9a/(j+1) < 1 for j >= N once 9a < N + 1, so
    # the series alternates with decreasing terms from N on and the remainder is bounded by |t_N|
    assert 9 * Q(a[1], 1 << SH) < N + 1
    rem = (-last[1], last[1])
    integral = add(total, rem)
    return div(mul(iv(3), mul(SQRT_2PI, integral)), Ci)


# ------------------------------------------------------------------ route B
def exp_neg(y, N=120):
    """e^(-y) for an interval y > 0: 1/e^y, e^y by Taylor with the tail <= t_N / (1 - y/(N+1))"""
    s, t = iv(1), iv(1)
    for n in range(1, N):
        t = div(mul(t, y), iv(n))
        s = add(s, t)
    tN = div(mul(t, y), iv(N))
    q = Q(y[1], 1 << SH) / (N + 1)
    assert q < 1
    tail = (0, ce(to_q(tN[1]) / (1 - q)))
    ey = add(s, tail)
    return div(iv(1), ey)


def F_route_b(C, N=150):
    Ci = iv(C)
    a = div(mul(PI, PI), mul(iv(2), mul(Ci, Ci)))           # a = pi^2 / (2 C^2)
    E = exp_neg(mul(iv(9), a))
    x = mul(iv(18), a)
    s, t = iv(1), iv(1)
    for n in range(1, N):
        t = div(mul(t, x), iv(2 * n + 1))
        s = add(s, t)
    q = Q(x[1], 1 << SH) / (2 * N + 1)
    assert q < 1
    tN = div(mul(t, x), iv(2 * N + 1))
    s = add(s, (0, ce(to_q(tN[1]) / (1 - q))))
    J0 = mul(iv(3), mul(E, s))
    one = iv(1)
    I = sub(sub(mul(J0, add(one, div(one, mul(iv(18), a)))), div(sub(one, E), mul(iv(3), a))), div(E, mul(iv(6), a)))
    return div(mul(iv(3), mul(SQRT_2PI, I)), Ci)


def sign_minus_one(F):
    if F[0] > (1 << SH):
        return 1
    if F[1] < (1 << SH):
        return -1
    return 0


def fmt(x, digits=40):
    q = to_q(x)
    s = "-" if q < 0 else ""
    q = abs(q)
    ip = int(q)
    frac = int((q - ip) * 10 ** digits)
    return f"{s}{ip}.{frac:0{digits}d}"


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        raise SystemExit(__doc__)
    h, tex = read_equation(args[0])
    same = tex == re.sub(r"\s+", "", TRANSCRIPTION)
    print(f"komlos-paper.tex sha256 {h}")
    print(f"(eq:cosine-constant) in the source: {tex}")
    print(f"our transcription equals it: {same}")
    print(f"pi in [{fmt(PI[0], 50)}, {fmt(PI[1], 50)}] (Gauss's arctangent formula)")
    ok = same
    points = ["6", "6.9012574", "6.9012575", "6.9013"]
    signs = {}
    for c in points:
        fa, fb = F_route_a(Q(c)), F_route_b(Q(c))
        sa, sb = sign_minus_one(fa), sign_minus_one(fb)
        agree = sa == sb and sa != 0
        ok &= agree
        signs[c] = sa
        print(f"F({c}) - 1: route A in [{fmt(fa[0] - (1 << SH), 30)}, {fmt(fa[1] - (1 << SH), 30)}], "
              f"route B in [{fmt(fb[0] - (1 << SH), 30)}, {fmt(fb[1] - (1 << SH), 30)}]; sign {sa:+d}, routes agree: {agree}")
    enclosure = signs["6.9012574"] == 1 and signs["6.9012575"] == -1
    bound = signs["6.9013"] == -1
    print(f"the paper's enclosure C_cos in (6.9012574, 6.9012575): {enclosure}")
    print(f"the row's bound C_cos < 6.9013 (F decreasing, F(6.9013) < 1): {bound}")
    ok &= enclosure and bound and signs["6"] == 1
    # bisection with route A, checked by route B at the ends
    lo, hi = Q("6.9012574"), Q("6.9012575")
    for _ in range(100):
        mid = (lo + hi) / 2
        s = sign_minus_one(F_route_a(mid))
        if s == 0:
            break
        if s > 0:
            lo = mid
        else:
            hi = mid
    ends = sign_minus_one(F_route_b(lo)) == 1 and sign_minus_one(F_route_b(hi)) == -1
    ok &= ends
    print(f"C_cos in [{float(lo):.17g}, {float(hi):.17g}] (width {float(hi - lo):.1e}): "
          f"{str(lo.numerator * 10 ** 30 // lo.denominator)[:1]}.{str(lo.numerator * 10 ** 30 // lo.denominator)[1:]}...; "
          f"route B agrees at both ends: {ends}")
    if "--controls" in sys.argv:
        print("\nred controls (each must be refused or move the root past the printed bound):")
        c1 = sign_minus_one(F_route_a(Q("6.9012"))) == 1 and sign_minus_one(F_route_b(Q("6.9012"))) == 1
        print(f"  the bound tightened to 6.9012: F(6.9012) > 1 in both routes, so C_cos < 6.9012 is refused: {c1}")
        f2 = F_route_a(Q("6.9013"), p=1)
        c2 = sign_minus_one(f2) == 1
        print(f"  weight (1 - r/3) instead of (1 - r/3)^2: F(6.9013) - 1 in [{fmt(f2[0] - (1 << SH), 12)}, ...] > 0, "
              f"so 6.9013 would not bound that root: {c2}")
        f3a, f3b = F_route_a(Q("6.9012574"), half=Q(1)), F_route_a(Q("6.9012575"), half=Q(1))
        c3 = not (sign_minus_one(f3a) == 1 and sign_minus_one(f3b) == -1)
        print(f"  the 2 dropped from 2 C^2 in the exponent: the enclosure (6.9012574, 6.9012575) is refused: {c3}")
        f4 = F_route_a(Q("6.9013"), N=3)
        c4 = sign_minus_one(f4) == 0
        print(f"  route A cut at three terms with its remainder bound: undecided at 6.9013 (width "
              f"{float(to_q(f4[1] - f4[0])):.1e}), not a pass: {c4}")
        cok = c1 and c2 and c3 and c4
        print(f"red controls all as expected: {cok}")
        ok &= cok
    print(f"RESULT: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
