#!/usr/bin/env python3
"""83a, program A: rigorous bounds for the Wirsing constant C_83 = |lambda_2| of the Gauss-Kuzmin-Wirsing operator.

Rows: [Schweiger1980] C_83 <= 2/(3 + sqrt 5) = 0.38197..., [MR1987] 0.30366299 <= C_83 <= 0.30366327 (eq. (5.20)).
A value from its definition, not a replay of any source's computation.

Wirsing's reduction (Acta Arith. 24 (1974); Knuth, TAOCP vol. 2, sec. 4.5.3): with phi = (1 + x) f, the operator
L f(x) = sum_k (x+k)^-2 f(1/(x+k)) becomes T phi(x) = sum_k p_k(x) phi(1/(x+k)), p_k = (1+x)/((x+k)(x+k+1)),
sum_k p_k = 1; Abel summation with P_k = p_1 + ... + p_k = k/(x+k+1) gives (T phi)' = -V(phi') with the POSITIVE
operator
    V psi(x) = sum_{k>=1} [ k/(x+k+1)^2 int_{1/(x+k+1)}^{1/(x+k)} psi + (1+x)/((x+k)^3 (x+k+1)) psi(1/(x+k)) ].
The spectral radius of V on C[0,1] is |lambda_2| (Wirsing; not re-proved here). For psi > 0 on [0, 1] and constants
a, b with a psi <= V psi <= b psi on [0, 1], positivity gives a^n psi <= V^n psi <= b^n psi, hence a <= rho(V) <= b
(Collatz-Wielandt). psi is the polynomial of psi.json (found by c83a_search.py; not trusted).

This program (route A):
  * V psi exactly in closed form: with y = x + k, X = 1 + x and k = (y + 1) - X, each summand is a sum of terms
    (alpha + beta X) y^-a (y+1)^-b; partial fractions in y and y + 1 (constant coefficients); sum over k:
    sum y^-i = zeta(i, X), sum (y+1)^-i = zeta(i, X) - X^-i; the i = 1 parts must cancel to P_1(X)/X (checked).
    So V psi(x) = sum_{i>=2} R_i(X) zeta(i, X) + E(X), R_i linear, E a Laurent polynomial in X (exact rationals).
  * Hurwitz zeta(n, X0) at rational X0 by Euler-Maclaurin (N = 30 terms, M = 25 Bernoulli corrections) with the
    remainder bound 4 zeta(2M+2) (n)_{2M+1} / ((2 pi)^{2M+2} (X0+N)^{n+2M+1}).
  * [0, 1] cut into cells; on each, the Taylor polynomial of V psi - c psi in s = x - m to order J + 1, with Lagrange
    remainders (zeta(i, .) and X^-i are completely monotone: the (J+1)-th coefficient is at most
    C(i+J, J+1) zeta(i+J+1, X0 - r)); a cell passes if b_0 - sum_{j>=1} |b_j| r^j - remainder > 0.
Arithmetic: exact rationals, every interval end rounded outward to a multiple of 2^-300. Standard library only.
Usage: c83a_a.py [--psi psi.json] [--claim-lower Q] [--claim-upper Q] [--perturb J:Q] [--variant wrong-weight]
                 [--cells N] [--order J] [--selftest-only]
"""
import sys, json, os
from fractions import Fraction as Q
from math import comb, floor

PREC = 300
def down(q):
    return Q(floor(q * (1 << PREC)), 1 << PREC)
def up(q):
    return -down(-q)

class I:
    """closed interval [lo, hi] of rationals, ends rounded outward"""
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None, exact=False):
        if hi is None:
            hi = lo
        if exact:
            self.lo, self.hi = Q(lo), Q(hi)
        else:
            self.lo, self.hi = down(Q(lo)), up(Q(hi))
    def __add__(self, o):
        o = o if isinstance(o, I) else I(o, exact=True)
        return I(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__
    def __neg__(self):
        return I(-self.hi, -self.lo, exact=True)
    def __sub__(self, o):
        return self + (-(o if isinstance(o, I) else I(o, exact=True)))
    def __mul__(self, o):
        if not isinstance(o, I):
            o = Q(o)
            return I(min(self.lo * o, self.hi * o), max(self.lo * o, self.hi * o))
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return I(min(p), max(p))
    __rmul__ = __mul__
    def mag(self):
        return max(abs(self.lo), abs(self.hi))

def bernoulli(n):
    B = [Q(0)] * (n + 1); B[0] = Q(1)
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, j) * B[j] for j in range(m)) / (m + 1)
    return B
EM_N, EM_M = 30, 25
BER = bernoulli(2 * EM_M + 2)
FACT = [1]
for _i in range(1, 2 * EM_M + 4):
    FACT.append(FACT[-1] * _i)
TWO_PI_LOW = Q(628, 100)              # 6.28 < 2 pi

def poch(s, n):
    r = 1
    for i in range(n):
        r *= s + i
    return r

def hurwitz(n, a):
    """zeta(n, a) for integer n >= 2 and rational a >= 1, as an interval (Euler-Maclaurin)"""
    a = Q(a)
    tot = sum(Q(1) / (a + k) ** n for k in range(EM_N))
    aN = a + EM_N
    tot += Q(1) / ((n - 1) * aN ** (n - 1)) + Q(1, 2) / aN ** n
    for j in range(1, EM_M + 1):
        tot += BER[2 * j] / FACT[2 * j] * poch(n, 2 * j - 1) / aN ** (n + 2 * j - 1)
    M1 = EM_M + 1
    zeta2m = 1 + Q(1, 2 * M1 - 1)                    # zeta(2M+2) <= 1 + 1/(2M+1)
    rem = 4 * zeta2m * poch(n, 2 * M1 - 1) / (TWO_PI_LOW ** (2 * M1) * aN ** (n + 2 * M1 - 1))
    return I(tot - rem, tot + rem)

def pf(a, b):
    """1/(y^a (y+1)^b) = sum A_i y^-i + sum B_i (y+1)^-i"""
    if b == 0:
        return [(("y", a), 1)]
    if a == 0:
        return [(("z", b), 1)]
    out = [(("y", i), (-1) ** (a - i) * comb(a + b - i - 1, b - 1)) for i in range(1, a + 1)]
    out += [(("z", i), (-1) ** a * comb(a + b - i - 1, a - 1)) for i in range(1, b + 1)]
    return out

def closed_form(c, variant=None):
    """V psi = sum_{i>=2} (al_i + be_i X) zeta(i, X) + sum_{i>=0} e_i X^-i, psi(t) = sum c_j t^j"""
    acc = {}
    def add(a, b, al, be):
        for key, co in pf(a, b):
            v = acc.setdefault(key, [Q(0), Q(0)])
            v[0] += al * co; v[1] += be * co
    for j, cj in enumerate(c):
        if cj == 0:
            continue
        d = cj / (j + 1)
        # k/(y+1)^2 = 1/(y+1) - X/(y+1)^2, times int = d (y^-(j+1) - (y+1)^-(j+1))
        add(j + 1, 1, d, 0); add(0, j + 2, -d, 0); add(j + 1, 2, 0, -d); add(0, j + 3, 0, d)
        # (1+x)/(y^3 (y+1)) psi(1/y)
        add(j + 3 if variant != "wrong-weight" else j + 2, 1, 0, cj)
    y1 = acc.pop(("y", 1), [Q(0), Q(0)]); z1 = acc.pop(("z", 1), [Q(0), Q(0)])
    if y1[0] + z1[0] != 0 or y1[1] + z1[1] != 0:
        raise SystemExit("the 1/y terms do not cancel")
    R = {}; E = {}
    def addE(p, v):
        E[p] = E.get(p, Q(0)) + v
    addE(1, y1[0]); addE(0, y1[1])                    # (al + be X)/X
    for (kind, i), (al, be) in acc.items():
        r = R.setdefault(i, [Q(0), Q(0)])
        r[0] += al; r[1] += be
        if kind == "z":                               # zeta(i, X + 1) = zeta(i, X) - X^-i
            addE(i, -al); addE(i - 1, -be)
    return R, E

def dec(q):
    """exact decimal string of a rational with a power-of-ten denominator, else a 20-digit approximation"""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 60
        return str(+(Decimal(q.numerator) / Decimal(q.denominator)))

def selftests(c):
    ok = True
    # partial fractions at a few rationals
    for (a, b) in [(1, 1), (3, 2), (5, 7), (20, 3), (0, 4), (6, 0)]:
        for y in (Q(7, 3), Q(11, 2)):
            lhs = Q(1) / (y ** a * (y + 1) ** b)
            rhs = sum(co * (y ** -i if k == "y" else (y + 1) ** -i) for (k, i), co in pf(a, b))
            ok &= lhs == rhs
    print("self-test partial fractions:", ok)
    # zeta against known constants (decimal strings of pi^2/6 and Apery's constant, 45 digits)
    z2 = hurwitz(2, 1); z3 = hurwitz(3, 1)
    k2 = Q("1.644934066848226436472415166646025189218949901")
    k3 = Q("1.202056903159594285399738161511449990764986292")
    t1 = z2.lo - Q(1, 10 ** 44) <= k2 <= z2.hi + Q(1, 10 ** 44) and z2.hi - z2.lo < Q(1, 10 ** 50)
    t2 = z3.lo - Q(1, 10 ** 44) <= k3 <= z3.hi + Q(1, 10 ** 44) and z3.hi - z3.lo < Q(1, 10 ** 50)
    # zeta(n, a) - zeta(n, a + 1) = a^-n
    d = hurwitz(5, Q(7, 4)) - hurwitz(5, Q(11, 4)) - Q(4, 7) ** 5
    t3 = d.lo <= 0 <= d.hi and d.hi - d.lo < Q(1, 10 ** 50)
    print("self-test zeta(2) = pi^2/6:", t1, "; zeta(3) = Apery:", t2, "; zeta(5, a) - zeta(5, a+1) = a^-5:", t3)
    ok &= t1 and t2 and t3
    # the closed form against direct summation (floating point, k <= 200000 plus the integral tail), x = 0.37
    from math import fsum
    x = 0.37
    def psi(t):
        r = 0.0
        for v in reversed(c):
            r = r * t + float(v)
        return r
    def Psi(t):
        r = 0.0
        for j in range(len(c) - 1, -1, -1):
            r = r * t + float(c[j]) / (j + 1)
        return r * t
    K = 200000
    s = fsum(k / (x + k + 1) ** 2 * (Psi(1 / (x + k)) - Psi(1 / (x + k + 1)))
             + (1 + x) / ((x + k) ** 3 * (x + k + 1)) * psi(1 / (x + k)) for k in range(1, K + 1))
    s += psi(0) / (2 * (K + 0.5) ** 2)               # tail ~ psi(0) sum k^-3
    R, E = closed_form(c)
    X = Q(137, 100)
    cf = sum((al + be * X) * hurwitz(i, X) for i, (al, be) in R.items()) + sum(v * X ** -p for p, v in E.items())
    t4 = abs(float(cf.lo) - s) < 1e-9
    print("self-test closed form = direct sum at x = 0.37: %s (%.15f vs %.15f)" % (t4, float(cf.lo), s))
    ok &= t4
    return ok

def taylor_shift(c, m):
    """coefficients of p(m + s) in s"""
    n = len(c)
    out = [Q(0)] * n
    for j in range(n - 1, -1, -1):          # Horner: out = out*(m+s) + c_j
        new = [Q(0)] * n
        for i in range(n - 1):
            new[i + 1] += out[i]
        for i in range(n):
            new[i] += out[i] * m
        new[0] += c[j]
        out = new
    return out

def cell(R, E, c, m, r, J):
    """Taylor coefficients (intervals) of V psi at x = m + s, |s| <= r, and the remainder bound"""
    X0 = 1 + m
    deg = J + 1
    b = [I(0, exact=True) for _ in range(deg + 1)]
    rem = Q(0)
    imax = max(R)
    Z = {n: hurwitz(n, X0) for n in range(2, imax + J + 2)}
    for i, (al, be) in R.items():
        lin0 = al + be * X0
        zs = [Z[i + j] * ((-1) ** j * comb(i + j - 1, j)) for j in range(J + 1)]
        for j in range(J + 1):
            b[j] = b[j] + zs[j] * lin0
            b[j + 1] = b[j + 1] + zs[j] * be
        n = i + J + 1
        zb = 1 + Q(1, n - 1)                       # zeta(n, xi) <= zeta(n, 1) <= 1 + 1/(n-1), xi >= 1
        rem += (abs(lin0) + abs(be) * r) * comb(i + J, J + 1) * zb * r ** (J + 1)
    for p, v in E.items():
        if p == 0:
            b[0] = b[0] + v
            continue
        for j in range(J + 1):
            b[j] = b[j] + v * (-1) ** j * comb(p + j - 1, j) / X0 ** (p + j)
        rem += abs(v) * comb(p + J, J + 1) * r ** (J + 1)    # (X0 - r)^-(p+J+1) <= 1
    return b, rem

def lower_bound(b, rem, r):
    return b[0].lo - sum(b[j].mag() * r ** j for j in range(1, len(b))) - rem

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    a = sys.argv[1:]
    opt = {"--psi": os.path.join(here, "psi.json"), "--cells": "32", "--order": "24"}
    flags = set()
    i = 0
    while i < len(a):
        if a[i] in ("--selftest-only",):
            flags.add(a[i]); i += 1
        else:
            opt[a[i]] = a[i + 1]; i += 2
    data = json.load(open(opt["--psi"]))
    c = [Q(v) for v in data["psi_monomial_coefficients"]]
    lo = Q(opt.get("--claim-lower", data["claimed_lower"]))
    hi = Q(opt.get("--claim-upper", data["claimed_upper"]))
    if "--perturb" in opt:
        j, d = opt["--perturb"].split(":")
        c[int(j)] += Q(d)
        print("CONTROL: coefficient %s of psi changed by %s" % (j, d))
    variant = opt.get("--variant")
    if variant:
        print("CONTROL: operator variant", variant, "(the point term with (x+k)^2 in place of (x+k)^3)")
    print("83a program A: psi of degree %d, claimed %s <= |lambda_2| <= %s" % (len(c) - 1, dec(lo), dec(hi)))
    if not selftests(c):
        print("VERDICT FAILED (self-test)"); return
    if "--selftest-only" in flags:
        return
    R, E = closed_form(c, variant)
    print("closed form: zeta(i, X) terms for i = %d..%d, Laurent terms X^-p for p = %d..%d; max |R_i| coefficient %.3e"
          % (min(R), max(R), min(E), max(E), float(max(max(abs(u), abs(v)) for u, v in R.values()))))
    NC, J = int(opt["--cells"]), int(opt["--order"])
    r = Q(1, 2 * NC)
    worst = {"psi": None, "lo": None, "hi": None}
    fail = {"psi": [], "lo": [], "hi": []}
    for k in range(NC):
        m = Q(2 * k + 1, 2 * NC)
        vb, rem = cell(R, E, c, m, r, J)
        pb = taylor_shift(c, m) + [Q(0)] * (J + 2 - len(c))
        pb = pb[:J + 2] if len(pb) >= J + 2 else pb
        psib = [I(v, exact=True) for v in pb]
        L_psi = lower_bound(psib, Q(0), r) if len(c) <= J + 2 else None
        g_lo = [vb[j] - psib[j] * lo for j in range(J + 2)]
        g_hi = [psib[j] * hi - vb[j] for j in range(J + 2)]
        L_lo = lower_bound(g_lo, rem, r); L_hi = lower_bound(g_hi, rem, r)
        for key, L, g in (("psi", L_psi, psib), ("lo", L_lo, g_lo), ("hi", L_hi, g_hi)):
            if worst[key] is None or L < worst[key][0]:
                worst[key] = (L, m)
            if L <= 0:
                fail[key].append((m, g[0]))
    print("cells: %d of width %s, Taylor order %d" % (NC, 2 * r, J + 1))
    print("psi > 0 on [0,1]: %s (least cell bound %.6e at x = %s)" % (not fail["psi"], float(worst["psi"][0]), worst["psi"][1]))
    print("V psi - %s psi >= 0 on [0,1]: %s (least cell bound %.6e at x = %s)"
          % (dec(lo), not fail["lo"], float(worst["lo"][0]), worst["lo"][1]))
    print("%s psi - V psi >= 0 on [0,1]: %s (least cell bound %.6e at x = %s)"
          % (dec(hi), not fail["hi"], float(worst["hi"][0]), worst["hi"][1]))
    for key in ("lo", "hi"):
        for m, v in fail[key][:3]:
            if v.hi < 0:
                print("  the inequality FAILS at x = %s: value in [%.6e, %.6e]" % (m, float(v.lo), float(v.hi)))
            else:
                print("  not proved on the cell at x = %s (value at the centre in [%.6e, %.6e])" % (m, float(v.lo), float(v.hi)))
    if fail["psi"] or fail["lo"] or fail["hi"]:
        print("VERDICT REFUSED: the claimed bounds are not proved")
        return
    print("PROVED: %s <= |lambda_2| <= %s" % (dec(lo), dec(hi)))
    # the rows: 2/(3 + sqrt 5) > hi  <=>  sqrt 5 < 2/hi - 3  (both sides positive)
    s5 = 2 / hi - 3
    row_s = s5 > 0 and s5 * s5 > 5
    row_mr_up = hi <= Q(30366327, 10 ** 8)
    row_mr_lo = lo >= Q(30366299, 10 ** 8)
    print("row [Schweiger1980] C_83 <= 2/(3 + sqrt 5) = 0.381966...:", row_s)
    print("row [MR1987] C_83 <= 0.30366327:", row_mr_up)
    print("row [MR1987] C_83 >= 0.30366299:", row_mr_lo)
    centre = Q("0.30366300289873265859744812190155623311087735225365")
    print("information: the 175-digit centre of [GKW2026] (first 50 digits) lies in [lo, hi]:", lo <= centre <= hi)
    print("VERDICT PASS" if (row_s and row_mr_up and row_mr_lo) else "VERDICT REFUSED: a row does not follow")

if __name__ == "__main__":
    main()
