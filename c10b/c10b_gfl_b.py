#!/usr/bin/env python3
"""Program B for the 10b lower bound [GFL26] (arXiv 2609.07000v1, K_G^C > 1.35584631827168, recorded in PR #196): the
finite conditions of the paper's Theorem 3.1 for its Section 4 parameters, by a second route that shares no code with
c10b_gfl_a.py. Interval arithmetic on Python's decimal module at 100 digits, every operation rounded outward.

Route (ours, different from program A's):
  * radial factors in the Laguerre form the paper's Appendix states (a >= b):
        e_{a,b}(r) = (-1)^b sqrt(b!/a!) r^(a-b) L_b^(a-b)(r^2),  L_j^(al)(s) = sum_k (-1)^k C(j+al, j-k) s^k / k!,
    and the mirror form for b > a; the square-root factors of a product are combined into one interval square root;
  * with Delta = (s + d1)(s + d2) and N = Delta W = n0 + n1 s + n2 s^2 + n3 s^3, the moments
        U_k = int_0^inf s^k e^-s / N(s) ds
    for k = 0, 1, 2 from the roots of N (s^k / N = sum_i (-lambda_i)^k / (N'(-lambda_i)(s + lambda_i)) for k <= 2)
    and for k >= 3 from the recurrence n0 U_k + n1 U_{k+1} + n2 U_{k+2} + n3 U_{k+3} = k!; then
        int s^k e^-s / W(s) ds = sum_j Delta_j U_{k+j};  the same with N + 2 rho Delta for 2 rho + W;
  * F(x) = e^x E1(x) = int_0^inf e^-t/(x + t) dt by the continued fraction 1/(x + 1/(1 + 1/(x + 2/(1 + 2/(x + ...)))))
    (Abramowitz-Stegun 5.1.22), whose successive convergents lie on alternate sides of F(x) for x > 0; each
    convergent is evaluated in interval arithmetic, and the hull of two successive ones encloses F(x);
  * roots of the cubics by interval bisection in the brackets (-inf, -d2), (-d2, -d1), (-d1, 0);
  * Gram blocks: interval LDL without pivoting in increasing index order (to compare pivots with program A and the
    paper); S as the last pivot of the LDL of the bordered matrix [[B^-1 + G, b], [b^T, R]].

  c10b_gfl_b.py [--set NAME VALUE]... [--claim-ratio DEC] [--claim-mu DEC]
Standard library only; run with python3 -I.
"""
import sys
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING
from math import comb, factorial

P = 100
DN = Context(prec=P, rounding=ROUND_FLOOR)
UP = Context(prec=P, rounding=ROUND_CEILING)

PARAM_TEXT = """rho 0.012039377730497942
beta1 0.4283744445791419
beta2 0.1506077793447579
beta3 0.08328691492582467
beta4 0.05294104124824108
beta5 0.03632381026351961
beta6 0.026165030774068086
beta7 0.01945831930582338
beta8 0.01490133597851163
beta9 0.011682554694871422
beta10 0.014696581084273385
alpha 0.16212090047056932029265504
eta 0.2995274614735440217489782
zeta1 0.1705006229073915451327806
delta1 0.6122395153199937
zeta2 1.3819822296622917537966546
delta2 6.340370329139007"""


class Iv:
    __slots__ = ("a", "b")

    def __init__(self, a, b=None):
        self.a = D(a)
        self.b = D(a) if b is None else D(b)
        if self.a > self.b:
            raise ValueError("empty interval")

    @staticmethod
    def of(x):
        return x if isinstance(x, Iv) else Iv(x)

    def __add__(s, o):
        o = Iv.of(o); return Iv(DN.add(s.a, o.a), UP.add(s.b, o.b))
    __radd__ = __add__

    def __sub__(s, o):
        o = Iv.of(o); return Iv(DN.subtract(s.a, o.b), UP.subtract(s.b, o.a))

    def __rsub__(s, o):
        return Iv.of(o) - s

    def __neg__(s):
        return Iv(s.b.copy_negate(), s.a.copy_negate())     # exact (unary minus would round)

    def __mul__(s, o):
        o = Iv.of(o)
        c = [(x, y) for x in (s.a, s.b) for y in (o.a, o.b)]
        return Iv(min(DN.multiply(x, y) for x, y in c), max(UP.multiply(x, y) for x, y in c))
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = Iv.of(o)
        if o.a <= 0 <= o.b:
            raise ZeroDivisionError("interval contains 0")
        c = [(x, y) for x in (s.a, s.b) for y in (o.a, o.b)]
        return Iv(min(DN.divide(x, y) for x, y in c), max(UP.divide(x, y) for x, y in c))

    def __rtruediv__(s, o):
        return Iv.of(o) / s

    def sqrt(s):
        return Iv(DN.sqrt(s.a), UP.sqrt(s.b))

    def hull(s, o):
        return Iv(min(s.a, o.a), max(s.b, o.b))

    def width(s):
        return UP.subtract(s.b, s.a)


def frac(n, d):
    """The rational n/d (integers) as a tight interval."""
    return Iv(DN.divide(n, d), UP.divide(n, d))


def F_cf(x):
    """e^x E1(x) for an interval x > 0 by the continued fraction; doubles the depth until two successive convergents
    agree to 1e-85, and returns their hull."""
    def conv(n):
        # 1/(x + 1/(1 + 1/(x + 2/(1 + 2/(x + 3/(1 + ...)))))): level j has the partial denominator x (j odd) or 1
        # (j even) and the partial numerator 1 (j = 1) or floor(j/2) (j >= 2); the tail below level n is cut to 0.
        t = Iv(0)
        for j in range(n, 0, -1):
            den = x if j % 2 == 1 else Iv(1)
            num = 1 if j == 1 else j // 2
            t = Iv(num) / (den + t)
        return t
    n = 64
    while True:
        c1, c2 = conv(n), conv(n + 1)
        h = c1.hull(c2)
        if h.width() < D("1e-85"):
            # alternation check on a few more convergents: c(n+2) must lie inside the hull of c(n), c(n+1)
            c3 = conv(n + 2)
            assert c3.a >= h.a and c3.b <= h.b, "convergents do not alternate"
            return h
        n *= 2


def poly_mul(p, q):
    out = [Iv(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] = out[i + j] + a * b
    return out


def peval(c, x):
    acc = Iv(0)
    for a in reversed(c):
        acc = acc * x + a
    return acc


def roots(N, d1, d2):
    out = []
    for lo, hi in ((D(-1000), d2.copy_negate()), (d2.copy_negate(), d1.copy_negate()), (d1.copy_negate(), D(0))):
        flo, fhi = peval(N, Iv(lo)), peval(N, Iv(hi))
        neg_at_lo = flo.b < 0
        assert (flo.b < 0 < fhi.a) or (fhi.b < 0 < flo.a), "no sign change"
        while hi - lo > D("1e-90"):
            mid = DN.divide(DN.add(lo, hi), 2)
            fm = peval(N, Iv(mid))
            if fm.b < 0:
                if neg_at_lo: lo = mid
                else: hi = mid
            elif fm.a > 0:
                if neg_at_lo: hi = mid
                else: lo = mid
            else:
                break
        out.append(Iv(lo, hi))
    return out


def moments(N, Delta, d1, d2, kmax):
    """int s^k e^-s Delta/N ds for k = 0..kmax."""
    rts = roots(N, d1, d2)
    dN = [N[1], N[2] * 2, N[3] * 3]
    Fv = [F_cf(-r) for r in rts]
    U = []
    for k in range(3):
        u = Iv(0)
        for r, f in zip(rts, Fv):
            pw = Iv(1)
            for _ in range(k):
                pw = pw * r
            u = u + pw * f / peval(dN, r)
        U.append(u)
    for k in range(0, kmax + 2):
        U.append((Iv(factorial(k)) - N[0] * U[k] - N[1] * U[k + 1] - N[2] * U[k + 2]) / N[3])
    return [Delta[0] * U[k] + Delta[1] * U[k + 1] + Delta[2] * U[k + 2] for k in range(kmax + 1)], rts


def laguerre_e(a, b):
    """e_{a,b}(r) as (sign, ratio, poly) with e = sign * sqrt(ratio) * poly(r): poly with Fraction-free rational coefs
    represented as (numerator list over a common factorial denominator)."""
    if a >= b:
        al, j, sign, num, den = a - b, b, (-1) ** b, factorial(b), factorial(a)
    else:
        al, j, sign, num, den = b - a, a, (-1) ** a, factorial(a), factorial(b)
    # r^al L_j^(al)(r^2) = sum_k (-1)^k C(j+al, j-k) r^(al+2k) / k!; common denominator j!
    coefs = {}
    for k in range(j + 1):
        coefs[al + 2 * k] = (-1) ** k * comb(j + al, j - k) * (factorial(j) // factorial(k))
    return sign, (num, den), coefs, factorial(j)


def product_integral(e1, e2, mom):
    """int e1 e2 (weight) d nu, as an interval, from the moment list mom[k] = int s^k e^-s (weight) ds."""
    s1, (n1, d1), c1, q1 = e1
    s2, (n2, d2), c2, q2 = e2
    acc = Iv(0)
    prod = {}
    for i, a in c1.items():
        for j, b in c2.items():
            prod[i + j] = prod.get(i + j, 0) + a * b
    for deg, c in prod.items():
        assert deg % 2 == 0
        if c:
            acc = acc + Iv(c) * mom[deg // 2]
    root = frac(n1 * n2, d1 * d2).sqrt()
    return acc * root * (s1 * s2) / Iv(q1 * q2)


def ldl(A, stop=True):
    """Pivots of an interval LDL^T without pivoting; with stop, it stops at the first pivot not certified positive."""
    n = len(A)
    L = [[None] * n for _ in range(n)]
    d = []
    for j in range(n):
        v = A[j][j]
        for k in range(j):
            v = v - L[j][k] * L[j][k] * d[k]
        d.append(v)
        if v.a <= 0:
            if stop or j < n - 1:
                return d, False
        for i in range(j + 1, n):
            w = A[i][j]
            for k in range(j):
                w = w - L[i][k] * L[j][k] * d[k]
            L[i][j] = w / v
    return d, True


def show(x, d=25):
    q = D(10) ** -d
    c1 = Context(prec=300, rounding=ROUND_FLOOR); c2 = Context(prec=300, rounding=ROUND_CEILING)
    return "[%s, %s]" % (c1.quantize(x.a, q), c2.quantize(x.b, q))


def main():
    p = {}
    for line in PARAM_TEXT.splitlines():
        k, v = line.split()
        p[k] = D(v)
    claim_ratio = claim_mu = None
    a = sys.argv[1:]
    while a:
        if a[0] == "--set":
            p[a[1]] = D(a[2]); print("CONTROL %s set to %s" % (a[1], a[2])); a = a[3:]
        elif a[0] == "--claim-ratio":
            claim_ratio = D(a[1]); a = a[2:]
        elif a[0] == "--claim-mu":
            claim_mu = D(a[1]); a = a[2:]
        else:
            raise SystemExit("unknown option " + a[0])
    rho = Iv(p["rho"])
    d1, d2 = p["delta1"], p["delta2"]
    # self-test: the Laguerre form against the moments of nu (orthonormality), exactly enough at 100 digits
    for aa in range(6):
        for bb in range(6):
            for cc in range(6):
                dd = cc - (aa - bb)
                if dd < 0:
                    continue
                v = product_integral(laguerre_e(aa, bb), laguerre_e(cc, dd), [Iv(factorial(k)) for k in range(30)])
                target = 1 if (aa, bb) == (cc, dd) else 0
                assert v.a <= target <= v.b and v.width() < D("1e-80"), (aa, bb, cc, dd)
    print("self-test: the Laguerre forms of e_{a,b} (a, b <= 5) are orthonormal in L2(nu)")
    print("e^1 E1(1) in %s" % show(F_cf(Iv(1)), 30))
    Delta = [Iv(d1) * Iv(d2), Iv(d1) + Iv(d2), Iv(1)]
    Wlin = [Iv(p["alpha"]), Iv(p["eta"])]
    N = poly_mul(Wlin, Delta)
    t1 = poly_mul([Iv(0), Iv(1)], [Iv(d2), Iv(1)])
    t2 = poly_mul([Iv(0), Iv(1)], [Iv(d1), Iv(1)])
    N = [N[i] + Iv(p["zeta1"]) * (t1[i] if i < 3 else 0) + Iv(p["zeta2"]) * (t2[i] if i < 3 else 0) for i in range(4)]
    ND = [N[i] + (Iv(2) * rho * Delta[i] if i < 3 else Iv(0)) for i in range(4)]
    momW, rW = moments(N, Delta, d1, d2, 21)
    momD, rD = moments(ND, Delta, d1, d2, 21)
    print("poles of 1/W at -lambda, lambda in %s" % ", ".join(show(-r, 12) for r in rW))
    print("poles of 1/(2 rho + W) at -lambda, lambda in %s" % ", ".join(show(-r, 12) for r in rD))
    print("int e^-s / W(s) ds in %s" % show(momW[0]))
    beta = {l: Iv(p["beta%d" % l]) for l in range(1, 11)}
    minpiv, minblk, npiv, nblk, allok = None, None, 0, 0, True
    for pp in range(12):
        for qq in range(11):
            I = [l for l in range(1, 11) if pp <= l + 1 and qq <= l]
            if not I:
                continue
            nblk += 1
            es = {l: laguerre_e(l + 1 - pp, l - qq) for l in I}
            A = [[(Iv(1) if l == m else Iv(0)) - (beta[l] * beta[m]).sqrt() * product_integral(es[l], es[m], momW)
                  for m in I] for l in I]
            d, ok = ldl(A)
            allok &= ok
            npiv += len(d)
            for v in d:
                if minpiv is None or v.a < minpiv.a:
                    minpiv, minblk = v, (pp, qq)
    print("Gram blocks: %d nonempty, %d pivots, all positive: %s; least pivot in %s, block (p,q) = %s"
          % (nblk, npiv, allok, show(minpiv, 15), minblk))
    phi = {l: laguerre_e(l + 1, l) for l in range(1, 11)}
    rfun = laguerre_e(1, 0)                         # e_{1,0}(r) = r
    R = product_integral(rfun, rfun, momD)
    bv = [product_integral(rfun, phi[l], momD) for l in range(1, 11)]
    Bord = [[product_integral(phi[l], phi[m], momD) + ((Iv(1) / beta[l]) if l == m else Iv(0)) for m in range(1, 11)]
            + [bv[l - 1]] for l in range(1, 11)]
    Bord.append(bv + [R])
    d, ok = ldl([row[:10] for row in Bord[:10]])
    if not ok:
        print("B^-1 + G not certified positive definite"); print("VERDICT REFUSED"); return
    dfull, _ = ldl(Bord, stop=False)              # its 11th pivot is R - b^T (B^-1 + G)^-1 b = S
    S = dfull[10]
    print("B^-1 + G positive definite (LDL): True")
    print("R in %s" % show(R))
    print("S in %s" % show(S))
    print("1 - S in %s" % show(Iv(1) - S, 20))
    mu = Iv(p["alpha"]) + Iv(p["eta"])
    for j in (1, 2):
        dl = Iv(p["delta%d" % j])
        mu = mu + Iv(p["zeta%d" % j]) * (Iv(1) - dl * F_cf(dl))
    ratio = (Iv(1) - rho) / (rho + mu)
    print("mu_W in %s" % show(mu))
    print("ratio (1 - rho)/(rho + mu_W) in %s" % show(ratio))
    okall = allok and (Iv(1) - S).a > 0
    print("conditions of Theorem 3.1 (all Gram blocks, S < 1): %s" % ("hold" if okall else "FAIL"))
    if claim_mu is not None:
        h = mu.b < claim_mu
        print("claimed mu_W < %s: %s" % (claim_mu, "holds" if h else "FAILS")); okall &= h
    if claim_ratio is not None:
        h = ratio.a > claim_ratio
        print("claimed ratio > %s: %s" % (claim_ratio, "holds" if h else "FAILS")); okall &= h
    print("VERDICT %s" % ("PASS" if okall else "REFUSED"))


if __name__ == "__main__":
    main()
