#!/usr/bin/env python3
"""8a, pull request #373: the finite inputs of the deduction C_8 <= 2/log(gamma_1) ~ 0.755107*, program B (written
separately from program A).

Same three statements as program A: (1) Hardy's Z changes sign on (14.13472514173469378, 14.13472514173469380);
(2) N(T0) = 1 at T0 = 89/5, from N(T) = theta(T)/pi + 1 + S(T) (the Riemann-von Mangoldt identity, read) and
Re zeta(sigma + i T0) > 0 for 1/2 <= sigma <= 2, which gives |S(T0)| < 1/2; (3) 2/log gamma_1 <= 0.755107,
gamma_1^4 < 4 * 10^4 and 8/log(3 * 10^12) < 0.279. [OAI2026]'s theorem and [PT2021]'s verification are not checked.

Different methods from A throughout:
- zeta(s) = eta(s)/(1 - 2^(1-s)) with eta summed by Algorithm 1 of Cohen, Rodriguez Villegas and Zagier (n = 90,
  exact rational weights), sigma taken as a decimal interval ((k+1)^-sigma = exp(-sigma ln(k+1)), monotone). Error
  bound for every sigma >= 1/2, derived here: the algorithm's error is at most 2 Gamma(sigma)/(|Gamma(s)| (3+sqrt 8)^n)
  (see the c8a docstring), and from Gamma(x+iy)/Gamma(x) = prod_{n>=0} (1 + y^2/(x+n)^2)^(-1/2) each factor is
  smallest at x = 1/2 when x >= 1/2, so Gamma(sigma)/|Gamma(sigma+it)| <= Gamma(1/2)/|Gamma(1/2+it)| = sqrt(cosh(pi t)):
  |error| <= 2 sqrt(cosh(pi t))/(3+sqrt 8)^n <= 2 exp(pi t/2)/5.828^n, the same bound as on the critical line.
- Re zeta on [1/2, 2] at height T0: a different first split (40 pieces) and bisection until each lower end is
  positive. theta(T0) by Stirling's series at w = z + 40 with 15 terms (shift 30 for Z in (1)); Bernoulli numbers by Akiyama-Tanigawa.
- Arithmetic: decimal intervals at 110 digits with directed rounding; ln, exp and sqrt from the decimal module
  (correctly rounded, then moved outward by one unit), cos/sin/atan of our own with explicit remainders, pi by BBP.
Standard library only.
"""
import sys
from fractions import Fraction as Q
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING, getcontext

PR = 110
getcontext().prec = PR
LO = Context(prec=PR, rounding=ROUND_FLOOR)
HI = Context(prec=PR, rounding=ROUND_CEILING)
TINY = D(10) ** -(PR - 8)


def dq(q, ctx):
    q = Q(q)
    return ctx.divide(D(q.numerator), D(q.denominator))


class V:
    """decimal interval [a, b]"""
    __slots__ = ("a", "b")

    def __init__(self, a, b=None):
        self.a = a; self.b = a if b is None else b

    @staticmethod
    def of(q):
        return V(dq(q, LO), dq(q, HI))

    def __add__(x, y):
        y = y if isinstance(y, V) else V.of(y)
        return V(LO.add(x.a, y.a), HI.add(x.b, y.b))

    def __sub__(x, y):
        y = y if isinstance(y, V) else V.of(y)
        return V(LO.subtract(x.a, y.b), HI.subtract(x.b, y.a))

    def __neg__(x): return V(x.b.copy_negate(), x.a.copy_negate())

    def __mul__(x, y):
        y = y if isinstance(y, V) else V.of(y)
        lo = min(LO.multiply(x.a, y.a), LO.multiply(x.a, y.b), LO.multiply(x.b, y.a), LO.multiply(x.b, y.b))
        hi = max(HI.multiply(x.a, y.a), HI.multiply(x.a, y.b), HI.multiply(x.b, y.a), HI.multiply(x.b, y.b))
        return V(lo, hi)

    def __truediv__(x, y):
        y = y if isinstance(y, V) else V.of(y)
        assert y.a > 0 or y.b < 0
        c = [LO.divide(x.a, y.a), LO.divide(x.a, y.b), LO.divide(x.b, y.a), LO.divide(x.b, y.b)]
        e = [HI.divide(x.a, y.a), HI.divide(x.a, y.b), HI.divide(x.b, y.a), HI.divide(x.b, y.b)]
        return V(min(c), max(e))

    def mag(x): return max(x.a.copy_abs(), x.b.copy_abs())
    def pad(x, e): return V(LO.subtract(x.a, e), HI.add(x.b, e))
    def has0(x): return x.a <= 0 <= x.b
    def width(x): return HI.subtract(x.b, x.a)


def v_ln(x):    # x > 0 interval; ln is monotone and correctly rounded by the decimal module
    return V(LO.next_minus(x.a.ln()), HI.next_plus(x.b.ln()))


def v_sqrt(x):
    return V(LO.next_minus(x.a.sqrt()), HI.next_plus(x.b.sqrt()))


def v_pi():
    s = Q(0)
    n = 120
    for k in range(n):
        s += Q(1, 16 ** k) * (Q(4, 8 * k + 1) - Q(2, 8 * k + 4) - Q(1, 8 * k + 5) - Q(1, 8 * k + 6))
    return V(dq(s, LO), dq(s + Q(4, 16 ** n) * Q(16, 15), HI))


PI = v_pi()


def v_cossin(x):
    k = int((x.a / (2 * PI.a)).to_integral_value())
    y = x - PI * V.of(2 * k)
    c = V.of(1); s = V.of(0); term = V.of(1); n = 0
    while True:
        n += 1
        term = term * y / V.of(n)
        r = n % 4
        if r == 1: s = s + term
        elif r == 2: c = c - term
        elif r == 3: s = s - term
        else: c = c + term
        if n > 8 and term.mag() < TINY:
            break
    e = HI.multiply(term.mag(), y.mag())
    return c.pad(e), s.pad(e)


def v_atan(x):   # x > 0 interval; halving atan(x) = 2 atan(x/(1+sqrt(1+x^2))) until x < 0.2, then Taylor
    m = 1
    while x.b >= D("0.2"):
        x = V(LO.divide(x.a, HI.add(1, HI.next_plus(HI.add(1, HI.multiply(x.a, x.a)).sqrt()))),
              HI.divide(x.b, LO.add(1, LO.next_minus(LO.add(1, LO.multiply(x.b, x.b)).sqrt()))))
        m *= 2
    s = V.of(0); p = x; x2 = x * x; j = 0
    while True:
        t = p / V.of(2 * j + 1)
        s = s + t if j % 2 == 0 else s - t
        j += 1
        p = p * x2
        if p.mag() < TINY:
            break
    return (s.pad(p.mag())) * V.of(m)


def bern_at(n):
    """B_0..B_n by Akiyama-Tanigawa (B_1 = +1/2 there; only even indices are used)"""
    out = []
    a = [Q(0)] * (n + 1)
    for m in range(n + 1):
        a[m] = Q(1, m + 1)
        for j in range(m, 0, -1):
            a[j - 1] = j * (a[j - 1] - a[j])
        out.append(a[0])
    return out


BB = bern_at(40)


class Z2:
    __slots__ = ("r", "i")
    def __init__(self, r, i): self.r, self.i = r, i
    def __add__(x, y): return Z2(x.r + y.r, x.i + y.i)
    def __sub__(x, y): return Z2(x.r - y.r, x.i - y.i)
    def __mul__(x, y):
        if isinstance(y, Z2):
            return Z2(x.r * y.r - x.i * y.i, x.r * y.i + x.i * y.r)
        return Z2(x.r * y, x.i * y)
    def recip(x):
        n = x.r * x.r + x.i * x.i
        return Z2(x.r / n, -x.i / n)


def crvz_weights(n):
    """Algorithm 1 of Cohen-Rodriguez Villegas-Zagier: S_n = sum_k w_k a_k / d with exact rationals"""
    T0, T1 = 1, 3
    for _ in range(n - 1):
        T0, T1 = T1, 6 * T1 - T0
    d = T1 if n >= 1 else 1
    b = Q(-1); c = Q(-d); w = []
    for k in range(n):
        c = b - c
        w.append(c)
        b = Q((k + n) * (k - n)) * b / ((Q(k) + Q(1, 2)) * (k + 1))
    return w, d


def eta_half(t, n):
    """eta(1/2 + it) by CRVZ with the error bound of the docstring"""
    w, d = crvz_weights(n)
    tt = V.of(t)
    acc = Z2(V.of(0), V.of(0))
    for k in range(n):
        m = k + 1
        lnm = v_ln(V(D(m)))
        c, s = v_cossin(lnm * tt)
        r = V.of(1) / v_sqrt(V(D(m)))
        acc = acc + Z2(c * r, -(s * r)) * V.of(Q(w[k], d))
    # error <= 2 sqrt(cosh(pi t)) / (3 + sqrt 8)^n <= 2 exp(pi t / 2) / 5.828^n
    e_pt = HI.next_plus(HI.multiply(PI.b, HI.divide(D(t.numerator), D(2 * t.denominator))).exp())
    bound = HI.divide(HI.multiply(2, e_pt), LO.power(D("5.828"), n))
    return Z2(acc.r.pad(bound), acc.i.pad(bound)), bound


def im_lngamma(zr, zi, K=30, m=15):
    """Im log Gamma(zr + i zi) up to a multiple of 2 pi; zr > 0"""
    def arg(x, y):     # 0 < y/x: atan(y/x)
        return v_atan(V.of(Q(y) / Q(x)))
    tot = V.of(0)
    for j in range(K):
        tot = tot - arg(zr + j, zi)
    wr, wi = zr + K, zi
    W = Z2(V.of(wr), V.of(wi))
    lnabs = v_ln(V.of(wr * wr + wi * wi)) / V.of(2)
    aw = arg(wr, wi)
    # Im[(w - 1/2) log w - w + log(2 pi)/2] = (wr - 1/2) arg w + wi ln|w| - wi
    tot = tot + aw * V.of(wr - Q(1, 2)) + lnabs * V.of(wi) - V.of(wi)
    inv = W.recip(); p = inv; inv2 = inv * inv
    for k in range(1, m + 1):
        tot = tot + (p * V.of(BB[2 * k] / (2 * k * (2 * k - 1)))).i
        p = p * inv2
    # remainder |B_{2m+2}|/((2m+2)(2m+1)|w|^(2m+1)) sec^(2m+2)(arg w/2); |w| >= wr, arg w < 0.3, cos(0.15) > 0.988
    assert Q(wi) / Q(wr) < Q(3, 10)
    rem = abs(BB[2 * m + 2]) / ((2 * m + 2) * (2 * m + 1)) / Q(wr) ** (2 * m + 1) / Q(988, 1000) ** (2 * m + 2)
    return tot.pad(dq(rem, HI))


def hardy_z(t, n=90, ctrl=None):
    th = im_lngamma(Q(1, 4), t / 2) - v_ln(PI) * V.of(t / 2)
    if ctrl == "theta-sign":
        th = -th
    c, s = v_cossin(th)
    eta, bnd = eta_half(t, n)
    # 1 - 2^(1-s) = 1 - sqrt 2 (cos(t ln 2) - i sin(t ln 2))
    c2, s2 = v_cossin(v_ln(V(D(2))) * V.of(t))
    r2 = v_sqrt(V(D(2)))
    den = Z2(V.of(1) - r2 * c2, r2 * s2)
    zeta = eta * den.recip()
    zz = Z2(c, s) * zeta
    return zz.r, zz.i, bnd




def v_exp_neg(slo, shi, lnm):
    """exp(-sigma ln m) for sigma in [slo, shi] (decimal), ln m an interval: decreasing in sigma"""
    a = LO.next_minus(LO.multiply(-shi, lnm.b).exp(LO))
    b = HI.next_plus(HI.multiply(-slo, lnm.a).exp(HI))
    return V(a, b)


LNCACHE = {}
CSCACHE = {}


def zeta_strip(slo, shi, t, n=90):
    """zeta(sigma + it), sigma in [slo, shi] (Decimals, slo >= 1/2), t a Fraction: CRVZ for eta, then the divisor"""
    w, d = crvz_weights(n)
    tt = V.of(t)
    acc = Z2(V.of(0), V.of(0))
    for k in range(n):
        m = k + 1
        if m not in LNCACHE:
            LNCACHE[m] = v_ln(V(D(m)))
        lnm = LNCACHE[m]
        if (m, t) not in CSCACHE:
            CSCACHE[(m, t)] = v_cossin(lnm * tt)
        c, s = CSCACHE[(m, t)]
        r = v_exp_neg(slo, shi, lnm) if m > 1 else V.of(1)
        acc = acc + Z2(c * r, -(s * r)) * V.of(Q(w[k], d))
    e_pt = HI.next_plus(HI.multiply(PI.b, HI.divide(D(t.numerator), D(2 * t.denominator))).exp())
    bound = HI.divide(HI.multiply(2, e_pt), LO.power(D("5.828"), n))
    eta = Z2(acc.r.pad(bound), acc.i.pad(bound))
    # 1 - 2^(1-s) = 1 - 2^(1-sigma) (cos(t ln 2) - i sin(t ln 2)); 2^(1-sigma) = 2 exp(-sigma ln 2)
    ln2 = v_ln(V(D(2)))
    c2, s2 = v_cossin(ln2 * tt)
    p2 = v_exp_neg(slo, shi, ln2) * V.of(2)
    den = Z2(V.of(1) - p2 * c2, p2 * s2)
    return eta * den.recip(), bound


def theta(t):
    return im_lngamma(Q(1, 4), t / 2, K=40) - v_ln(PI) * V.of(t / 2)


def re_positive_on_segment(t, pieces=40, maxdepth=14):
    lo0, hi0 = Q(1, 2), Q(2)
    queue = [(lo0 + (hi0 - lo0) * k / pieces, lo0 + (hi0 - lo0) * (k + 1) / pieces, 0) for k in range(pieces)]
    count = 0; least = None
    while queue:
        a, b, dep = queue.pop(0)
        z, bnd = zeta_strip(dq(a, LO), dq(b, HI), t)
        if z.r.a > 0:
            count += 1
            least = z.r.a if least is None or z.r.a < least else least
            continue
        if z.r.b < 0 or dep >= maxdepth:
            return False, (a, b, z.r.a, z.r.b), count
        mid = (a + b) / 2
        queue += [(a, mid, dep + 1), (mid, b, dep + 1)]
    return True, least, count


def main():
    a = sys.argv[1:]
    ctrl = a[a.index("--control") + 1] if "--control" in a else None
    T0 = Q(89, 5)
    claim_N = 1
    printed = D("0.755107")
    t_lo, t_hi = Q("14.13472514173469378"), Q("14.13472514173469380")
    if ctrl == "height22":
        T0 = Q(22)
    if ctrl == "printed":
        printed = D("0.755106")
    if ctrl == "bracket":
        t_lo, t_hi = Q("14.13472514173469381"), Q("14.13472514173469383")
    print("program B (CRVZ for eta with sigma as an interval, n = 90; Stirling shift 30, 15 terms; decimal %d digits)%s"
          % (PR, (" CONTROL " + ctrl) if ctrl else ""))
    w, d = crvz_weights(40)
    l2 = sum(Q(w[k], d) / (k + 1) for k in range(40))
    ln2 = v_ln(V(D(2)))
    gap = max(abs(dq(l2, LO) - ln2.a), abs(dq(l2, HI) - ln2.b))
    assert gap < D("5.828") ** -40 * 2, "CRVZ on log 2"
    assert D("3.14159265358979323846264338327") < PI.a and PI.b < D("3.14159265358979323846264338328")
    z, _ = zeta_strip(D(2), D(2), Q(1, 10 ** 30))
    pi2 = PI * PI / V.of(6)
    assert z.r.a - D(10) ** -60 <= pi2.b and pi2.a <= z.r.b + D(10) ** -60, "zeta(2) self-test"
    print("self-tests: CRVZ on log 2, pi, zeta(2 + 10^-30 i) = pi^2/6 within 10^-60: ok")
    ok = True
    vals = []
    for t in (t_lo, t_hi):
        zr, zi, bnd = hardy_z(t)
        ok = ok and zi.has0()
        vals.append(zr)
    change = (vals[0].b < 0 < vals[1].a) or (vals[0].a > 0 > vals[1].b)
    print("(1) Z changes sign on (%s, %s): %s" % (D(t_lo.numerator) / D(t_lo.denominator),
                                                  D(t_hi.numerator) / D(t_hi.denominator), change))
    ok = ok and change
    good, info, count = re_positive_on_segment(T0)
    if good:
        print("(2) Re zeta(sigma + %s i) > 0 on [1/2, 2]: yes (%d pieces, least lower end %s), so |S(T0)| < 1/2"
              % (D(T0.numerator) / D(T0.denominator), count, str(info)[:10]))
    else:
        print("(2) Re zeta > 0 on [1/2, 2]: NOT SHOWN at [%s, %s]" % (info[0], info[1]))
    x = theta(T0) / PI + V.of(1)
    cands = [k for k in range(-2, 10) if x.a - D("0.5") < k < x.b + D("0.5")]
    nT = cands[0] if good and len(cands) == 1 else None
    print("    theta(T0)/pi + 1 in [%s, %s]; integers within 1/2: %s; N(T0) = %s; the claim N = %d: %s"
          % (str(x.a)[:16], str(x.b)[:16], cands, nT, claim_N, "holds" if nT == claim_N else "FAILS"))
    ok = ok and good and nT == claim_N
    lo = LO.divide(2, HI.next_plus(dq(t_hi, HI).ln(HI)))
    hi = HI.divide(2, LO.next_minus(dq(t_lo, LO).ln(LO)))
    print("(3) 2/log gamma_1 in [%s, %s]; printed upper bound %s: %s" % (str(lo)[:24], str(hi)[:24], printed,
          "holds" if hi <= printed else "FAILS"))
    ok = ok and hi <= printed
    g4 = HI.power(dq(t_hi, HI), 4)
    e8 = HI.divide(8, LO.next_minus(D(3 * 10 ** 12).ln(LO)))
    print("    gamma_1^4 <= %s < 40000: %s; 8/log(3*10^12) <= %s < 0.279 < 2/log gamma_1: %s"
          % (str(g4)[:12], g4 < 40000, str(e8)[:14], e8 < D("0.279") < lo))
    ok = ok and g4 < 40000 and e8 < D("0.279") < lo
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if ctrl else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
