"""Rigorous interval arithmetic for the 22b replay (our own code, standard library only).

Intervals are pairs of binary64 floats. Every +, -, *, / and sqrt is rounded to nearest by IEEE 754 and then widened
outward by one ulp (math.nextafter), which encloses the exact result. sin and cos of a float point are computed with
Python's decimal module at 60 digits (Taylor series after reduction by pi/2, with pi enclosed by Machin's formula in
exact rationals) and an error allowance of 1e-45, far above every rounding and truncation error; over an interval
they follow from sin(m + d) = sin m cos d + cos m sin d with |d| <= r, cos d in [1 - r^2/2, 1], sin d in [-r, r].
So no libm function is trusted. Jets (truncated Taylor series with interval coefficients) give enclosures of a
function and its first derivatives over an interval of the parameter (natural interval extension of the Taylor
recurrences).
"""
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import lru_cache
from math import inf, nextafter, sqrt as fsqrt

getcontext().prec = 60


def dn(x):
    return nextafter(x, -inf)


def up(x):
    return nextafter(x, inf)


def pi_bounds():
    """exact rational bounds for pi from Machin's formula pi = 16 atan(1/5) - 4 atan(1/239), alternating series."""
    def atan_bounds(x, n):   # partial sums S_n, S_{n+1} bracket atan(x) for 0 < x < 1
        s, terms = Fraction(0), []
        for k in range(n + 2):
            terms.append(Fraction((-1) ** k, 2 * k + 1) * x ** (2 * k + 1))
        s1 = sum(terms[:n + 1])
        s2 = s1 + terms[n + 1]
        return min(s1, s2), max(s1, s2)
    a_lo, a_hi = atan_bounds(Fraction(1, 5), 45)
    b_lo, b_hi = atan_bounds(Fraction(1, 239), 15)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo


PI_LO, PI_HI = pi_bounds()
assert PI_HI - PI_LO < Fraction(1, 10 ** 60)
assert Fraction(314159265358979323846, 10 ** 20) < PI_LO < PI_HI < Fraction(314159265358979323847, 10 ** 20)
HALF_PI_D = (Decimal(PI_LO.numerator) / Decimal(PI_LO.denominator)) / 2   # error below 1e-59


class Iv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi

    def __repr__(self):
        return f"[{self.lo!r}, {self.hi!r}]"

    def __add__(self, o):
        if isinstance(o, Iv):
            return Iv(dn(self.lo + o.lo), up(self.hi + o.hi))
        return Iv(dn(self.lo + o), up(self.hi + o))   # o: a float that is exact data (used only with exact floats)

    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        if isinstance(o, Iv):
            return Iv(dn(self.lo - o.hi), up(self.hi - o.lo))
        return Iv(dn(self.lo - o), up(self.hi - o))

    def __rsub__(self, o):
        return Iv(dn(o - self.hi), up(o - self.lo))

    def __mul__(self, o):
        if isinstance(o, Iv):
            a, b, c, d = self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi
            return Iv(dn(min(a, b, c, d)), up(max(a, b, c, d)))
        a, b = self.lo * o, self.hi * o
        return Iv(dn(min(a, b)), up(max(a, b)))

    __rmul__ = __mul__

    def inv(self):
        if self.lo <= 0.0 <= self.hi:
            raise ZeroDivisionError("interval contains 0")
        a, b = 1.0 / self.lo, 1.0 / self.hi
        return Iv(dn(min(a, b)), up(max(a, b)))

    def __truediv__(self, o):
        if isinstance(o, Iv):
            if o.lo <= 0.0 <= o.hi:
                raise ZeroDivisionError("interval contains 0")
            a, b, c, d = self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi
            return Iv(dn(min(a, b, c, d)), up(max(a, b, c, d)))
        a, b = self.lo / o, self.hi / o
        return Iv(dn(min(a, b)), up(max(a, b)))

    def sqr(self):
        a, b = self.lo * self.lo, self.hi * self.hi
        if self.lo >= 0.0:
            return Iv(dn(a), up(b))
        if self.hi <= 0.0:
            return Iv(dn(b), up(a))
        return Iv(0.0, up(max(a, b)))

    def sqrt(self):
        if self.lo < 0.0:
            raise ValueError("sqrt of an interval with a negative part")
        return Iv(dn(fsqrt(self.lo)) if self.lo > 0.0 else 0.0, up(fsqrt(self.hi)))

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def mid(self):
        return 0.5 * (self.lo + self.hi)

    def width(self):
        return self.hi - self.lo

    def hull(self, o):
        return Iv(min(self.lo, o.lo), max(self.hi, o.hi))

    def meets(self, o):
        return self.lo <= o.hi and o.lo <= self.hi


ZERO = Iv(0.0)
ONE = Iv(1.0)


def ivq(q):
    """the tightest float interval containing the rational q"""
    q = Fraction(q)
    x = float(q)
    fx = Fraction(x)
    return Iv(x if fx <= q else dn(x), x if fx >= q else up(x))


def ivq_hull(q1, q2):
    a, b = ivq(min(q1, q2)), ivq(max(q1, q2))
    return Iv(a.lo, b.hi)


PI = ivq_hull(PI_LO, PI_HI)
ERR = Decimal("1e-45")


@lru_cache(maxsize=1 << 18)
def sincos_point(m):
    """rigorous float intervals for sin(m) and cos(m), m a float (exact)"""
    d = Decimal(m)
    k = int(round(float(d / HALF_PI_D)))
    if abs(k) > 10 ** 6:
        raise ValueError("argument too large")
    r = d - k * HALF_PI_D
    r2 = r * r
    s, t, j = r, r, 1
    while True:   # sin r = sum (-1)^j r^(2j+1)/(2j+1)!
        t = -t * r2 / ((2 * j) * (2 * j + 1))
        s += t
        j += 1
        if abs(t) < Decimal("1e-58"):
            break
    c, t, j = Decimal(1), Decimal(1), 1
    while True:
        t = -t * r2 / ((2 * j - 1) * (2 * j))
        c += t
        j += 1
        if abs(t) < Decimal("1e-58"):
            break
    q = k % 4
    if q == 0:
        sv, cv = s, c
    elif q == 1:
        sv, cv = c, -s
    elif q == 2:
        sv, cv = -s, -c
    else:
        sv, cv = -c, s

    def iv(x):
        return Iv(dn(float(x - ERR)), up(float(x + ERR)))
    return iv(sv), iv(cv)


def sincos(x):
    """enclosures of sin and cos over the interval x"""
    if x.lo == x.hi:
        s, c = sincos_point(x.lo)
    else:
        m = 0.5 * (x.lo + x.hi)
        r = up(max(m - x.lo, x.hi - m))
        if r > 0.5:
            return Iv(-1.0, 1.0), Iv(-1.0, 1.0)
        s0, c0 = sincos_point(m)
        sd = Iv(-r, r)
        cd = Iv(dn(1.0 - up(0.5 * up(r * r))), 1.0)
        s = s0 * cd + c0 * sd
        c = c0 * cd - s0 * sd
    return Iv(max(s.lo, -1.0), min(s.hi, 1.0)), Iv(max(c.lo, -1.0), min(c.hi, 1.0))


class Jet:
    """truncated Taylor series sum_k c[k] e^k with interval coefficients (c[k] encloses f^(k)/k! over the base)"""
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = c

    @staticmethod
    def const(x, K):
        return Jet([x if isinstance(x, Iv) else ivq(x)] + [ZERO] * K)

    @staticmethod
    def var(base, K):
        c = [base, ONE] + [ZERO] * (K - 1)
        return Jet(c[:K + 1])

    def K(self):
        return len(self.c) - 1

    def __add__(self, o):
        if isinstance(o, Jet):
            return Jet([a + b for a, b in zip(self.c, o.c)])
        o = o if isinstance(o, Iv) else ivq(o)
        return Jet([self.c[0] + o] + self.c[1:])

    __radd__ = __add__

    def __neg__(self):
        return Jet([-a for a in self.c])

    def __sub__(self, o):
        if isinstance(o, Jet):
            return Jet([a - b for a, b in zip(self.c, o.c)])
        o = o if isinstance(o, Iv) else ivq(o)
        return Jet([self.c[0] - o] + self.c[1:])

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if isinstance(o, Jet):
            a, b, K = self.c, o.c, len(self.c) - 1
            out = []
            for k in range(K + 1):
                s = a[0] * b[k]
                for i in range(1, k + 1):
                    s = s + a[i] * b[k - i]
                out.append(s)
            return Jet(out)
        o = o if isinstance(o, Iv) else ivq(o)
        return Jet([a * o for a in self.c])

    __rmul__ = __mul__

    def inv(self):
        a, K = self.c, len(self.c) - 1
        b0 = a[0].inv()
        b = [b0]
        for k in range(1, K + 1):
            s = a[1] * b[k - 1]
            for j in range(2, k + 1):
                s = s + a[j] * b[k - j]
            b.append(-(s * b0))
        return Jet(b)

    def __truediv__(self, o):
        if isinstance(o, Jet):
            return self * o.inv()
        o = o if isinstance(o, Iv) else ivq(o)
        return Jet([a / o for a in self.c])

    def sqrt(self):
        a, K = self.c, len(self.c) - 1
        b0 = a[0].sqrt()
        inv2b0 = (b0 * 2.0).inv()
        b = [b0]
        for k in range(1, K + 1):
            s = a[k]
            for j in range(1, k):
                s = s - b[j] * b[k - j]
            b.append(s * inv2b0)
        return Jet(b)

    def sincos(self):
        a, K = self.c, len(self.c) - 1
        s0, c0 = sincos(a[0])
        s, c = [s0], [c0]
        for k in range(1, K + 1):
            ss, cc = ZERO, ZERO
            for j in range(1, k + 1):
                ja = a[j] * float(j)
                ss = ss + ja * c[k - j]
                cc = cc + ja * s[k - j]
            s.append(ss / float(k))
            c.append(-(cc / float(k)))
        return Jet(s), Jet(c)

    def deriv(self, k):
        """enclosure of the k-th derivative over the base interval"""
        f = 1
        for i in range(2, k + 1):
            f *= i
        return self.c[k] * float(f)


def poly_jet(coeffs, x):
    """sum_e coeffs[e] x^e for exact rational coefficients (dict e -> Fraction), by Horner on jets"""
    deg = max(coeffs)
    acc = Jet.const(coeffs.get(deg, 0), x.K())
    for e in range(deg - 1, -1, -1):
        acc = acc * x
        ce = coeffs.get(e, 0)
        if ce:
            acc = acc + ivq(ce)
    return acc
