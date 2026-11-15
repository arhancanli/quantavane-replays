#!/usr/bin/env python3
"""55a, program B (written separately from program A): the computational inputs of [KLS2026]
chi'_a(G) <= ceil(3.142 (Delta - 1)) + 1 (arXiv 2602.14859v1): rho >= 0.6677 for the generating function T(z) of rooted
unordered trees whose internal nodes have even outdegree >= 4, and (1 - e^(-1/1.142))^(3/4) < 0.6677.

Route B:
  * t_n for n <= 100 from the functional equation of the paper's Lemma 9 itself, coefficient by coefficient:
    T = z [ (exp(S+) + exp(S-)) / 2 - (T(z)^2 + T(z^2)) / 2 ], S+ = sum_j T(z^j)/j, S- = sum_j (-1)^j T(z^j)/j, the
    exponentials by f_m = (1/m) sum_k k s_k f_(m-k), in exact rationals (every t_n must come out an integer).
  * Decimal interval arithmetic (60 digits, directed rounding; exp and ln of the decimal module, correctly rounded,
    widened by one unit in the last place).
  * Supersolution test (Phi has nonnegative coefficients in every argument, so Phi(x; Y, W_2, ...) <= Y with
    W_j >= T(x^j) gives T(x) <= Y): along v = 0.25 (crude bound T(v^j) <= v^(j-1) Y), 0.45 and 0.62 (Cauchy tails
    from the previous point), giving T(0.62) <= Y1 and t_n <= Y1 / 0.62^n; then at x = 0.6677 with W_j = L(x^j) + sum_(n>100) Y1 (x^j / 0.62)^n for j <= 30 and, for
    j > 30, T(x^j) <= x^j Y1 / 0.62.
  * gamma by logarithms: (3/4) ln(1 - e^(-1/gamma)) < ln 0.6677.
Usage: c55a_b.py [--x0 X] [--gamma G] [--variant with-degree-2]
"""
import sys
from fractions import Fraction
from decimal import Decimal as D, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

getcontext().prec = 60

def rd(f):                     # lower decimal of a Fraction
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        return D(f.numerator) / D(f.denominator)
def ru(f):
    with localcontext() as c:
        c.rounding = ROUND_CEILING
        return D(f.numerator) / D(f.denominator)

def op(a, b, f, rnd):
    with localcontext() as c:
        c.rounding = rnd
        return f(a, b)

class Iv:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = lo, (lo if hi is None else hi)
        assert self.lo <= self.hi, (lo, hi)
    def __add__(a, b):
        b = b if isinstance(b, Iv) else Iv(b)
        return Iv(op(a.lo, b.lo, lambda x, y: x + y, ROUND_FLOOR), op(a.hi, b.hi, lambda x, y: x + y, ROUND_CEILING))
    def __sub__(a, b):
        b = b if isinstance(b, Iv) else Iv(b)
        return Iv(op(a.lo, b.hi, lambda x, y: x - y, ROUND_FLOOR), op(a.hi, b.lo, lambda x, y: x - y, ROUND_CEILING))
    def __mul__(a, b):
        b = b if isinstance(b, Iv) else Iv(b)
        lo = min(op(p, q, lambda x, y: x * y, ROUND_FLOOR) for p in (a.lo, a.hi) for q in (b.lo, b.hi))
        hi = max(op(p, q, lambda x, y: x * y, ROUND_CEILING) for p in (a.lo, a.hi) for q in (b.lo, b.hi))
        return Iv(lo, hi)
    def scale(a, num, den):     # a * num / den for positive integers
        return Iv(op(op(a.lo, D(num), lambda x, y: x * y, ROUND_FLOOR), D(den), lambda x, y: x / y, ROUND_FLOOR),
                  op(op(a.hi, D(num), lambda x, y: x * y, ROUND_CEILING), D(den), lambda x, y: x / y, ROUND_CEILING))

def widen(x, up):
    return x.next_plus() if up else x.next_minus()

def ex(a):
    with localcontext() as c:
        c.rounding = ROUND_HALF_EVEN
        lo, hi = a.lo.exp(), a.hi.exp()
    return Iv(widen(lo, False), widen(hi, True))

def ln(a):
    with localcontext() as c:
        c.rounding = ROUND_HALF_EVEN
        lo, hi = a.lo.ln(), a.hi.ln()
    return Iv(widen(lo, False), widen(hi, True))

# ------------------------------------------------------------------ arguments
x0 = D("0.6677"); gam = D("1.142"); variant = None
argv = sys.argv[1:]
while argv:
    k = argv.pop(0)
    if k == "--x0":
        x0 = D(argv.pop(0))
    elif k == "--gamma":
        gam = D(argv.pop(0))
    elif k == "--variant":
        variant = argv.pop(0)
    else:
        raise SystemExit("bad argument")
DROP2 = (variant != "with-degree-2")        # the paper's class removes outdegree 2

# ------------------------------------------------------------------ coefficients from the functional equation
NB = 100
def coefficients(NB):
    t = [Fraction(0)] * (NB + 1)
    sp = [Fraction(0)] * (NB + 1)          # S+ coefficients
    sm = [Fraction(0)] * (NB + 1)          # S- coefficients
    fp = [Fraction(1)] + [Fraction(0)] * NB
    fm = [Fraction(1)] + [Fraction(0)] * NB
    for n in range(1, NB + 1):
        m = n - 1
        if m >= 1:
            fp[m] = sum(k * sp[k] * fp[m - k] for k in range(1, m + 1)) / m
            fm[m] = sum(k * sm[k] * fm[m - k] for k in range(1, m + 1)) / m
        sq = sum((t[i] * t[m - i] for i in range(1, m)), Fraction(0))       # [z^m] T^2
        t2 = t[m // 2] if m % 2 == 0 else Fraction(0)                      # [z^m] T(z^2)
        val = (fp[m] + fm[m]) / 2 - ((sq + t2) / 2 if DROP2 else 0)
        assert val.denominator == 1
        t[n] = val
        for j in range(1, NB // n + 1):                                    # T(z^j)/j gains t_n / j at z^(n j)
            sp[n * j] += t[n] / j
            sm[n * j] += (t[n] / j) * (1 if j % 2 == 0 else -1)
    return [int(c) for c in t]

print("program B (c55a_b.py): coefficients from the functional equation; decimal intervals")
if variant:
    print("CONTROL variant:", variant)
t = coefficients(NB)
print("t_n, n = 1..41:", t[1:42])
print("t_99 =", t[99])

def Lsum(x):                                   # x an Iv, sum_(n<=NB) t_n x^n
    s, p = Iv(D(0)), Iv(D(1))
    for n in range(1, NB + 1):
        p = p * x
        if t[n]:
            s = s + p * Iv(D(t[n]))
    return s

def Phi(x, Y, W, tail_all):
    """x Iv, Y Iv, W: dict j -> Iv for 2 <= j <= J, tail_all: upper bound of sum_(j > J) T(x^j)/j"""
    Sp = Y + Iv(D(0), tail_all)
    Sm = Iv(D(0)) - Y + Iv(D(0), tail_all)
    for j, w in W.items():
        Sp = Sp + w.scale(1, j)
        Sm = (Sm + w.scale(1, j)) if j % 2 == 0 else (Sm - w.scale(1, j))
    half = (ex(Sp) + ex(Sm)).scale(1, 2)
    if DROP2:
        half = half - (Y * Y + W[2]).scale(1, 2)
    return x * half

ok = True
# supersolutions along v = 0.25 (crude bound T(v^j) <= v^(j-1) T(v)), v = 0.45, v = 0.62 (Cauchy tails from the
# previous point for j <= 30, and T(v^j) <= v^j T(v_prev) / v_prev for j > 30)
def level(vs, prev, Y):
    v = Iv(D(vs))
    W = {}
    if prev is None:
        p = Iv(D(1))
        for j in range(2, 81):
            p = p * v
            W[j] = p * Y
        with localcontext() as c:
            c.rounding = ROUND_CEILING
            tail = Y.hi * D(vs) ** 80 / 81 / (1 - D(vs))
        return Phi(v, Y, W, tail)
    pv, pY = prev
    for j in range(2, 31):
        with localcontext() as c:
            c.rounding = ROUND_CEILING
            vj_hi = D(vs) ** j
            q = vj_hi / D(pv)
            tail = pY * q ** (NB + 1) / (1 - q)
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            vj_lo = D(vs) ** j
        W[j] = Lsum(Iv(vj_lo, vj_hi)) + Iv(D(0), tail)
    with localcontext() as c:
        c.rounding = ROUND_CEILING
        tail_all = (pY / D(pv)) * D(vs) ** 31 / 31 / (1 - D(vs))
    return Phi(v, Y, W, tail_all)
prev = None
for vs in ("0.25", "0.45", "0.62"):
    Yok = None
    k = int(D(vs) * 500)
    while k <= 1500:
        y = D(k) / 500
        r = level(vs, prev, Iv(y))
        if r.hi <= y:
            Yok = y
            break
        k += 1
    if Yok is None:
        print("no supersolution at v =", vs)
        ok = False
        break
    print("level v = %s: Phi(v; Y, W) <= Y at Y = %s (upper end %s), so T(v) <= %s" % (vs, Yok, format(r.hi, ".12f"), Yok))
    prev = (vs, Yok)
chain_ok = prev is not None and prev[0] == "0.62"
ok = ok and chain_ok
Y1 = Iv(prev[1] if chain_ok else D(10))
# level 0
x = Iv(x0)
J = 30
W = {}
for j in range(2, J + 1):
    xj = Iv(op(x0, j, lambda a, b: a ** b, ROUND_FLOOR), op(x0, j, lambda a, b: a ** b, ROUND_CEILING))
    with localcontext() as c:
        c.rounding = ROUND_CEILING
        q = xj.hi / D("0.62")      # Cauchy: t_n <= Y1 / 0.62^n
        tail = Y1.hi * q ** (NB + 1) / (1 - q)
    W[j] = Lsum(xj) + Iv(D(0), tail)
with localcontext() as c:
    c.rounding = ROUND_CEILING
    tail_all = (Y1.hi / D("0.62")) * x0 ** (J + 1) / (J + 1) / (1 - x0)     # sum_(j>J) x^j Y1/0.62 / j
# Y0: ternary search in plain decimals, then the exact test
def g(y):
    r = Phi(x, Iv(y), W, tail_all)
    return r.hi - y
a, b = D("0.5"), D("3")
for _ in range(120):
    m1, m2 = a + (b - a) / 3, b - (b - a) / 3
    if g(m1) < g(m2):
        b = m2
    else:
        a = m1
Y0 = ((a + b) / 2).quantize(D("1e-12"))
r0 = Phi(x, Iv(Y0), W, tail_all)
sup_ok = r0.hi <= Y0
print("level 0: x = %s, Y0 = %s, Phi <= %s, margin %.3e" % (x0, Y0, format(r0.hi, ".15f"), float(Y0 - r0.hi)))
print("supersolution at x = %s: %s, so rho >= %s" % (x0, sup_ok, x0 if sup_ok else "?"))
ok = ok and sup_ok

# gamma by logarithms
with localcontext() as c:
    c.rounding = ROUND_FLOOR
    inv_lo = D(1) / gam
    c.rounding = ROUND_CEILING
    inv_hi = D(1) / gam
e = ex(Iv(-inv_hi, -inv_lo))
lhs = ln(Iv(D(1)) - e).scale(3, 4)
rhs = ln(Iv(D("0.6677")))
g_ok = lhs.hi < rhs.lo
print("gamma = %s: (3/4) ln(1 - e^(-1/gamma)) <= %s < ln 0.6677 >= %s: %s"
      % (gam, format(lhs.hi, ".15f"), format(rhs.lo, ".15f"), g_ok))
ok = ok and g_ok
print("colours: (2 + gamma)(Delta - 1) + 1 = %s (Delta - 1) + 1" % (2 + gam))
# [FLM2020] 3.569 and [GKPT2017] 3.74: min_(x>0) phi(x)/x < 1 is shown by one x where the decimal interval of
# phi(x)/x lies below 1; phi(x) = (1/gamma) a^3 (x+1)^p / (1 - a^2 (x+1)^2), a = 1 - e^(-1/gamma)
def row(gs, p):
    g = D(gs)
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        ilo = D(1) / g
        c.rounding = ROUND_CEILING
        ihi = D(1) / g
    a = Iv(D(1)) - ex(Iv(-ihi, -ilo))
    # golden-section search for x in plain decimals
    amid = a.lo
    lo_, hi_ = D("0.000001"), D(1) / amid - 1 - D("1e-12")
    def fv(x):
        u = x + 1
        return a.lo ** 3 * u ** p / (1 - a.lo ** 2 * u ** 2) / g / x
    for _ in range(150):
        m1, m2 = lo_ + (hi_ - lo_) / 3, hi_ - (hi_ - lo_) / 3
        if fv(m1) < fv(m2):
            hi_ = m2
        else:
            lo_ = m1
    x = ((lo_ + hi_) / 2).quantize(D("1e-9"))
    u = Iv(x + 1)
    num = a * a * a
    for _ in range(p):
        num = num * u
    den = (Iv(D(1)) - a * a * u * u) * Iv(x) * Iv(g)
    assert den.lo > 0
    with localcontext() as c:
        c.rounding = ROUND_CEILING
        hi = num.hi / den.lo
    return x, hi
if gam == D("1.142") and DROP2:
    for name, gs, p in (("FLM2020", "1.569", 4), ("GKPT2017", "1.74", 6)):
        x, hi = row(gs, p)
        good = hi < 1
        print("[%s] gamma = %s, p = %d: phi(x)/x <= %s at x = %s: rho_gamma < 1: %s" % (name, gs, p, format(hi, ".12f"), x, good))
        ok = ok and good
print("VERDICT", "PASS" if ok else "REFUSED")
