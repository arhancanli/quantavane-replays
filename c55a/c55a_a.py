#!/usr/bin/env python3
"""55a, program A: the computational inputs of [KLS2026] chi'_a(G) <= ceil(3.142 (Delta - 1)) + 1 (Kirousis, Livieratos,
Singh, arXiv 2602.14859v1), from the paper's text alone.

The bound rests on two numbers. Lemma 8 (boundtrees): t_n, the number of rooted unordered trees with n nodes whose
internal nodes all have even outdegree >= 4, has limsup t_n^(1/n) = 1/rho with rho >= 0.6677 (Theorem in the
Appendix: 0.6677 <= rho <= 0.6678). And with gamma = 1.142: (1 - e^(-1/gamma))^(3/4) < 0.6677, so that
(2 + gamma)(Delta - 1) + 1 = 3.142 (Delta - 1) + 1 colours suffice.

This program (route A):
  * t_n for n <= 120 from the multiset table: p[m][k] = number of multisets of k trees with m nodes in all, built
    one tree size at a time; t_(n) = sum over k in {0, 4, 6, 8, ...} of p[n-1][k].
  * T(x) = sum t_n x^n satisfies T = Phi(x; T(x), T(x^2), T(x^3), ...) with
    Phi(x; w) = x [ (exp(sum_j w_j/j) + exp(sum_j (-1)^j w_j/j)) / 2 - (w_1^2 + w_2)/2 ],
    and Phi has nonnegative coefficients in every w_j (it is x times a sum of cycle indices of symmetric groups).
    Hence if Phi(x; Y, W_2, W_3, ...) <= Y with W_j >= T(x^j), the monotone iteration from 0 stays below Y and
    T(x) <= Y < infinity, so rho >= x ("supersolution").
  * A chain of supersolutions: v = 3/10 with W_j = v^(j-1) Y (T(v)/v is increasing, so T(v^j) <= v^(j-1) T(v));
    then v = 1/2 and v = 3/5 with W_j = L(v^j) + Y_prev sum_(n>200) (v^j / v_prev)^n (Cauchy: t_n <= Y_prev/v_prev^n).
  * Level 0, x = 0.6677: W_j = L(x^j) + Y1 sum_(n>200) (x^j / (3/5))^n for 2 <= j <= 40 (L = the series to n = 200),
    W_j <= x^(j-2) W_2 beyond; the supersolution Y0 is found in floating point and checked in exact arithmetic.
  * For information, an upper bound: with L in place of T (coefficientwise smaller), Phi(x1; y, L(x1^2), ...) - y
    is convex in y and shown positive for all y >= 0 at x1 = 0.66771, so T has radius <= 0.66771.
  * gamma: (1 - e^(-1/1.142))^3 < 0.6677^4 exactly, and the least such gamma for 0.6677, for information.
Arithmetic: exact rationals; exp by its Taylor series with a remainder bound; every interval end rounded outward to
a multiple of 2^-300. Standard library only. Usage: c55a_a.py [--x0 X] [--gamma G] [--claim-rho R]
"""
import sys
from fractions import Fraction as Q
from math import comb, floor, exp as fexp

PREC = 300
def down(q):
    return Q(floor(q * (1 << PREC)), 1 << PREC)
def up(q):
    return -down(-q)

class I:
    """closed interval [lo, hi] of rationals, ends rounded outward"""
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = Q(lo if hi is None else hi)
        assert self.lo <= self.hi
    def __add__(s, o):
        o = o if isinstance(o, I) else I(o)
        return I(down(s.lo + o.lo), up(s.hi + o.hi))
    __radd__ = __add__
    def __sub__(s, o):
        o = o if isinstance(o, I) else I(o)
        return I(down(s.lo - o.hi), up(s.hi - o.lo))
    def __rsub__(s, o):
        return I(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, I) else I(o)
        c = [s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi]
        return I(down(min(c)), up(max(c)))
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = o if isinstance(o, I) else I(o)
        assert o.lo > 0 or o.hi < 0
        return s * I(down(1 / o.hi), up(1 / o.lo)) if o.lo > 0 else s * I(down(1 / o.hi), up(1 / o.lo))

def iexp(x):
    """e^x for an interval x, |x| <= 40"""
    def exp_pos(q, upper):           # q >= 0 rational: Taylor sum and remainder
        k, term, s = 0, Q(1), Q(1)
        K = int(3 * q) + 60
        for k in range(1, K + 1):
            term = up(term * q / k) if upper else down(term * q / k)
            s += term
        rem = term * q / (K + 1) / (1 - q / (K + 2))
        return up(s + rem) if upper else down(s)
    def e_at(q, upper):
        if q >= 0:
            return exp_pos(q, upper)
        v = exp_pos(-q, not upper)
        return up(1 / v) if upper else down(1 / v)
    return I(e_at(x.lo, False), e_at(x.hi, True))

# ------------------------------------------------------------------ the tree counts
N = 200
def counts(N):
    p = [[0] * (N + 1) for _ in range(N + 1)]
    p[0][0] = 1
    t = [0] * (N + 1)
    for n in range(1, N + 1):
        t[n] = sum(p[n - 1][k] for k in range(0, n) if k % 2 == 0 and (k != 2 or not DROP2))
        if t[n] == 0:
            continue
        new = [row[:] for row in p]
        for m in range(n, N + 1):
            for k in range(1, N + 1):
                s, i = 0, 1
                while i * n <= m and i <= k:
                    s += comb(t[n] + i - 1, i) * p[m - n * i][k - i]
                    i += 1
                new[m][k] += s
        p = new
    return t

args = sys.argv[1:]
X0 = Q("0.6677"); GAMMA = Q("1.142"); CLAIM = None; DROP2 = True
while args:
    a = args.pop(0)
    if a == "--x0":
        X0 = Q(args.pop(0))
    elif a == "--gamma":
        GAMMA = Q(args.pop(0))
    elif a == "--claim-rho":
        CLAIM = Q(args.pop(0))
    elif a == "--variant" and args and args[0] == "with-degree-2":
        args.pop(0); DROP2 = False; print("CONTROL variant: with-degree-2")
    else:
        sys.exit("bad argument")

print("program A (c55a_a.py): multiset table, supersolutions in exact interval arithmetic")
t = counts(N)
print("t_n, n = 1..41:", t[1:42])
print("t_99 =", t[99])

def L(xI):                      # sum_(n <= N) t_n x^n, x an interval >= 0
    s, pw = I(0), I(1)
    for n in range(1, N + 1):
        pw = pw * xI
        if t[n]:
            s = s + t[n] * pw
    return s

def phi(x, Y, W):               # W[j] for j = 2..J (intervals), and a tail bound for j > J (two sums)
    Wl, tailp, tailm = W
    Sp = Y + tailp
    Sm = I(0) - Y + tailm
    for j, wj in Wl.items():
        Sp = Sp + wj / j
        Sm = Sm + (wj / j if j % 2 == 0 else I(0) - wj / j)
    return x * ((iexp(Sp) + iexp(Sm)) / 2 - ((Y * Y + Wl[2]) / 2 if DROP2 else 0))

ok = True
# a chain of supersolutions: v = 3/10 with the crude bound T(v^j) <= v^(j-1) T(v); then v = 1/2 and v = 3/5, where
# T(v^j) <= L(v^j) + Y_prev sum_(n > N) (v^j / v_prev)^n (Cauchy: t_n <= Y_prev / v_prev^n) for 2 <= j <= 40 and
# T(v^j) <= v^(j-2) T(v^2) beyond
def level_bounds(vq, prev):
    """W_j (j >= 2) and the tail bound for the sum over j beyond them, at the point v = vq"""
    v = I(vq)
    Wl = {}
    if prev is None:
        return None
    pv, pY = prev
    for j in range(2, 41):
        q = vq ** j / pv
        assert q < 1
        Wl[j] = L(I(vq ** j)) + I(0, up(pY * q ** (N + 1) / (1 - q)))
    tl = Wl[2] * I(vq ** 39) / I(41) / (I(1) - v)
    return (Wl, I(0, tl.hi))
def level_phi(vq, bounds, Y):
    v = I(vq)
    if bounds is None:                                  # crude: W_j = v^(j-1) Y
        Wl, vp = {}, I(1)
        for j in range(2, 61):
            vp = vp * v
            Wl[j] = vp * Y
        tl = Y * I(vq ** 60) / I(61) / (I(1) - v)
        return phi(v, Y, (Wl, I(0, tl.hi), I(0, tl.hi)))
    Wl, tl = bounds
    return phi(v, Y, (Wl, tl, tl))
def float_min(fun, lo_, hi_):
    for _ in range(100):
        m1, m2 = lo_ + (hi_ - lo_) / 3, hi_ - (hi_ - lo_) / 3
        if fun(m1) < fun(m2):
            hi_ = m2
        else:
            lo_ = m1
    return (lo_ + hi_) / 2
chain = []
prev = None
for vq in (Q(3, 10), Q(1, 2), Q(3, 5)):
    bounds = level_bounds(vq, prev)
    # choose Y where Phi(v; Y) - Y is least (floating point on the interval upper ends), then test it exactly
    gap = lambda y: float(level_phi(vq, bounds, I(Q(round(y * 10 ** 6), 10 ** 6))).hi) - y
    yf = float_min(gap, float(vq), 3.0)
    if gap(yf) < 0:                    # Phi(Y) - Y is convex: take the left end of {gap <= 0}, plus 1/1000
        lo_, hi_ = float(vq), yf
        for _ in range(40):
            mid = (lo_ + hi_) / 2
            if gap(mid) <= 0:
                hi_ = mid
            else:
                lo_ = mid
        yf = hi_
    Yq = Q(round(yf * 1000) + 1, 1000)
    P = level_phi(vq, bounds, I(Yq))
    if P.hi > Yq:
        print("no supersolution at v = %s (least float gap at Y = %.6f)" % (vq, yf))
        break
    print("level v = %s: Phi(v; Y, W) <= Y at Y = %s (upper end %.12f), so T(v) <= %s" % (vq, Yq, float(P.hi), Yq))
    chain.append((vq, Yq))
    prev = (vq, Yq)
chain_ok = len(chain) == 3
ok &= chain_ok
V1, Y1 = chain[-1] if chain_ok else (Q(3, 10), Q(10))
# level 0 at X0
x = I(X0)
J = 40
Wl = {}
for j in range(2, J + 1):
    xj = I(X0 ** j)
    q = X0 ** j / V1                     # ratio for the Cauchy tail (t_n <= Y1 / V1^n)
    q = q if q < 1 else Q(99, 100)       # (only when the chain failed; the verdict is then REFUSED anyway)
    tail = I(0, up(Y1 * q ** (N + 1) / (1 - q)))
    Wl[j] = L(xj) + tail
# j > J: W_j <= x^(j-2) W_2, so sum_(j > J) W_j / j <= W_2 x^(J-1) / ((J+1)(1 - x))
tl = Wl[2] * I(X0 ** (J - 1)) / I(J + 1) / (I(1) - x)
tl = I(0, tl.hi)
# find Y0 numerically (floating point), then check exactly
def fphi(y):
    Sp = y + sum(float(w.hi) / j for j, w in Wl.items()) + float(tl.hi)
    Sm = -y + sum((float(w.hi) if j % 2 == 0 else -float(w.hi)) / j for j, w in Wl.items()) + float(tl.hi)
    return float(X0) * ((fexp(Sp) + fexp(Sm)) / 2 - (y * y + float(Wl[2].hi)) / 2) - y
lo_, hi_ = 0.5, 3.0
for _ in range(200):                       # minimise the convex function fphi on [0.5, 3]
    m1, m2 = lo_ + (hi_ - lo_) / 3, hi_ - (hi_ - lo_) / 3
    if fphi(m1) < fphi(m2):
        hi_ = m2
    else:
        lo_ = m1
Y0 = Q(round((lo_ + hi_) / 2 * 10 ** 12), 10 ** 12)
P0 = phi(x, I(Y0), (Wl, tl, tl))
margin = Y0 - P0.hi
print("level 0: x = %s, Y0 = %s, Phi(x; Y0, W_2, ...) <= %.15f, margin %.3e" % (X0, Y0, float(P0.hi), float(margin)))
super_ok = P0.hi <= Y0
print("supersolution at x = %s: %s, so rho >= %s" % (X0, super_ok, X0 if super_ok else "?"))
ok &= super_ok

# information: upper bound with the truncated system at x1
X1 = Q("0.66771")
x1 = I(X1)
W1 = {j: L(I(X1 ** j)) for j in range(2, N + 1)}        # lower values (truncations), no tails
zero = I(0)
def g_and_dg(yq):
    # f(y) = Phi_trunc(x1; y) - y and f'(y), as intervals, at a rational y
    y = I(yq)
    Sp = y; Sm = zero - y
    for j, w in W1.items():
        Sp = Sp + w / j
        Sm = Sm + (w / j if j % 2 == 0 else zero - w / j)
    ep, em = iexp(Sp), iexp(Sm)
    f = x1 * ((ep + em) / 2 - ((y * y + W1[2]) / 2 if DROP2 else 0)) - y
    df = x1 * ((ep - em) / 2 - (y if DROP2 else 0)) - 1
    return f, df
# convex in y (second derivative x1 (ep + em)/2 - x1 >= 0 since Sp + Sm >= 0); tangents at ya < yb bound it below
ya, yb = Q("1.3420"), Q("1.3450")
fa, da = g_and_dg(ya)
fb, db = g_and_dg(yb)
lower_mid = None
if da.hi < 0 < db.lo:
    # f convex: on [ya, yb], f(y) >= A(y) = fa.lo + da.lo (y - ya) and f(y) >= B(y) = fb.lo + db.hi (y - yb)
    # (tangent lines, with the slope ends that make them lower); max(A, B) >= its value where A = B.
    # Below ya f decreases to f(ya) and above yb it increases from f(yb), both >= their lower ends.
    ys = (fb.lo - fa.lo + da.lo * ya - db.hi * yb) / (da.lo - db.hi)
    lower_mid = fa.lo + da.lo * (ys - ya)
upper_ok = lower_mid is not None and lower_mid > 0 and fa.lo > 0 and fb.lo > 0
print("information: truncated system at x1 = %s: f(y) > 0 for all y >= 0 (tangent bound %.3e): %s, so rho <= %s"
      % (X1, float(lower_mid) if lower_mid is not None else float("nan"), upper_ok, X1))

# gamma
g = GAMMA
e = iexp(I(Q(-1) / g))
lhs = (I(1) - e) * (I(1) - e) * (I(1) - e)
rhs = I(Q("0.6677") ** 4)
gam_ok = lhs.hi < rhs.lo
print("gamma = %s: (1 - e^(-1/gamma))^3 in [%.15f, %.15f] < 0.6677^4 = %.15f: %s"
      % (g, float(lhs.lo), float(lhs.hi), float(rhs.lo), gam_ok))
# least gamma (information): bisection on (1 - e^(-1/g))^3 < 0.6677^4
glo, ghi = Q(1), Q(2)
for _ in range(40):
    gm = (glo + ghi) / 2
    em = iexp(I(Q(-1) / gm)); lm = (I(1) - em) * (I(1) - em) * (I(1) - em)
    if lm.hi < rhs.lo:
        ghi = gm
    else:
        glo = gm
print("information: the least such gamma lies in [%.9f, %.9f]" % (float(glo), float(ghi)))
ok &= gam_ok
print("colours: (2 + gamma)(Delta - 1) + 1 = %s (Delta - 1) + 1" % (2 + g))
ok &= (2 + g == Q("3.142")) if GAMMA == Q("1.142") else True
# [FLM2020] 3.569 and [GKPT2017] 3.74: rho_gamma = min_(x>0) phi(x)/x < 1, with a = 1 - e^(-1/gamma) and
# phi(x) = (1/gamma) a^3 (x+1)^p / (1 - a^2 (x+1)^2), p = 4 (FLM2020, gamma = 1.569) or p = 6 (GKPT2017, gamma = 1.74).
# Any x in (0, 1/a - 1) with phi(x)/x < 1 shows rho_gamma < 1; x is taken from a floating-point search.
def rho_row(gq, p):
    aI = I(1) - iexp(I(Q(-1) / gq))
    def ratio_hi(xq):
        xI = I(xq)
        u = xI + 1
        num = aI * aI * aI
        for _ in range(p):
            num = num * u
        den = (I(1) - aI * aI * u * u) * xI * gq
        assert den.lo > 0
        return (num / den).hi
    af = 1 - fexp(-1 / float(gq))
    f = lambda x: (af ** 3 * (x + 1) ** p / (1 - af * af * (x + 1) ** 2)) / float(gq) / x
    lo_, hi_ = 1e-6, 1 / af - 1 - 1e-9
    for _ in range(200):
        m1, m2 = lo_ + (hi_ - lo_) / 3, hi_ - (hi_ - lo_) / 3
        if f(m1) < f(m2):
            hi_ = m2
        else:
            lo_ = m1
    xs = Q(round((lo_ + hi_) / 2 * 10 ** 9), 10 ** 9)
    return xs, ratio_hi(xs)
for name, gq, p in (("FLM2020", Q("1.569"), 4), ("GKPT2017", Q("1.74"), 6)):
    if GAMMA != Q("1.142"):
        break
    xs, rh = rho_row(gq, p)
    good = rh < 1
    print("[%s] gamma = %s, p = %d: phi(x)/x <= %.12f at x = %s: rho_gamma < 1: %s; colours (2 + gamma)(Delta - 1) = %s (Delta - 1)"
          % (name, gq, p, float(rh), xs, good, 2 + gq))
    ok &= good
if CLAIM is not None:
    c = CLAIM <= X0 and super_ok
    print("CONTROL claim rho >=", CLAIM, ":", c)
    ok &= c
print("VERDICT", "PASS" if ok else "REFUSED")
