#!/usr/bin/env python3
"""33a, program B (written separately from program A): the finite inputs of [DM2013] A(2) >= 0.316999...
(Duursma and Mak, arXiv 1102.4127v4, Proposition 5.5 and Theorem 1.3 in its numbering; the criterion is its Theorem 1.1), from the paper's text alone.

Same objects as program A (H: y^2 + (x^3+x+1) y = x^2 + x over F_2; k = compositum of v^2 + (x^5+x^2+1) v =
(x^2+x) h, h in span{1, x, x^2, y, y^2}), different route:
  * F_(2^d) by logarithm tables for the numerically LARGEST primitive polynomial of degree d; the trace as a linear
    form (traces of the basis vectors); w^2 + w = c by Gaussian elimination over F_2.
  * Closed points of H of each degree e <= 10 (points over F_(2^e) in no proper subfield, divided by e); for an
    unramified closed point P the Frobenius of k/H is the vector (Tr_(F_(2^e)/F_2) f_h(P))_h; P has 32 places of
    degree e above it if that vector is 0, otherwise 16 places of degree 2e.  Places at infinity from valuations
    (v(x) = -1; v(y) = -3 or +1 by the Newton polygon; v(f_h) >= 2 > 0 there).
  * Ramification at the zeros of A with the local parameter A itself: x - x0 as a power series in A (series
    reversion), then y, then f_h as a Laurent series in A.
  * The genus independently, by the zeta functions: for each of the 31 characters the L-polynomial
    (degree 2(g(C) - g(H)) = 22 if the conductor has degree 20) is solved from the point counts over F_(2^d),
    d <= 11, with the functional equation, and must predict the counts for d = 12, 13, 14; g(H) = 2 the same way.
  * Theorem 1.1's conditions and the bound, in exact rationals.
Usage: c33a_b.py [--claim X] [--variant NAME] [--S "5:2,6:16,8:15,10:4"]
"""
import sys
from fractions import Fraction

def clmul(a, b):
    out = 0
    i = 0
    while (b >> i):
        if (b >> i) & 1:
            out ^= a << i
        i += 1
    return out

def reduce_mod(a, m):
    dm = m.bit_length()
    while a.bit_length() >= dm:
        a ^= m << (a.bit_length() - dm)
    return a

class Field:
    """F_(2^d) with log tables over the largest primitive polynomial of degree d."""
    def __init__(self, d):
        self.d, self.size = d, 1 << d
        for m in range((1 << (d + 1)) - 1, 1 << d, -2):
            exp, z = [], 1
            for _ in range(self.size - 1):
                exp.append(z)
                z = reduce_mod(z << 1, m)
            if z == 1 and len(set(exp)) == self.size - 1:
                self.modulus = m
                break
        self.exp = exp + exp
        self.log = [None] * self.size
        for i, z in enumerate(exp):
            self.log[z] = i
        # trace as a linear form: Tr(e_i) for the basis 1, X, ..., X^(d-1)
        self.trmask = 0
        for i in range(d):
            t, z = 0, 1 << i
            for _ in range(d):
                t ^= z
                z = self.m(z, z)
            assert t in (0, 1)
            self.trmask |= t << i
        # w -> w^2 + w as a d x d matrix over F_2 (columns = images of basis vectors)
        self.cols = [self.m(1 << i, 1 << i) ^ (1 << i) for i in range(d)]
    def m(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]
    def div(self, a, b):
        if a == 0:
            return 0
        return self.exp[self.log[a] - self.log[b] + self.size - 1]
    def tr(self, a):
        return bin(a & self.trmask).count("1") & 1
    def sq(self, a):
        return self.m(a, a)
    def root(self, a):            # square root
        if a == 0:
            return 0
        l = self.log[a]
        return self.exp[(l if l % 2 == 0 else l + self.size - 1) // 2]
    def solve_as(self, c):
        """one w with w^2 + w = c, or None; Gaussian elimination on [cols | c]."""
        d = self.d
        rows = []
        for r in range(d):
            bits = 0
            for j in range(d):
                bits |= ((self.cols[j] >> r) & 1) << j
            rows.append((bits, (c >> r) & 1))
        piv = []
        for col in range(d):
            sel = next((k for k in range(len(rows)) if (rows[k][0] >> col) & 1), None)
            if sel is None:
                continue
            pr = rows.pop(sel)
            rows = [((b ^ pr[0], v ^ pr[1]) if (b >> col) & 1 else (b, v)) for (b, v) in rows]
            piv = [((b ^ pr[0], v ^ pr[1]) if (b >> col) & 1 else (b, v)) for (b, v) in piv]
            piv.append(pr)
        if any(b == 0 and v for (b, v) in rows):
            return None
        w = 0
        for (b, v) in piv:
            col = (b & -b).bit_length() - 1
            if v:
                w |= 1 << col
        assert self.m(w, w) ^ w == c
        return w
    def poly(self, coeffs, x):     # coeffs: list of 0/1, low degree first
        r = 0
        for c in reversed(coeffs):
            r = self.m(r, x) ^ c
        return r

# --------------------------------------------------------------- data
a_poly = [1, 1, 0, 1]          # 1 + x + x^3
b_poly = [0, 1, 1]             # x + x^2
A_poly = [1, 0, 1, 0, 0, 1]    # 1 + x^2 + x^5
basis = ["1", "x", "xx", "y", "yy"]
paper_aH = [6, 0, 1, 1, 6]
paper_ak = [192, 0, 0, 0, 2, 16, 0, 16, 0, 64]
paper_gk = 343
S = {5: 2, 6: 16, 8: 15, 10: 4}
t_size = 192
claim, variant = None, None
argv = sys.argv[1:]
while argv:
    k = argv.pop(0)
    if k == "--claim":
        claim = Fraction(argv.pop(0))
    elif k == "--variant":
        variant = argv.pop(0)
    elif k == "--S":
        S = {}
        for item in argv.pop(0).split(","):
            f, c = item.split(":")
            S[int(f)] = int(c)
    else:
        raise SystemExit("bad argument " + k)
if variant == "h-x3":
    basis = ["1", "x", "xx", "y", "xxx"]
elif variant == "A-x5x3":
    A_poly = [1, 0, 0, 1, 0, 1]
elif variant:
    raise SystemExit("unknown variant")

def basis_value(F, nm, x, y):
    v = 1
    for ch in nm:
        if ch == "x":
            v = F.m(v, x)
        elif ch == "y":
            v = F.m(v, y)
    return v

# --------------------------------------------------------------- points
def affine_points(F):
    out = []
    for x in range(F.size):
        a, b = F.poly(a_poly, x), F.poly(b_poly, x)
        if a == 0:
            out.append((x, F.root(b)))
        else:
            w = F.solve_as(F.div(b, F.m(a, a)))
            if w is not None:
                out.append((x, F.m(a, w)))
                out.append((x, F.m(a, w ^ 1)))
    return out

def infinity():
    """valuations at the places over x = infinity: v(a) = -3, v(b) = -2, so the Newton polygon of
    Y^2 + aY - b gives two places (v(y) = -3 and v(y) = v(b) - v(a) = +1), each with v(x) = -1, rational
    (residue equations w^2 + w = 0 and w = b/a leading coefficient)."""
    va, vb = -(len(a_poly) - 1), -(len(b_poly) - 1)
    places_inf = []
    for vy in (va, vb - va):
        vals = {"1": 0}
        vh = min(nm.count("x") * -1 + nm.count("y") * vy for nm in basis)   # lower bound for v(h)
        vf = (-2) + vh - 2 * (-(len(A_poly) - 1))                             # v(x^2+x) = -2, v(A) = -5
        places_inf.append(vf)
    return places_inf

def frob_vector(F, x, y):
    A = F.poly(A_poly, x)
    base = F.div(F.poly(b_poly, x), F.m(A, A))
    v = 0
    for j, nm in enumerate(basis):
        v |= F.tr(F.m(base, basis_value(F, nm, x, y))) << j
    return v

fields = {}
def field(d):
    if d not in fields:
        fields[d] = Field(d)
    return fields[d]

def in_subfield(F, z, e):
    """z^(2^e) == z"""
    w = z
    for _ in range(e):
        w = F.sq(w)
    return w == z

# --------------------------------------------------------------- closed points and decomposition
print("program B (c33a_b.py): F_(2^d) by log tables over the largest primitive polynomial; closed points")
if variant:
    print("CONTROL variant:", variant)
inf_vals = infinity()
print("places of H at infinity: 2 rational, v(f_h) >= %s there (> 0, so the Frobenius is trivial)" % inf_vals)
inf_ok = all(v > 0 for v in inf_vals)
DMAX = 10
closed_H = [0] * (DMAX + 1)
ak = [0] * (2 * DMAX + 1)
ram_points = []
for e in range(1, DMAX + 1):
    F = field(e)
    pts = affine_points(F)
    exact = [(x, y) for (x, y) in pts if not any(e % s == 0 and in_subfield(F, x, s) and in_subfield(F, y, s)
                                                  for s in range(1, e))]
    assert len(exact) % e == 0
    seen = set()
    for (x, y) in exact:
        if (x, y) in seen:
            continue
        orbit = [(x, y)]
        xx, yy = F.sq(x), F.sq(y)
        while (xx, yy) != (x, y):
            orbit.append((xx, yy))
            xx, yy = F.sq(xx), F.sq(yy)
        assert len(orbit) == e
        seen.update(orbit)
        closed_H[e] += 1
        if F.poly(A_poly, x) == 0:
            ram_points.append((e, x, y))
            continue
        if frob_vector(F, x, y) == 0:
            ak[e] += 32
        else:
            ak[2 * e] += 16
    if e == 1:
        closed_H[1] += 2
        ak[1] += 2 * 32 if inf_ok else 0
print("closed points of H by degree, e = 1..10:", closed_H[1:])
print("closed points of H on A = 0:", [(e) for (e, _, _) in ram_points])

# --------------------------------------------------------------- ramification with the local parameter A
def ramification(e, x0, y0, nterms=6):
    F = field(e)
    def smul(s, t):
        r = [0] * nterms
        for i in range(nterms):
            if s[i]:
                for j in range(nterms - i):
                    r[i + j] ^= F.m(s[i], t[j])
        return r
    def scompose_poly(coeffs, s):   # coeffs(x0 + s(A)) as a series in A, s without constant term
        r = [0] * nterms
        pw = [1] + [0] * (nterms - 1)
        xs = [x0] + s[1:]
        for c in coeffs:
            if c:
                r = [u ^ v for u, v in zip(r, pw)]
            pw = smul(pw, xs)
        return r
    # x = x0 + s(A): solve A(x0 + s) = A by fixed-point iteration on the coefficients
    s = [0] * nterms
    dA = F.poly([A_poly[i] * i % 2 for i in range(1, len(A_poly))], x0)   # A'(x0)
    assert dA != 0
    for _ in range(nterms + 1):
        val = scompose_poly(A_poly, s)          # should equal [0, 1, 0, 0, ...]
        target = [0, 1] + [0] * (nterms - 2)
        err = [u ^ v for u, v in zip(val, target)]
        s = [u ^ F.div(v, dA) for u, v in zip(s, err)]
    assert scompose_poly(A_poly, s) == [0, 1] + [0] * (nterms - 2)
    xser = [x0] + s[1:]
    aser, bser = scompose_poly(a_poly, s), scompose_poly(b_poly, s)
    # y as a series: iterate y <- (b + y^2) / a
    ainv0 = F.div(1, aser[0])
    yser = [y0] + [0] * (nterms - 1)
    for _ in range(nterms + 1):
        rhs = [u ^ v for u, v in zip(bser, smul(yser, yser))]
        # solve a * y = rhs as series
        ynew = [0] * nterms
        for k2 in range(nterms):
            acc = rhs[k2]
            for i in range(1, k2 + 1):
                acc ^= F.m(aser[i], ynew[k2 - i])
            ynew[k2] = F.m(acc, ainv0)
        yser = ynew
    chk = [u ^ v ^ w for u, v, w in zip(smul(yser, yser), smul(aser, yser), bser)]
    assert not any(chk)
    out = 0
    for h in range(1, 32):
        hser = [0] * nterms
        for j, nm in enumerate(basis):
            if (h >> j) & 1:
                term = [1] + [0] * (nterms - 1)
                for ch in nm:
                    if ch in "xy":
                        term = smul(term, xser if ch == "x" else yser)
                hser = [u ^ v for u, v in zip(hser, term)]
        num = smul(bser, hser)                  # f_h = num / A^2: Laurent coefficient of A^-2 is num[0], of A^-1 num[1]
        if num[1] ^ F.root(num[0]):
            out += 1
    return out

ram_total = 0
for (e, x0, y0) in ram_points:
    r = ramification(e, x0, y0)
    ram_total += r
    if r == 31:
        ak[e] += 1          # every character ramified: the place is totally ramified in k
cond_ok = (len(ram_points) == 2 and all(e == 5 for (e, _, _) in ram_points) and ram_total == 62)
print("characters ramified (pole order 1 after reduction, exponent 2) at the two places on A = 0: %d of 62" % ram_total)
print("a_d(k), d = 1..10:", ak[1:DMAX + 1])

# --------------------------------------------------------------- the genus from the zeta functions
def count_all(d):
    F = field(d)
    pts = affine_points(F)
    nh = len(pts) + 2
    nc = [0] * 32
    for (x, y) in pts:
        if F.poly(A_poly, x) == 0:
            for h in range(1, 32):
                nc[h] += 1
            continue
        fv = frob_vector(F, x, y)
        for h in range(1, 32):
            if bin(fv & h).count("1") % 2 == 0:
                nc[h] += 2
    for h in range(1, 32):
        nc[h] += 4
    return nh, nc

def l_poly_predicts(S_list, deg, q=2):
    """S_list[d-1] = sum omega_i^d for d = 1..len. Solve e_1..e_(deg/2) by Newton, complete by the functional
    equation, predict the remaining power sums; returns (ok, coefficients)."""
    half = deg // 2
    e = [Fraction(1)]
    for k in range(1, half + 1):
        e.append(sum((-1) ** (i - 1) * e[k - i] * S_list[i - 1] for i in range(1, k + 1)) / k)
    if any(x.denominator != 1 for x in e):
        return False, None
    l = [(-1) ** k * e[k] for k in range(half + 1)]
    l += [q ** (half - (deg - k)) * l[deg - k] for k in range(half + 1, deg + 1)]   # l_(deg-k) = q^(half-k) l_k
    ee = [(-1) ** k * l[k] for k in range(deg + 1)]
    ok = True
    for d in range(half + 1, len(S_list) + 1):
        # Newton: S_d = sum_(i=1)^(d-1) (-1)^(i-1) e_i S_(d-i) + (-1)^(d-1) d e_d   (e_i = 0 for i > deg)
        pred = sum((-1) ** (i - 1) * (ee[i] if i <= deg else 0) * S_list[d - i - 1] for i in range(1, d))
        pred += (-1) ** (d - 1) * d * (ee[d] if d <= deg else 0)
        if pred != S_list[d - 1]:
            ok = False
    return ok, l

DZ = 14
NH_all, NC_all = {}, {}
for d in range(1, DZ + 1):
    NH_all[d], NC_all[d] = count_all(d)
# H over F_2(x): L-polynomial = numerator of Z_H, degree 2 g(H); S_d = q^d + 1 - N_d
SH = [2 ** d + 1 - NH_all[d] for d in range(1, DZ + 1)]
gH = None
for g in range(1, 6):
    if l_poly_predicts(SH, 2 * g)[0]:
        gH = g
        break
print("g(H) from its zeta function (L-polynomial of degree 2g predicts #H(F_(2^d)) for d <= %d): %s" % (DZ, gH))
deg_ok = 0
for h in range(1, 32):
    Sd = [NH_all[d] - NC_all[d][h] for d in range(1, DZ + 1)]
    ok22, l = l_poly_predicts(Sd, 22)
    ok20, _ = l_poly_predicts(Sd, 20)
    ok24, _ = l_poly_predicts(Sd, 24)
    if ok22 and not ok20 and not ok24:
        deg_ok += 1
print("characters whose L-polynomial has degree 22 (not 20 or 24; counts to d = %d): %d of 31" % (DZ, deg_ok))
gk = None
if cond_ok and deg_ok == 31 and gH == 2:
    gk = (2 * gH + 31 * 22) // 2               # Z_k = Z_H prod_chi L_chi: 2 g(k) = 2 g(H) + sum_chi deg L_chi
print("g(k) =", gk)
ok = inf_ok and cond_ok
if closed_H[1:6] != paper_aH or ak[1:DMAX + 1] != paper_ak or gk != paper_gk:
    ok = False
    print("MISMATCH with the paper's data")
split = (ak[1] == 32 * closed_H[1])
print("all %d rational places of H split completely: %s" % (closed_H[1], split))
ok = ok and split

# --------------------------------------------------------------- Theorem 1.1 (q = 2, nu = 2 on S)
for f, c in S.items():
    if c > ak[f]:
        ok = False
        print("S needs %d places of degree %d, k has %d" % (c, f, ak[f]))
sum1 = sum(c * f for f, c in S.items())                  # e f (nu - 1 - [(nu-1)/2]) = f
sum2 = sum(c * f * (f + 1) for f, c in S.items())        # e f (nu - 1) (e f (nu - 1) + 1)
ineq = (1 + sum1 - t_size) ** 2 - 2 * sum2 - 4 * sum1
print("t = %d <= %d; inequality value %d" % (t_size, sum1, ineq))
ok = ok and t_size <= sum1 and ineq >= 0 and t_size <= ak[1]
g_used = gk if gk else paper_gk
avg = sum(Fraction(2 * c * f) * (1 - Fraction(1, 2 ** f)) for f, c in S.items())
val = Fraction(t_size) / (g_used - 1 + avg / 2)
digits = (val.numerator * 10 ** 20) // val.denominator
print("A(2) >= %s = 0.%020d" % (val, digits))
if claim is not None:
    good = val >= claim
    print("CONTROL claim A(2) >= %s: %s" % (claim, good))
    ok = ok and good
else:
    ok = ok and digits // 10 ** 14 == 316999
print("VERDICT", "PASS" if ok else "REFUSED")
