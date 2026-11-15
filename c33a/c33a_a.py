#!/usr/bin/env python3
"""33a, program A: the finite inputs of [DM2013] A(2) >= 0.316999... (Duursma and Mak, arXiv 1102.4127v4,
Proposition 5.5 and Theorem 1.3 in its numbering; the criterion is its Theorem 1.1), from the paper's text alone.

The base field k is the compositum, over H = F_2(x, y) with y^2 + (x^3 + x + 1) y = x^2 + x, of the 32 Artin-Schreier
extensions v^2 + (x^5 + x^2 + 1) v = (x^2 + x) h, h in the F_2-span of {1, x, x^2, y, y^2}.  Substituting v = A w,
A = x^5 + x^2 + 1, gives w^2 + w = f_h with f_h = (x^2 + x) h / A^2.

This program (route A):
  * F_(2^d) as integers modulo the numerically smallest irreducible polynomial of degree d; carry-less products;
    the trace by d - 1 squarings; w^2 + w = c solved from a table of all w.
  * #H(F_(2^d)) from the affine model plus the two places at infinity (u = 1/x, w = y u^3: w^2 + w = 0 at u = 0).
  * #C_h(F_(2^d)) for the 31 nonzero h by the trace test, and #k(F_(2^d)) = sum_h #C_h - 30 #H (k/H has group
    (Z/2)^5, so Z_k = Z_H prod_(chi != 1) L(chi) and Z_(C_chi) = Z_H L(chi)).  Places by Moebius inversion.
  * The conductor of every nontrivial character at the zeros of A: Laurent series of f_h in t = x - x0 (x - x0 is a
    local parameter there because a(x0) != 0); the pole is reduced by g = sqrt(c_-2)/t; ramified iff
    c_-1 + sqrt(c_-2) != 0, and then the conductor exponent is 2.  g(k) by the conductor-discriminant formula.
  * Theorem 1.1's numerical conditions and the bound, in exact rationals.
Standard library only.  Usage: c33a_a.py [--claim X] [--variant NAME] [--S "5:2,6:16,8:15,10:4"]
"""
import sys
from fractions import Fraction as Fr

# ---------------------------------------------------------------- polynomials over F_2 as Python integers
def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r

def pmod(a, m):
    dm = m.bit_length() - 1
    while a and a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a

def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a

def irreducible(m):
    """Rabin-style test: m of degree n is irreducible iff gcd(m, x^(2^i) - x) = 1 for 1 <= i <= n/2."""
    n = m.bit_length() - 1
    if n < 1:
        return False
    xp = 2  # the polynomial x
    for i in range(1, n // 2 + 1):
        xp = pmod(pmul(xp, xp), m)   # x^(2^i) mod m
        if pgcd(m, xp ^ 2) != 1:
            return False
    return True

def smallest_irreducible(n):
    m = (1 << n) | 1
    while not irreducible(m):
        m += 2
    return m

class GF:
    def __init__(self, d):
        self.d = d
        self.m = smallest_irreducible(d)
        self.q = 1 << d
        self.half = {}
        for w in range(self.q):            # c -> one root w of w^2 + w = c (the other is w + 1)
            self.half.setdefault(self.mul(w, w) ^ w, w)
    def mul(self, a, b):
        return pmod(pmul(a, b), self.m)
    def pw(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r
    def inv(self, a):
        assert a
        return self.pw(a, self.q - 2)
    def tr(self, a):
        s, t = a, a
        for _ in range(self.d - 1):
            t = self.mul(t, t)
            s ^= t
        assert s in (0, 1)
        return s
    def sqrt(self, a):
        return self.pw(a, self.q // 2)
    def ev(self, poly, x):  # poly over F_2 (integer bit mask) evaluated at x
        r = 0
        for i in range(poly.bit_length() - 1, -1, -1):
            r = self.mul(r, x) ^ ((poly >> i) & 1)
        return r

# ---------------------------------------------------------------- the data of the paper (and red-control variants)
P_a = 0b1011          # x^3 + x + 1
P_b = 0b110           # x^2 + x
P_A = 0b100101        # x^5 + x^2 + 1
H_BASIS = ["1", "x", "x^2", "y", "y^2"]
CLAIM_aH = [6, 0, 1, 1, 6]
CLAIM_ak = [192, 0, 0, 0, 2, 16, 0, 16, 0, 64]
CLAIM_gH, CLAIM_gk = 2, 343
S_DEFAULT = {5: 2, 6: 16, 8: 15, 10: 4}
T_SIZE = 192

args = sys.argv[1:]
claim = None
variant = None
S = dict(S_DEFAULT)
i = 0
while i < len(args):
    if args[i] == "--claim":
        claim = Fr(args[i + 1]); i += 2
    elif args[i] == "--variant":
        variant = args[i + 1]; i += 2
    elif args[i] == "--S":
        S = {int(a): int(b) for a, b in (p.split(":") for p in args[i + 1].split(","))}; i += 2
    else:
        sys.exit("unknown argument " + args[i])
if variant == "h-x3":          # y^2 replaced by x^3 in the span
    H_BASIS = ["1", "x", "x^2", "y", "x^3"]
elif variant == "A-x5x3":      # A = x^5 + x^3 + 1 (also irreducible) in place of x^5 + x^2 + 1
    P_A = 0b101001
elif variant is not None:
    sys.exit("unknown variant")

def hval(F, name, x, y):
    return {"1": 1, "x": x, "x^2": F.mul(x, x), "x^3": F.mul(F.mul(x, x), x),
            "y": y, "y^2": F.mul(y, y)}[name]

# ---------------------------------------------------------------- points over F_(2^d)
def points_H(F):
    """affine points (x, y) of y^2 + a(x) y = b(x), each checked nonsingular; plus the count at infinity."""
    pts = []
    for x in range(F.q):
        a, b = F.ev(P_a, x), F.ev(P_b, x)
        if a:
            c = F.mul(b, F.inv(F.mul(a, a)))
            if c in F.half:
                w = F.half[c]
                pts += [(x, F.mul(a, w)), (x, F.mul(a, w ^ 1))]
        else:
            y = F.sqrt(b)
            # F_y = a(x) = 0 here, so smoothness needs F_x = a'(x) y + b'(x) != 0; a' = x^2 + 1, b' = 1
            assert (F.mul(F.mul(x, x) ^ 1, y) ^ 1) != 0, "singular affine point"
            pts.append((x, y))
    for (x, y) in pts:
        assert F.mul(y, y) ^ F.mul(F.ev(P_a, x), y) ^ F.ev(P_b, x) == 0
    n_inf = sum(1 for w in range(F.q) if F.mul(w, w) ^ w == 0)   # u = 0: w^2 + (1 + u^2 + u^3) w = u^4 + u^5
    return pts, n_inf

def trace_vector(F, x, y):
    """(Tr f_h(P))_(h in basis) for an affine point with A(x) != 0, as a 5-bit integer."""
    A = F.ev(P_A, x)
    base = F.mul(F.ev(P_b, x), F.inv(F.mul(A, A)))
    v = 0
    for j, name in enumerate(H_BASIS):
        v |= F.tr(F.mul(base, hval(F, name, x, y))) << j
    return v

def counts(d):
    F = GF(d)
    pts, n_inf = points_H(F)
    NH = len(pts) + n_inf
    NC = [0] * 32       # NC[h] for the 31 nonzero h (bit j = coefficient of H_BASIS[j])
    ram = 0
    for (x, y) in pts:
        if F.ev(P_A, x) == 0:
            ram += 1
            continue
        tv = trace_vector(F, x, y)
        for h in range(1, 32):
            if bin(h & tv).count("1") % 2 == 0:
                NC[h] += 2
    for h in range(1, 32):
        NC[h] += ram + 2 * n_inf      # ramified points (checked below) and the split points at infinity
        # at infinity f_h = (1 + u) u^8 h / (1 + u^3 + u^5)^2 with u^6 h polynomial in (u, w): f_h(u = 0) = 0
    Nk = sum(NC[1:]) - 30 * NH
    return NH, NC, Nk, ram

# ---------------------------------------------------------------- local analysis at the zeros of A
def series_mul(F, s, t, n):
    r = [0] * n
    for i, a in enumerate(s[:n]):
        if a:
            for j, b in enumerate(t[:n - i]):
                r[i + j] ^= F.mul(a, b)
    return r

def series_inv(F, s, n):
    inv0 = F.inv(s[0])
    r = [0] * n
    r[0] = inv0
    for k in range(1, n):
        acc = 0
        for j in range(1, k + 1):
            if j < len(s):
                acc ^= F.mul(s[j], r[k - j])
        r[k] = F.mul(acc, inv0)
    return r

def poly_shift(F, poly, x0, n):
    """coefficients of poly(x0 + t) in t, up to t^(n-1)."""
    c = [0] * n
    xs = [1, 0]          # (x0 + t)^0
    cur = [1] + [0] * (n - 1)
    for k in range(poly.bit_length()):
        if (poly >> k) & 1:
            for i in range(n):
                c[i] ^= cur[i]
        cur = series_mul(F, cur, [x0, 1] + [0] * (n - 2), n)
    return c

def conductor_check():
    """returns (number of zeros of A on H counted as points, number of (point, h) pairs ramified with exponent 2)."""
    F = GF(5)
    assert irreducible(P_A) and P_A.bit_length() - 1 == 5
    n = 6
    found, ram_pairs, places = 0, 0, set()
    for x0 in range(F.q):
        if F.ev(P_A, x0):
            continue
        a0, b0 = F.ev(P_a, x0), F.ev(P_b, x0)
        assert a0 != 0                    # H -> P^1 unramified over x0, so t = x - x0 is a local parameter
        c = F.mul(b0, F.inv(F.mul(a0, a0)))
        if c not in F.half:               # this zero of A is inert in H: no point of H over F_32 above it
            print("zero of A inert in H (the paper's two degree-5 places P_5, P_5' are absent)")
            continue
        for y0 in (F.mul(a0, F.half[c]), F.mul(a0, F.half[c] ^ 1)):
            found += 1
            aa, bb = poly_shift(F, P_a, x0, n), poly_shift(F, P_b, x0, n)
            # y(t) from y^2 + a y = b: a_0 y_k = b_k + [k even] y_(k/2)^2 + sum_(i>=1) a_i y_(k-i)
            ys = [y0] + [0] * (n - 1)
            ia0 = F.inv(aa[0])
            for k in range(1, n):
                acc = bb[k]
                if k % 2 == 0:
                    acc ^= F.mul(ys[k // 2], ys[k // 2])
                for i2 in range(1, k + 1):
                    acc ^= F.mul(aa[i2], ys[k - i2])
                ys[k] = F.mul(acc, ia0)
            xs = [x0, 1] + [0] * (n - 2)
            AA = poly_shift(F, P_A, x0, n + 2)
            assert AA[0] == 0 and AA[1] != 0          # simple zero of A
            q = series_inv(F, AA[1:] , n)            # A = t * (AA[1] + AA[2] t + ...)
            inv_A2 = series_mul(F, q, q, n)          # 1/A^2 = t^-2 * inv_A2
            num_base = series_mul(F, poly_shift(F, P_b, x0, n), [1] + [0] * (n - 1), n)
            hs = {"1": [1] + [0] * (n - 1), "x": xs, "x^2": series_mul(F, xs, xs, n),
                  "x^3": series_mul(F, series_mul(F, xs, xs, n), xs, n), "y": ys, "y^2": series_mul(F, ys, ys, n)}
            basis_series = [series_mul(F, series_mul(F, num_base, hs[nm], n), inv_A2, n) for nm in H_BASIS]
            for h in range(1, 32):
                s = [0] * n
                for j in range(5):
                    if (h >> j) & 1:
                        s = [u ^ v for u, v in zip(s, basis_series[j])]
                cm2, cm1 = s[0], s[1]                 # f_h = cm2 t^-2 + cm1 t^-1 + ...
                if cm1 ^ F.sqrt(cm2):
                    ram_pairs += 1
    return found, ram_pairs

def moebius(n):
    r, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    return -r if m > 1 else r

def places(N, dmax):
    a = []
    for d in range(1, dmax + 1):
        s = sum(moebius(d // e) * N[e] for e in range(1, d + 1) if d % e == 0)
        assert s % d == 0
        a.append(s // d)
    return a

# ---------------------------------------------------------------- run
print("program A (c33a_a.py): F_(2^d) by carry-less products modulo the smallest irreducible; trace by squarings")
if variant:
    print("CONTROL variant:", variant)
ok = True
NH, Nk = {}, {}
for d in range(1, 11):
    nh, nc, nk, ram = counts(d)
    NH[d], Nk[d] = nh, nk
    print("d = %2d: #H = %5d  #k = %7d  ramified points over A = 0: %d" % (d, nh, nk, ram))
aH = places(NH, 10)
ak = places(Nk, 10)
print("a_d(H), d = 1..10:", aH)
print("a_d(k), d = 1..10:", ak)
# genus of H: smooth affine model y^2 + a y = b with deg a = 3 and deg b = 2 <= 6, two rational places at infinity
gH = (P_a.bit_length() - 1) - 1
print("g(H) =", gH, "(deg a - 1 for a smooth model y^2 + a(x) y = b(x) with deg a = 3 > deg b / 2)")
if aH[:5] != CLAIM_aH or gH != CLAIM_gH:
    ok = False; print("MISMATCH with the paper's a_d(H) = (6, 0, 1, 1, 6), g(H) = 2")
split = (Nk[1] == 32 * NH[1])
print("all %d rational places of H split completely in k: %s" % (NH[1], split))
ok &= split
found, ram_pairs = conductor_check()
print("zeros of A on H (points over F_32): %d; (point, character) pairs ramified with conductor exponent 2: %d of %d"
      % (found, ram_pairs, found * 31))
cond_ok = (found == 10 and ram_pairs == found * 31)
if cond_ok:
    deg_f = 2 * 5 + 2 * 5                   # 2 P_5 + 2 P_5' for every nontrivial character
    gk = (32 * (2 * gH - 2) + 31 * deg_f + 2) // 2
else:
    gk = None
print("every nontrivial character has conductor 2P_5 + 2P_5' (degree 20):", cond_ok, "; g(k) =", gk)
if ak != CLAIM_ak or gk != CLAIM_gk:
    ok = False; print("MISMATCH with the paper's a_d(k) = (192,0,0,0,2,16,0,16,0,64), g(k) = 343")
# Theorem 1.1 with q = 2 (e = 1, p = 2), nu = 2 at every place of S
e, p, nu = 1, 2, 2
for f, cnt in S.items():
    if cnt > ak[f - 1]:
        ok = False; print("S asks for %d places of degree %d; k has %d" % (cnt, f, ak[f - 1]))
if T_SIZE > ak[0]:
    ok = False
sig = sum(cnt * e * f * (nu - 1 - (nu - 1) // p) for f, cnt in S.items())
lhs = (1 + sig - T_SIZE) ** 2 - 2 * sum(cnt * e * f * (nu - 1) * (e * f * (nu - 1) + 1) for f, cnt in S.items()) - 4 * sig
print("S:", ", ".join("%d of degree %d" % (c, f) for f, c in sorted(S.items())), "; t =", T_SIZE,
      "<= sum e f (nu-1-[(nu-1)/p]) =", sig, ":", T_SIZE <= sig)
print("inequality (1.1): (1 + %d - %d)^2 - ... = %d >= 0: %s" % (sig, T_SIZE, lhs, lhs >= 0))
ok &= (T_SIZE <= sig and lhs >= 0)
fprime = sum(cnt * 2 * f * (1 - Fr(1, 2 ** f)) for f, cnt in S.items())
gk_used = gk if gk is not None else CLAIM_gk
bound = Fr(T_SIZE) / (gk_used - 1 + fprime / 2)
plain = Fr(T_SIZE) / (gk_used - 1 + sum(cnt * f * nu for f, cnt in S.items()) / Fr(2))
def dec(x, n=20):
    s = str((x.numerator * 10 ** n) // x.denominator).rjust(n + 1, "0")
    return s[:-n] + "." + s[-n:]
print("Corollary 1.2 (all characters at full conductor): A(2) >=", plain, "=", dec(plain))
print("Theorem 1.1 with f' = 2 sum f (1 - 2^-f):         A(2) >=", bound, "=", dec(bound))
if claim is not None:
    print("CONTROL claim A(2) >=", claim, ":", bound >= claim)
    ok &= bound >= claim
else:
    ok &= (dec(bound, 6) == "0.316999")
print("VERDICT", "PASS" if ok else "REFUSED")
