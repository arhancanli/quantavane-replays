#!/usr/bin/env python3
"""Program 2 for C_82 certificates (here PR #280's update of 14 November 2026, v6), written separately from
../c82/c82_check.py (program 1), from the 82a page and the certificate format alone. The two programs share no code,
and no program of the claimant was read or run.

  c82v6_second.py selftest
  c82v6_second.py roots CERT APPROX.json                  root approximations (certified again wherever they are used)
  c82v6_second.py part CERT I N OUT.json pos|jensen [--approx APPROX.json]   the integral over the panels of part I of N
  c82v6_second.py combine CERT PRINTED SHA256 OUT.json ... [--approx APPROX.json]   hypotheses, log M, B, red controls

The certificate: integer polynomials R_j (denominators, multiplicities m_j) and S_i (numerators, exponents q_i >= 0),
T(z) = z(1 - z), u = exp(2 pi i t), D = sum m_j deg R_j,
    g(t) = sum_i q_i log|S_i(T(u))| - sum_j m_j log|R_j(T(u))|,
    B    = [ sum_j m_j log M(R_j o T) + int_0^1 max(0, g(t)) dt ] / D,
and C_82 <= B when every R_j is monic up to sign, m_j >= 1, q_i >= 0, deg A = sum q_i deg S_i < D and S_i, R_j are
coprime whenever q_i > 0 (the certificates' notes; read, not re-derived here).

How each part is decided (no floating-point number decides anything; floats only give starting points):
- Coprimality: Euclid's algorithm modulo a prime (10^9 + 7, else 998244353, else 10^9 + 9) dividing neither leading
  coefficient; a trivial gcd modulo p proves coprimality over Q.
- Roots. P = F o T has N = 2 deg F roots. Starting points: the roots of F (see below), lifted by
  z = (1 +- sqrt(1 - 4w))/2; Newton's method on P in 640-bit fixed point; the results rounded to Gaussian rationals
  z = (x + iy)/2^256. Certification by the Newton inclusion (P'/P = sum 1/(z - root), so some root lies within
  N |P(z)/P'(z)| of z), with P(z) and P'(z) = F'(T(z))(1 - 2z) evaluated EXACTLY as Gaussian integers by
  homogeneous Horner at T(z) = ((x+iy) 2^256 - (x+iy)^2)/2^512; the N discs are checked pairwise disjoint, so each
  holds exactly one root. (Starting points: the roots of F by the Weierstrass iteration, first in floating point on F
  scaled by Fujiwara's bound, then in 256-bit fixed point; if Newton from the lifted points does not separate the
  roots, the Weierstrass iteration on P itself in 640-bit fixed point.)
- log M(R o T) = sum over the discs of log max(1, |z|) (|lc(R o T)| = 1); logarithms from Python's decimal module at
  60 digits (correctly rounded by its specification), widened by one unit in the last place.
- The integral, over [0, 1/2] and doubled (g(1 - t) = g(t), real coefficients), on dyadic panels t = c + h eta,
  |eta| <= 1, starting from 1024 panels:
  * g(c) by Horner's rule on the exact coefficients at T(u_c), in complex ball arithmetic at 2^-192, with
    u_c = exp(2 pi i c) from our own pi (Euler: 20 arctan(1/7) + 8 arctan(3/79)) and Taylor series with bounds;
  * the variation inside the panel from the roots r: log|u - r| = log|u_c - r| + Re log(1 + b_r V(eta)),
    b_r = 2 pi h u_c/(u_c - r), V(eta) = (exp(2 pi i h eta) - 1)/(2 pi h). A root is far when
    x_r = |b_r| Vmax <= 1/2, Vmax = sup_{|eta| <= 4} |V| <= 4 (1 + 8 pi h). For the far roots, the degree-28 Taylor
    polynomial in eta is Re sum_j (-1)^(j+1) B_j V(eta)^j / j with power sums B_j = sum_r c_r b_r^j (V^j truncated);
    the remainder is bounded by Cauchy's estimate on |eta| <= 4, sum |c_r| 2 x_r 4^-29 (4/3), plus, for each root
    whose powers stop at J < 28, (y^(J+1)/((J+1)(1-y))) |c_r| with y = |b_r| (1 + 2 pi h). The polynomial is
    integrated exactly in eta.
  * A panel is negative (sup g <= 0), positive (inf g > 0: the integral of the polynomial, plus the remainder) or
    mixed. Near roots (not far) enter through their distance to the arc only: a numerator root through the largest
    distance, a denominator root through the smallest. A mixed panel is bisected until 2h sup g <= 2^-24, and then
    contributes [0, 2h sup g].
- Self-test of the integrator (mode jensen): int_0^1 sum_j m_j log|R_j(T(u))| dt against sum_j m_j log M(R_j o T)
  from the root moduli (Jensen's formula); and closed forms in `selftest`.
"""
import cmath
import hashlib
import json
import math
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 60
W = 256            # scale of the certified root approximations
WN = 640           # working precision of Newton's method
S = 192            # fixed-point scale of the integral
K = 28             # degree of the Taylor polynomial in eta
GAMMA = 4          # radius of the Cauchy disc in eta
L0 = 11            # first level: panels of width 2^-11 on [0, 1/2]
LMAX = 60
CRUDE = Fraction(1, 1 << 24)
ONE = 1 << S
isqrt = math.isqrt


# ----------------------------------------------------------------------------------------------------------------
# the certificate and the hypotheses
# ----------------------------------------------------------------------------------------------------------------

def read_cert(path):
    with open(path) as fh:
        d = json.load(fh)
    den, num = [], []
    for f in d["denominator"]:
        c = f["coeffs"]
        if not all(type(v) is int for v in c) or type(f["multiplicity"]) is not int:
            raise ValueError("non-integer data in a denominator")
        den.append((f["name"], list(c), f["multiplicity"]))          # descending coefficients
    for f in d["numerator"]:
        c = f["coeffs"]
        if not all(type(v) is int for v in c) or not isinstance(f["q_rational"], str):
            raise ValueError("non-integer data in a numerator")
        num.append((f["name"], list(c), Fraction(f["q_rational"])))
    return den, num


def gcd_degree_mod(a, b, p):
    """Degree of gcd(a, b) over F_p (descending integer lists whose leading coefficients p does not divide)."""
    def norm(v):
        v = [x % p for x in v]
        i = 0
        while i < len(v) and v[i] == 0:
            i += 1
        return v[i:]
    a, b = norm(a), norm(b)
    while b:
        inv = pow(b[0], p - 2, p)
        while a and len(a) >= len(b):
            f = a[0] * inv % p
            a = norm([a[i] - f * b[i] if i < len(b) else a[i] for i in range(len(a))])
        a, b = b, a
    return len(a) - 1


def hypotheses(den, num, rep):
    ok = True
    for name, c, m in den:
        ok &= rep("H1 %s: deg >= 1, leading coefficient +-1, multiplicity a positive integer" % name,
                  len(c) >= 2 and c[0] in (1, -1) and m >= 1, quiet=True)
    for name, c, q in num:
        ok &= rep("H1 %s: deg >= 1, nonzero leading coefficient, q >= 0" % name,
                  len(c) >= 2 and c[0] != 0 and q >= 0, quiet=True)
    rep("H1 every denominator monic up to sign with a positive integer multiplicity, every numerator with q >= 0 "
        "(%d + %d factors)" % (len(den), len(num)), ok)
    D = sum(m * (len(c) - 1) for _, c, m in den)
    degA = sum(q * (len(c) - 1) for _, c, q in num)
    ok &= rep("H2 deg A < D", degA < D, "D = %d, deg A = %s = %.12f D" % (D, degA, float(degA / D)))
    bad, pairs = [], 0
    for nn, s, q in num:
        if q == 0:
            continue
        for dn, r, m in den:
            pairs += 1
            fine = False
            for p in (1000000007, 998244353, 1000000009):
                if s[0] % p and r[0] % p:
                    fine = gcd_degree_mod(s, r, p) == 0
                    break
            if not fine:
                bad.append("%s/%s" % (nn, dn))
    ok &= rep("H3 S_i, R_j coprime for all %d pairs with q_i > 0 (Euclid modulo 10^9+7)" % pairs, not bad,
              "; ".join(bad[:4]))
    return ok, D, degA


# ----------------------------------------------------------------------------------------------------------------
# roots: Newton's method in fixed point; the Newton inclusion with exact Gaussian-integer evaluation
# ----------------------------------------------------------------------------------------------------------------

def deriv_desc(c):
    n = len(c) - 1
    return [c[i] * (n - i) for i in range(n)]


def _mul(a, b, p):
    return ((a[0] * b[0] - a[1] * b[1]) >> p, (a[0] * b[1] + a[1] * b[0]) >> p)


def _horner(c, z, p):
    v = (c[0] << p, 0)
    for a in c[1:]:
        v = _mul(v, z, p)
        v = (v[0] + (a << p), v[1])
    return v


def _newton(c, dc, z, iters=80):
    one = 1 << WN
    for _ in range(iters):
        zz = _mul(z, z, WN)
        t = (z[0] - zz[0], z[1] - zz[1])
        pv = _horner(c, t, WN)
        dv = _mul(_horner(dc, t, WN), (one - 2 * z[0], -2 * z[1]), WN)
        dd = dv[0] * dv[0] + dv[1] * dv[1]
        if dd == 0:
            return z
        st = (((pv[0] * dv[0] + pv[1] * dv[1]) << WN) // dd, ((pv[1] * dv[0] - pv[0] * dv[1]) << WN) // dd)
        z = (z[0] - st[0], z[1] - st[1])
        if abs(st[0]) + abs(st[1]) < (1 << (WN - W - 16)):
            break
    return z


def float_dk(c, iters=1500):
    """Starting points only: the roots of F (descending coefficients) by the Weierstrass (Durand-Kerner) iteration
    in floating point, on G(s) = F(R s)/(lc R^n) with R Fujiwara's bound (so every root of G lies in |s| <= 1),
    from points on a circle."""
    n = len(c) - 1
    lc = c[0]
    # Fujiwara: every root w has |w| <= 2 max_k |c_k / c_0|^(1/k) (k = n with an extra factor 1/2); logs avoid overflow
    lg = max((math.log(abs(c[k])) - math.log(abs(lc))) / k for k in range(1, n + 1) if c[k]) if any(c[1:]) else 0.0
    R = 2.0 * math.exp(lg)
    b = [1.0 + 0j] + [complex(c[k] / lc) / (R ** k) if math.log(R) * k < 700 else 0j for k in range(1, n + 1)]
    z = [0.9 * cmath.exp(1j * (2 * math.pi * k / n + 0.4)) for k in range(n)]
    for _ in range(iters):
        big = 0.0
        for i in range(n):
            pv = 1.0 + 0j
            for v in b[1:]:
                pv = pv * z[i] + v
            q = 1.0 + 0j
            for j in range(n):
                if j != i:
                    q *= (z[i] - z[j])
            if q == 0:
                q = 1e-300
            st = pv / q
            z[i] -= st
            big = max(big, abs(st))
        if big < 1e-15 or big != big:
            break
    return [R * v for v in z]


def compose_T_desc(c):
    """Coefficients (descending) of F(z - z^2), exactly."""
    res = [c[0]]
    for a in c[1:]:
        # res * (z - z^2) + a, with res descending
        nxt = [0] * (len(res) + 2)
        for i, v in enumerate(res):                         # v z^d times (-z^2 + z): degrees shift by 2 and 1
            nxt[i] -= v
            nxt[i + 1] += v
        nxt[-1] += a
        res = nxt
    return res


def fx_dk(P, zs, iters=200):
    """Fallback: the Weierstrass iteration on P (descending, exact integers) in WN-bit fixed point, all roots at once."""
    N = len(P) - 1
    lc = P[0]
    z = list(zs)
    for _ in range(iters):
        big = 0
        for i in range(N):
            pv = (P[0] << WN, 0)
            for v in P[1:]:
                pv = _mul(pv, z[i], WN)
                pv = (pv[0] + (v << WN), pv[1])
            q = (lc << WN, 0)
            for j in range(N):
                if j != i:
                    q = _mul(q, (z[i][0] - z[j][0], z[i][1] - z[j][1]), WN)
            dd = q[0] * q[0] + q[1] * q[1]
            if dd == 0:
                continue
            st = (((pv[0] * q[0] + pv[1] * q[1]) << WN) // dd, ((pv[1] * q[0] - pv[0] * q[1]) << WN) // dd)
            z[i] = (z[i][0] - st[0], z[i][1] - st[1])
            big = max(big, abs(st[0]) + abs(st[1]))
        if big < (1 << (WN - W - 16)):
            break
    return z


WD = 256           # precision of the Weierstrass iteration on F


def fx_dk_F(c, ws, iters=600):
    """The roots of F (descending integer coefficients) by the Weierstrass iteration in WD-bit fixed point, from the
    floating-point points ws; returns floats (only starting points for Newton on P)."""
    n = len(c) - 1
    z = [(int(round(w.real * 2 ** 50)) << (WD - 50), int(round(w.imag * 2 ** 50)) << (WD - 50)) for w in ws]
    lc = c[0] << WD
    for _ in range(iters):
        big = 0
        for i in range(n):
            pv = _horner(c, z[i], WD)
            q = (lc, 0)
            zi = z[i]
            for j in range(n):
                if j != i:
                    q = _mul(q, (zi[0] - z[j][0], zi[1] - z[j][1]), WD)
            dd = q[0] * q[0] + q[1] * q[1]
            if dd == 0:
                q = (1, 0)
                dd = 1
            st = (((pv[0] * q[0] + pv[1] * q[1]) << WD) // dd, ((pv[1] * q[0] - pv[0] * q[1]) << WD) // dd)
            z[i] = (zi[0] - st[0], zi[1] - st[1])
            big = max(big, abs(st[0]) + abs(st[1]))
        if big < (1 << (WD - 160)):
            break
    return [complex(x / 2 ** WD, y / 2 ** WD) for x, y in z]


def approx_roots(c, fallback=False):
    """Approximations (x, y) at scale 2^-W of the 2n roots of P = F o T (F: descending integer coefficients)."""
    ws = fx_dk_F(c, float_dk(c)) if len(c) > 2 else [complex(-c[1] / c[0])]
    dc = deriv_desc(c)
    starts = []
    for w in ws:
        s = cmath.sqrt(1 - 4 * complex(w))
        for z0 in ((1 + s) / 2, (1 - s) / 2):
            starts.append((int(round(z0.real * 2 ** 52)) << (WN - 52), int(round(z0.imag * 2 ** 52)) << (WN - 52)))
    if fallback:
        zs = fx_dk(compose_T_desc(c), starts)
    else:
        zs = [_newton(c, dc, z) for z in starts]
    return [(z[0] >> (WN - W), z[1] >> (WN - W)) for z in zs]


def isqrt_ceil(v):
    r = isqrt(v)
    return r if r * r >= v else r + 1


def certify(c, zs):
    """[(x, y, R)] at scale 2^-W: N pairwise disjoint discs, each holding one root of P = F o T; or None."""
    n = len(c) - 1
    N = 2 * n
    if len(zs) != N:
        return None
    dc = deriv_desc(c)
    sc = 2 * W
    discs = []
    for x, y in zs:
        X = (x << W) - (x * x - y * y)                     # T(z) 2^(2W)
        Y = (y << W) - 2 * x * y
        vr, vi = c[0], 0                                   # F(T(z)) 2^(2W n), exactly
        for k, a in enumerate(c[1:], start=1):
            vr, vi = vr * X - vi * Y, vr * Y + vi * X
            vr += a << (sc * k)
        dr, di = dc[0], 0                                  # F'(T(z)) 2^(2W (n-1)), exactly
        for k, a in enumerate(dc[1:], start=1):
            dr, di = dr * X - di * Y, dr * Y + di * X
            dr += a << (sc * k)
        qr, qi = (1 << W) - 2 * x, -2 * y                  # T'(z) 2^W
        den = (dr * dr + di * di) * (qr * qr + qi * qi)
        if den == 0:
            return None
        # radius in units of 2^-W: N |P| / |P'| 2^W = N sqrt(|v|^2 / (|d|^2 |q|^2))
        t = -((-(N * N) * (vr * vr + vi * vi)) // den)
        discs.append((x, y, isqrt_ceil(t)))
    order = sorted(range(N), key=lambda i: discs[i][0])
    rmax = max(r for _, _, r in discs)
    for a in range(N):
        xa, ya, ra = discs[order[a]]
        for b in range(a + 1, N):
            xb, yb, rb = discs[order[b]]
            if xb - xa > ra + 2 * rmax:                    # every later disc is too far to the right
                break
            s = ra + rb
            if (xa - xb) ** 2 + (ya - yb) ** 2 <= s * s:
                return None
    return discs


APPROX = {}


def roots_of(factors):
    """Certified discs for every factor; the approximations come from APPROX (a cache written by this program's
    `roots` command in the same run) when present, and are certified again here in every case."""
    out = []
    for name, c, w in factors:
        key = name + ":" + hashlib.sha256(str(c).encode()).hexdigest()[:16]
        zs = [tuple(v) for v in APPROX[key]] if key in APPROX else approx_roots(c)
        d = certify(c, zs)
        if d is None:                                      # Newton from the starting points did not separate them
            zs = approx_roots(c, fallback=True)
            d = certify(c, zs)
        APPROX[key] = zs
        out.append((name, c, w, d))
    return out


def roots_digest(rs):
    h = hashlib.sha256()
    for name, c, w, d in rs:
        h.update(("%s|%s|" % (name, w)).encode())
        for x, y, r in (d or []):
            h.update(("%d,%d,%d;" % (x, y, r)).encode())
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------------------------
# logarithms as Fractions
# ----------------------------------------------------------------------------------------------------------------

def _ulp(dv):
    return Fraction(1, 10 ** (getcontext().prec - 1 - dv.adjusted())) if dv != 0 else Fraction(1, 10 ** 70)


_LN2 = []


def ln2():
    if not _LN2:
        v = Decimal(2).ln()
        _LN2.extend([Fraction(v) - _ulp(v), Fraction(v) + _ulp(v)])
    return _LN2


def ln_bounds(v):
    """(lo, hi) Fractions with lo <= ln v <= hi, v a positive integer."""
    sh = max(0, v.bit_length() - 120)
    a, b = (v >> sh, (v >> sh) + 1) if sh else (v, v)
    la, lb = Decimal(a).ln(), Decimal(b).ln()
    L2 = ln2()
    return Fraction(la) - _ulp(la) + sh * L2[0], Fraction(lb) + _ulp(lb) + sh * L2[1]


def log_mahler(discs):
    """Enclosure of sum log max(1, |z|) over discs at scale 2^-W."""
    lo, hi = Fraction(0), Fraction(0)
    L2 = ln2()
    one = 1 << W
    for x, y, R in discs:
        m = isqrt(x * x + y * y)
        a, b = m - R, m + 1 + R
        if b <= one:
            continue
        if a > one:
            lo += ln_bounds(a)[0] - W * L2[1]
        hi += ln_bounds(b)[1] - W * L2[0]
    return lo, hi


# ----------------------------------------------------------------------------------------------------------------
# complex balls at scale 2^-S: (x, y, r) = centre (x + iy) 2^-S, radius r 2^-S
# ----------------------------------------------------------------------------------------------------------------

def b_mul(a, b):
    na, nb = abs(a[0]) + abs(a[1]), abs(b[0]) + abs(b[1])
    return ((a[0] * b[0] - a[1] * b[1]) >> S, (a[0] * b[1] + a[1] * b[0]) >> S,
            ((na * b[2] + nb * a[2] + a[2] * b[2]) >> S) + 3)


def b_abs(a):
    m = isqrt(a[0] * a[0] + a[1] * a[1])
    return max(0, m - a[2]), m + 1 + a[2]


def b_frac(a, q):
    p, d = q.numerator, q.denominator
    return ((a[0] * p) // d, (a[1] * p) // d, -((-a[2] * abs(p)) // d) + 2)


# ----------------------------------------------------------------------------------------------------------------
# pi and exp(2 pi i t)
# ----------------------------------------------------------------------------------------------------------------

PP = S + 64
_PI = []


def pi_bounds():
    """(lo, hi) at scale 2^-PP: pi = 20 arctan(1/7) + 8 arctan(3/79), alternating series, terms floored."""
    if not _PI:
        def atan_frac(p, q):
            lo = hi = 0
            k, num, den = 0, p, q
            while True:
                term = (num << PP) // (den * (2 * k + 1))
                if term == 0:
                    break
                if k % 2 == 0:
                    lo, hi = lo + term, hi + term + 1
                else:
                    lo, hi = lo - term - 1, hi - term
                k, num, den = k + 1, num * p * p, den * q * q
            return lo - 1, hi + 1
        a1, a2 = atan_frac(1, 7), atan_frac(3, 79)
        _PI.extend([20 * a1[0] + 8 * a2[0], 20 * a1[1] + 8 * a2[1]])
    return _PI


def unit_point(num, den):
    """Ball (scale 2^-S) holding exp(2 pi i num/den), 0 <= num/den <= 1/2."""
    plo, phi = pi_bounds()
    th_lo = (2 * num * plo) // den
    th_hi = -((-2 * num * phi) // den)
    th = th_lo
    c, s, term, k = 1 << PP, 0, 1 << PP, 0
    while True:
        k += 1
        term = (term * th) // (k << PP)
        if term == 0:
            break
        r = k % 4
        if r == 1:
            s += term
        elif r == 2:
            c -= term
        elif r == 3:
            s -= term
        else:
            c += term
    # each computed term is within 5 units of the true one (|th| < 4), the true tail after the first zero term is
    # below 12 units, and cos, sin are 1-Lipschitz in theta
    err = 5 * k + 12 + (th_hi - th_lo)
    sh = PP - S
    return (c >> sh, s >> sh, ((2 * err) >> sh) + 3)


# ----------------------------------------------------------------------------------------------------------------
# the integrator
# ----------------------------------------------------------------------------------------------------------------

class Integrand:
    def __init__(self, rooted):
        """rooted: [(name, descending coeffs, weight (Fraction), discs at scale 2^-W)]; weight 0 factors dropped."""
        self.f = []
        for name, c, w, discs in rooted:
            if w:
                roots = [(x >> (W - S), y >> (W - S), (R >> (W - S)) + 2) for x, y, R in discs]
                self.f.append((name, c, Fraction(w), roots))
        self.cache = {}

    def level(self, L):
        """V^j truncated at degree K, for h = 2^-(L+1); 2 pi h as a ball; the bounds Vmaj, Vmax (scale 2^-S)."""
        if L in self.cache:
            return self.cache[L]
        plo, phi = pi_bounds()
        sh = L + 1 + PP - S
        tl, th = (2 * plo) >> sh, ((2 * phi) >> sh) + 1          # 2 pi h at scale 2^-S
        v = [None] * (K + 1)
        lo_m, hi_m = ONE, ONE
        for k in range(1, K + 1):
            if k > 1:
                lo_m, hi_m = (lo_m * tl) // (k << S), -((-hi_m * th) // (k << S))
            mid, rad = (lo_m + hi_m) // 2, (hi_m - lo_m) // 2 + 2
            v[k] = [(mid, 0, rad), (0, mid, rad), (-mid, 0, rad), (0, -mid, rad)][k % 4]   # i^k (2 pi h)^(k-1)/k!
        zero = (0, 0, 0)
        pw = [None, [zero] + v[1:]]
        for j in range(2, K + 1):
            prev, cur = pw[j - 1], [zero] * (K + 1)
            for k in range(j, K + 1):
                acc = zero
                for a in range(j - 1, k):
                    m = b_mul(prev[a], v[k - a])
                    acc = (acc[0] + m[0], acc[1] + m[1], acc[2] + m[2])
                cur[k] = acc
            pw.append(cur)
        vmaj = ONE + th + 1                                      # (e^x - 1)/x <= 1 + x, x = 2 pi h <= 1
        vmax = GAMMA * (ONE + GAMMA * th) + 1                    # (e^(4x) - 1)/x <= 4 (1 + 4x), 4x <= 1
        out = (pw, (tl + th) // 2, (th - tl) // 2 + 2, vmaj, vmax, th)
        self.cache[L] = out
        return out

    def g_centre(self, uc):
        """[lo, hi] (Fractions) for g at the panel centre: Horner on the exact coefficients at T(u_c)."""
        u2 = b_mul(uc, uc)
        T = (uc[0] - u2[0], uc[1] - u2[1], uc[2] + u2[2])
        lo = hi = Fraction(0)
        L2 = ln2()
        for name, c, w, roots in self.f:
            v = (c[0] << S, 0, 0)
            for a in c[1:]:
                v = b_mul(v, T)
                v = (v[0] + (a << S), v[1], v[2])
            alo, ahi = b_abs(v)
            if alo <= 0:
                return None
            llo, lhi = ln_bounds(alo)[0] - S * L2[1], ln_bounds(ahi)[1] - S * L2[0]
            if w > 0:
                lo, hi = lo + w * llo, hi + w * lhi
            else:
                lo, hi = lo + w * lhi, hi + w * llo
        return lo, hi

    def panel(self, L, k, mode):
        """Panel [k 2^-L, (k+1) 2^-L]: ('neg'|'pos'|'crude'|'split', lo, hi), the bounds for its integral of
        max(0, g) (mode 'pos') or of g (mode 'jensen')."""
        uc = unit_point(2 * k + 1, 1 << (L + 1))
        pw, bmid, brad, vmaj, vmax, tph = self.level(L)
        twopih = (bmid, 0, brad)
        chord = tph + 1                                          # |u - u_c| <= 2 pi h on the panel
        tolJ = 1 << (S - 70)
        half = ONE >> 1
        Bs = [(0, 0, 0)] * (K + 1)
        E = Fraction(0)
        near = []
        for name, c, w, roots in self.f:
            ps = [[0, 0, 0] for _ in range(K + 1)]
            ef = 0                                               # this factor's remainder, scale 2^-S
            top = 0
            for rt in roots:
                dx, dy, dr = uc[0] - rt[0], uc[1] - rt[1], uc[2] + rt[2]
                m2 = dx * dx + dy * dy
                m = isqrt(m2)
                if m <= dr + 1:
                    near.append((w, (dx, dy, dr)))
                    continue
                # inv = 1/(u_c - r) as a ball
                ix, iy = (dx << (2 * S)) // m2, (-dy << (2 * S)) // m2
                ir = -((-(dr + 1) << (2 * S)) // (m * (m - dr - 1))) + 3
                bb = b_mul(b_mul(uc, (ix, iy, ir)), twopih)
                bhi = b_abs(bb)[1]
                x = -((-bhi * vmax) >> S)
                if x > half:
                    near.append((w, (dx, dy, dr)))
                    continue
                ef += 2 * x                                      # -log(1 - x) <= 2x, times 4^-(K+1) 4/3 below
                y = -((-bhi * vmaj) >> S)
                p = bb
                J = 1
                acc = ps[1]
                acc[0] += p[0]; acc[1] += p[1]; acc[2] += p[2]
                yj = y
                while J < K:
                    # the terms j > J: sum_{j > J} y^j / j <= y^(J+1) / ((J+1)(1 - y))
                    tail = -((-yj * y) // ((J + 1) * (ONE - y)))
                    if tail <= tolJ:
                        ef += tail * (4 ** (K + 1)) * 3 // 4 + 1  # kept in the same units as the Cauchy terms
                        break
                    J += 1
                    p = b_mul(p, bb)
                    acc = ps[J]
                    acc[0] += p[0]; acc[1] += p[1]; acc[2] += p[2]
                    yj = -((-yj * y) >> S)
                top = max(top, J)
            for j in range(1, top + 1):
                if ps[j][0] or ps[j][1] or ps[j][2]:
                    t = b_frac(tuple(ps[j]), w)
                    Bs[j] = (Bs[j][0] + t[0], Bs[j][1] + t[1], Bs[j][2] + t[2])
            E += abs(w) * Fraction(ef, ONE)
        E = E * Fraction(4, 3) / (4 ** (K + 1))                  # Cauchy: M 4^-(K+1) / (1 - 1/4)
        # the polynomial: a_k = sum_{j <= k} (-1)^(j+1) B_j [V^j]_k / j; its real parts
        rlo, rhi = [0] * (K + 1), [0] * (K + 1)
        for kk in range(1, K + 1):
            ax = ay = ar = 0
            for j in range(1, kk + 1):
                if Bs[j] == (0, 0, 0):
                    continue
                t = b_mul(Bs[j], pw[j][kk])
                q = Fraction(1 if j % 2 else -1, j)
                t = b_frac(t, q)
                ax, ay, ar = ax + t[0], ay + t[1], ar + t[2]
            rlo[kk], rhi[kk] = ax - ar, ax + ar
        sumabs = Fraction(sum(max(abs(rlo[kk]), abs(rhi[kk])) for kk in range(1, K + 1)), ONE)
        h = Fraction(1, 1 << (L + 1))
        gc = self.g_centre(uc)
        if gc is None:
            return ("split", None, None)
        if near:
            if mode == "jensen":
                return ("split", None, None)
            sup = gc[1] + sumabs + E
            for w, d in near:
                dlo, dhi = b_abs(d)
                if w > 0:
                    if dlo <= 0:
                        return ("split", None, None)
                    sup += w * (ln_bounds(dhi + chord)[1] - ln_bounds(dlo)[0])
                else:
                    if dlo <= chord:
                        return ("split", None, None)
                    sup += (-w) * (ln_bounds(dhi)[1] - ln_bounds(dlo - chord)[0])
            if sup <= 0:
                return ("neg", Fraction(0), Fraction(0))
            if 2 * h * sup <= CRUDE:
                return ("crude", Fraction(0), 2 * h * sup)
            return ("split", None, None)
        evlo = sum(Fraction(2 * rlo[kk], (kk + 1) * ONE) for kk in range(2, K + 1, 2))
        evhi = sum(Fraction(2 * rhi[kk], (kk + 1) * ONE) for kk in range(2, K + 1, 2))
        ilo = h * (2 * gc[0] + evlo - 2 * E)
        ihi = h * (2 * gc[1] + evhi + 2 * E)
        if mode == "jensen":
            return ("pos", ilo, ihi)
        sup = gc[1] + sumabs + E
        inf = gc[0] - sumabs - E
        if sup <= 0:
            return ("neg", Fraction(0), Fraction(0))
        if inf > 0:
            return ("pos", ilo, ihi)
        if 2 * h * sup <= CRUDE:
            return ("crude", Fraction(0), 2 * h * sup)
        return ("split", None, None)

    def integrate(self, part, nparts, mode):
        lo = hi = Fraction(0)
        counts = {"neg": 0, "pos": 0, "crude": 0, "split": 0}
        fails, maxL = 0, L0
        crude_w = Fraction(0)
        t0 = time.time()
        for k0 in range(1 << (L0 - 1)):
            if k0 % nparts != part:
                continue
            stack = [(L0, k0)]
            while stack:
                L, k = stack.pop()
                kind, a, b = self.panel(L, k, mode)
                counts[kind] += 1
                maxL = max(maxL, L)
                if kind == "split":
                    if L >= LMAX:
                        fails += 1
                        continue
                    stack += [(L + 1, 2 * k + 1), (L + 1, 2 * k)]
                    continue
                lo, hi = lo + a, hi + b
                if kind == "crude":
                    crude_w += b - a
        return {"lo": lo, "hi": hi, "counts": counts, "fails": fails, "maxL": maxL, "crude": crude_w,
                "time": time.time() - t0}


# ----------------------------------------------------------------------------------------------------------------
# reporting, parts, combination, self-tests
# ----------------------------------------------------------------------------------------------------------------

class Rep:
    def __init__(self):
        self.ok = True

    def __call__(self, desc, cond, detail="", quiet=False):
        cond = bool(cond)
        self.ok &= cond
        if not quiet or not cond:
            print("  check %s: %s%s" % (desc, "PASS" if cond else "FAIL", (" (%s)" % detail) if detail else ""))
        return cond


def dec(fr, nd=16):
    s = "-" if fr < 0 else ""
    fr = abs(fr)
    q = fr.numerator * 10 ** nd // fr.denominator
    t = str(q).rjust(nd + 1, "0")
    return s + t[:-nd] + "." + t[-nd:]


def factors_for(den, num, mode):
    if mode == "jensen":
        return [(n, c, -m) for n, c, m in den]
    return [(n, c, -m) for n, c, m in den] + [(n, c, q) for n, c, q in num if q]


def run_part(cert, part, nparts, out, mode):
    den, num = read_cert(cert)
    t0 = time.time()
    rooted = roots_of(factors_for(den, num, mode))
    if any(d is None for _, _, _, d in rooted):
        bad = [n for n, _, _, d in rooted if d is None]
        json.dump({"error": "roots not certified: %s" % bad}, open(out, "w"))
        print("roots not certified:", bad)
        return 1
    troots = time.time() - t0
    res = Integrand(rooted).integrate(part, nparts, mode)
    js = {"mode": mode, "part": part, "nparts": nparts, "lo": str(res["lo"]), "hi": str(res["hi"]),
          "counts": res["counts"], "fails": res["fails"], "maxL": res["maxL"], "crude": str(res["crude"]),
          "roots_digest": roots_digest(rooted), "time_roots": round(troots, 1), "time": round(res["time"], 1),
          "cert_sha256": hashlib.sha256(open(cert, "rb").read()).hexdigest()}
    with open(out, "w") as fh:
        json.dump(js, fh, indent=1)
    print("part %d/%d (%s): [%s, %s], %s, fails %d, depth %d..%d, roots %.1f s, integral %.1f s" % (
        part, nparts, mode, dec(res["lo"], 8), dec(res["hi"], 8), res["counts"], res["fails"], L0, res["maxL"],
        troots, res["time"]))
    return 0 if res["fails"] == 0 else 1


def run_combine(cert, printed, sha, parts):
    t0 = time.time()
    rep = Rep()
    h = hashlib.sha256(open(cert, "rb").read()).hexdigest()
    rep("SHA-256 of the certificate is the stated %s..." % sha[:12], h == sha, h)
    den, num = read_cert(cert)
    print("  factors: %d denominators (%s), %d numerators" % (
        len(den), ", ".join("deg %d ^ %d" % (len(c) - 1, m) for _, c, m in den), len(num)))
    ok, D, degA = hypotheses(den, num, rep)
    rooted = roots_of(factors_for(den, num, "pos"))
    nd = sum(1 for _, _, _, d in rooted if d is not None)
    rep("E1 every root of every F o T in its own certified disc (Newton inclusion, %d factors, %d discs, pairwise "
        "disjoint per factor)" % (len(rooted), sum(len(d) for _, _, _, d in rooted if d)), nd == len(rooted),
        "largest radius 2^%d" % max((max(r for _, _, r in d).bit_length() - W) for _, _, _, d in rooted if d))
    dig_pos = roots_digest(rooted)
    dig_jen = roots_digest(rooted[:len(den)])                 # the same denominators, weights and discs
    lmlo = lmhi = Fraction(0)
    for (name, c, w, d), (_, _, m) in zip(rooted[:len(den)], den):
        a, b = log_mahler(d)
        lmlo, lmhi = lmlo + m * a, lmhi + m * b
    print("  sum m_j log M(R_j o T): [%s, %s]" % (dec(lmlo, 10), dec(lmhi, 10)))
    I = {"pos": [Fraction(0), Fraction(0), 0, set(), {}], "jensen": [Fraction(0), Fraction(0), 0, set(), {}]}
    for p in parts:
        js = json.load(open(p))
        if "error" in js:
            rep("part file %s" % p, False, js["error"])
            continue
        md = js["mode"]
        rep("part %d/%d (%s) used the same certificate and the same certified roots" % (js["part"], js["nparts"], md),
            js["cert_sha256"] == h and js["roots_digest"] == (dig_pos if md == "pos" else dig_jen))
        rep("part %d/%d (%s): every panel resolved by depth %d" % (js["part"], js["nparts"], md, LMAX),
            js["fails"] == 0, "%s, depth %d..%d, %.0f s" % (js["counts"], L0, js["maxL"], js["time"]))
        acc = I[md]
        acc[0] += Fraction(js["lo"])
        acc[1] += Fraction(js["hi"])
        acc[3].add((js["part"], js["nparts"]))
        for kk, v in js["counts"].items():
            acc[4][kk] = acc[4].get(kk, 0) + v
    for md in ("pos", "jensen"):
        got = I[md][3]
        n = {b for _, b in got}
        rep("the %s parts cover every panel" % md, len(n) == 1 and got == {(i, list(n)[0]) for i in range(list(n)[0])},
            str(sorted(got)))
    Ilo, Ihi = 2 * I["pos"][0], 2 * I["pos"][1]
    print("  int_0^1 max(0, g) dt: [%s, %s]  (%s)" % (dec(Ilo, 8), dec(Ihi, 8), I["pos"][4]))
    Jlo, Jhi = 2 * I["jensen"][0], 2 * I["jensen"][1]
    rep("self-test (Jensen): int_0^1 sum_j (-m_j) log|R_j(T(u))| dt overlaps -sum_j m_j log M(R_j o T)",
        Jlo <= -lmlo and -lmhi <= Jhi, "[%s, %s] vs [%s, %s]" % (dec(Jlo, 6), dec(Jhi, 6), dec(-lmhi, 6), dec(-lmlo, 6)))
    Blo, Bhi = (lmlo + Ilo) / D, (lmhi + Ihi) / D
    print("  B: [%s, %s]  (width %.2e)" % (dec(Blo, 16), dec(Bhi + Fraction(1, 10 ** 16), 16), float(Bhi - Blo)))
    pr = Fraction(printed)
    rep("B < %s, the printed upper bound" % printed, Bhi < pr, "margin %.3e" % float(pr - Bhi))
    # red controls: each must be refused
    low = Fraction(int(Blo * 10 ** 13) - 1, 10 ** 13)
    rep("red control 1: the bound %s, just below our enclosure, is not implied" % dec(low, 13), not Bhi < low)
    class Q:
        def __init__(self):
            self.ok = True
        def __call__(self, d, c, detail="", quiet=False):
            self.ok &= bool(c)
            return c
    n0, c0, m0 = den[0]
    q = Q(); hypotheses(den, num + [("copy", c0, Fraction(1, 1000))], q)
    rep("red control 2: %s added as a numerator (q = 1/1000) is refused by H3" % n0, not q.ok)
    q = Q(); hypotheses([(n0, [2 * v for v in c0], m0)] + den[1:], num, q)
    rep("red control 3: %s doubled (not monic) is refused by H1" % n0, not q.ok)
    sc = Fraction(D + 1) / degA
    q = Q(); hypotheses(den, [(n, c, qq * sc) for n, c, qq in num], q)
    rep("red control 4: every q_i scaled so that deg A = D + 1 is refused by H2", not q.ok)
    print("  run time (combine): %.1f s" % (time.time() - t0))
    print("  RESULT: %s" % ("ALL CHECKS PASS" if rep.ok else "FAILED"))
    return 0 if rep.ok else 1


def selftest():
    rep = Rep()
    plo, phi = pi_bounds()
    ref = Fraction(314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798, 10 ** 101)
    rep("pi enclosed (width %d units of 2^-%d)" % (phi - plo, PP), Fraction(plo, 1 << PP) <= ref <= Fraction(phi, 1 << PP)
        and phi - plo < 1 << 12)
    u = unit_point(1, 6 * 1)                                     # exp(2 pi i / 6) = (1 + i sqrt 3)/2: not dyadic, a test
    x, y, r = u
    rep("exp(2 pi i/6) = (1 + i sqrt(3))/2 inside the ball", abs(2 * x - ONE) <= 2 * r + 2 and
        abs(4 * y * y - 3 * ONE * ONE) <= 8 * (abs(y) + ONE) * (r + 1), "radius %d units" % r)
    u = unit_point(1, 4)
    rep("exp(2 pi i/4) = i inside the ball", abs(u[0]) <= u[2] and abs(u[1] - ONE) <= u[2])
    l2 = ln2()
    ref = Fraction(693147180559945309417232121458176568075500134360255254120680009, 10 ** 63)
    rep("ln 2 enclosed", l2[0] <= ref <= l2[1])
    # roots and Mahler measures with closed forms: F = x - 2 (|roots of z^2 - z + 2|^2 = 2), F = 2x - 1 (roots
    # (1 +- i)/2 inside), F = x^2 - x - 1
    for c, val, what in (([1, -2], (ln2()[0], ln2()[1]), "log M((x - 2) o T) = ln 2"),
                         ([2, -1], (Fraction(0), Fraction(0)), "the roots of (2x - 1) o T lie inside the circle")):
        d = certify(c, approx_roots(c))
        lm = log_mahler(d) if d else None
        rep(what, lm is not None and lm[0] <= val[1] and val[0] <= lm[1], "%s" % (lm and [dec(lm[0], 12), dec(lm[1], 12)]))
    # Jensen on closed forms: int_0^1 log|T(u) - 2| dt = ln 2, int_0^1 log|2T(u) - 1| dt = ln 2
    for c in ([1, -2], [2, -1]):
        d = certify(c, approx_roots(c))
        it = Integrand([("t", c, 1, d)]).integrate(0, 1, "jensen")
        lo, hi = 2 * it["lo"], 2 * it["hi"]
        rep("int_0^1 log|F(T(u))| dt = ln 2 for F = %s (integrator, Jensen)" % c, lo <= ln2()[1] and ln2()[0] <= hi
            and it["fails"] == 0, "[%s, %s], %s" % (dec(lo, 14), dec(hi, 14), it["counts"]))
    # a deliberately wrong disc is refused: shift one approximation by 2^-100
    c = [1, -3, 1, 7]
    zs = approx_roots(c)
    good = certify(c, zs) is not None
    zs2 = list(zs)
    zs2[0] = (zs2[0][0] + (1 << (W - 100)), zs2[0][1])
    d2 = certify(c, zs2)
    rep("certification accepts a cubic's roots, and a root moved by 2^-100 gives a disc of radius >= 2^-101",
        good and (d2 is None or max(r for _, _, r in d2) >= (1 << (W - 101))))
    zs3 = list(zs)
    zs3[1] = zs3[0]
    rep("two equal approximations are refused (discs not disjoint)", certify(c, zs3) is None)
    print("  RESULT: %s" % ("ALL CHECKS PASS" if rep.ok else "FAILED"))
    return 0 if rep.ok else 1


def main():
    a = sys.argv[1:]
    if "--approx" in a:
        i = a.index("--approx")
        path = a[i + 1]
        a = a[:i] + a[i + 2:]
        with open(path) as fh:
            APPROX.update(json.load(fh))
    if len(a) == 3 and a[0] == "roots":
        den, num = read_cert(a[1])
        t0 = time.time()
        rs = roots_of(factors_for(den, num, "pos"))
        with open(a[2], "w") as fh:
            json.dump({k: [list(z) for z in v] for k, v in APPROX.items()}, fh)
        print("approximations for %d factors, %d certified, %.1f s" % (len(rs), sum(1 for r in rs if r[3]),
                                                                       time.time() - t0))
        return 0 if all(r[3] for r in rs) else 1
    if a and a[0] == "selftest":
        print("self-tests of the arithmetic, the root certification and the integrator")
        return selftest()
    if len(a) == 6 and a[0] == "part":
        return run_part(a[1], int(a[2]), int(a[3]), a[4], a[5])
    if len(a) >= 5 and a[0] == "combine":
        print("PR #280 [Pat26] updated 14 Nov: C_82 < %s* (tao/certificates/c82/v6 at e3bae2e), program 2" % a[2])
        return run_combine(a[1], a[2], a[3], a[4:])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())
