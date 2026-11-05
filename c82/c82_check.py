#!/usr/bin/env python3
"""Independent replay of upper bounds for C_82, the essential minimum of the Zhang-Zagier height
h_Z(a) = h(a) + h(1 - a) (teorth/optimizationproblems, constant 82a), from the certificates' data files.

Written from the 82a page and the certificates' own notes (the [Pin26] README, the degree-44 PROOF.md and METHOD.md,
the [Epi26] README and result.json, the PR #280 README). The claimants' checkers were neither read nor run.
Standard library only (Python 3.9+); run with python3 -I. No floating-point number decides anything: floats only
pick starting points for root refinement and truncation orders; every decision is an exact integer comparison.

  c82_check.py ROW CERTS_ROOT [--controls] [--workers=N]       ROW in: pin epi pr280 pr280b pr280c pr301 pr301a doche

The certificate (every note states the same): integer polynomials R_j (denominators) with multiplicities m_j >= 1,
integer polynomials S_i (numerators) with rational exponents q_i >= 0, coefficient arrays in descending order.
With T(z) = z(1 - z), D = sum_j m_j deg R_j, deg A = sum_i q_i deg S_i, u = exp(2 pi i t),
  g(t) = sum_i q_i log|S_i(T(u))| - sum_j m_j log|R_j(T(u))|,
  B    = [ sum_j m_j log M(R_j o T) + int_0^1 max(0, g(t)) dt ] / D.
Conditions checked: H1 integer coefficients, every R_j monic up to sign, m_j positive integers, q_i >= 0 exact
rationals; H2 deg A < D; H3 gcd(S_i, R_j) = 1 over Q whenever q_i > 0. Then:
  E1 every root of every composed polynomial F o T is enclosed in a disc (all deg(F o T) of them, pairwise disjoint);
  E2 log M(Q o T) = sum_j m_j sum_roots log max(1, |r|), enclosed;
  E3 int_0^1 max(0, g) enclosed by an adaptive dyadic tiling of [0, 1/2] (g(1 - t) = g(t));
  E4 the enclosure of B against the printed bound (an upper bound: B_hi must be below it).

How each is decided
- Arithmetic. Fixed-point integers at scale 2^-128 for the integral, and at 2^-(Mb + 300) for root certification,
  Mb bounding the size of Horner's partial sums (about 300 to 520 bits); every rounding is directed or counted in a
  radius. Logarithms are our own: ln n = e ln 2 + ln(1 + j/256) + 2 atanh(z), z < 1/512,
  series rounded outward with the tail bounded. pi by Machin's formula with alternating-series bounds; exp(2 pi i c)
  by Taylor series with counted rounding and a Lagrange remainder.
- H3 by gcd modulo a prime p not dividing either leading coefficient (a common factor over Q survives reduction
  mod p, so a trivial gcd mod p proves coprimality); a reported failure is confirmed by an exact gcd over Q.
- E1 by Gershgorin's theorem: for approximations z_1..z_n (refined by our own Aberth iteration), the monic F o T is
  the characteristic polynomial of diag(z) - 1 w^T with w_i = F(z_i)/(lc prod_{j != i}(z_i - z_j)); its column discs
  D(z_i - w_i, (n - 1)|w_i|), enclosed in ball arithmetic, are checked pairwise disjoint, so each holds one root.
- E3. On a panel |t - c| <= h, u = u_c e^{2 pi i s}; with eps = e^{2 pi i s} - 1 (|eps| <= 2 pi h <= e, e dyadic)
  and b_r = e u_c/(u_c - r), every root r with e/|u_c - r| <= 1/4 ("far") gives
  log|u - r| = log|u_c - r| + Re sum_k (-1)^{k+1} (b_r eta)^k / k, eta = eps/e, |eta| <= 1,
  truncated at an order K_r per root with remainder <= |b_r|^{K_r+1} / ((K_r + 1)(1 - |b_r|)). The panel integral of
  eta^k is j_k * 2h with j_k a real series in theta = 2 pi h, enclosed once per level. Roots closer than that ("near")
  enter through distance bounds only. A panel is decided negative (sup g <= 0, contributes 0), positive (inf g > 0,
  contributes the enclosed integral of g) or is bisected; an undecided panel (a sign change) contributes
  [max(0, enclosed integral of g), 2h max(0, sup g)] once that width is below 2^-40, or at the depth limit 2^-44.
- Self-tests on every row: g at random points from the certified roots against Horner's rule on the exact
  coefficients; Jensen's formula int_0^1 g = sum w ln M(F o T) for the denominators and for the numerators with no
  root near the circle, the left side from the panel series over every panel, the right side from the root moduli.
  On [Pin26] also a second tiling (first level 2^-11, remainder target 2^-20) and the red controls R5-R7.
What rests on mathematics rather than computation: that C_82 <= B for every certificate meeting H1-H3. The
certificates cite Doche 2001b eq. (9); their notes also give a direct argument (Jensen's formula, multiplicativity
of the Mahler measure, continuity of roots, coprimality of the cyclotomic quotients G_p for distinct primes), which
we read and checked by hand; it is not re-derived by this program.
"""
import cmath
import hashlib
import json
import math
import os
import random
import re
import sys
import time
from fractions import Fraction
from math import isqrt

# ----------------------------------------------------------------------------------------------------------------
# Fixed point: an integer v stands for v / 2^P; real intervals are pairs (lo, hi)
# ----------------------------------------------------------------------------------------------------------------

P = 128
ONE = 1 << P
LP = P + 40                     # internal scale of the logarithms
SH = LP - P
PH = P + 48                     # internal scale of exp(2 pi i c)
ROOT_BITS = 300                 # root certification works at 2^-(Mb + ROOT_BITS), Mb bounding Horner's partials
LS = 72                         # panel contributions are summed at scale 2^-(P + LS)


def cdiv(a, b):
    """ceil(a / b) for b > 0."""
    return -((-a) // b)


def ivq(iv, p, q):
    """The interval iv times the exact rational p/q (q > 0), rounded outward."""
    lo, hi = iv
    if p >= 0:
        return (lo * p) // q, cdiv(hi * p, q)
    return (hi * p) // q, cdiv(lo * p, q)


def ivmul(a, b, scale=P):
    """Product of two intervals at the given scale, outward."""
    c = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(c) >> scale, cdiv(max(c), 1 << scale)


# ---- logarithms (rigorous) ----

def _atanh(a, b):
    """(lo, hi) with lo <= 2^LP atanh(a/b) <= hi, for integers 0 <= 2a <= b."""
    if a == 0:
        return 0, 0
    z = (a << LP) // b
    zl, zh = z, z + 1
    z2l = (zl * zl) >> LP
    z2h = cdiv(zh * zh, 1 << LP)
    lo, t, k = 0, zl, 1
    while t:
        lo += t // k
        t = (t * z2l) >> LP
        k += 2
    hi, t, k = 0, zh, 1
    while t > 1:
        hi += cdiv(t, k)
        t = cdiv(t * z2h, 1 << LP)
        k += 2
    return lo, hi + 2 * t + 2           # the omitted terms sum to at most (4/3) t / k, plus the last roundings


_L2 = _atanh(1, 3)
LN2_LP = (2 * _L2[0], 2 * _L2[1])
_TAB = [tuple(2 * v for v in _atanh(j, 512 + j)) for j in range(256)]   # ln(1 + j/256) = 2 atanh(j/(512 + j))


def ln_int_lp(n):
    """ln n for an integer n >= 1, as an interval at scale 2^-LP."""
    e = n.bit_length() - 1
    j = ((n << 8) >> e) - 256
    a = (n << 8) - ((256 + j) << e)
    b = (n << 8) + ((256 + j) << e)
    zl, zh = _atanh(a, b)
    return e * LN2_LP[0] + _TAB[j][0] + 2 * zl, e * LN2_LP[1] + _TAB[j][1] + 2 * zh


def ln_pos(mlo, mhi, ex):
    """ln of the real interval [mlo 2^ex, mhi 2^ex] (integers 1 <= mlo <= mhi), at scale 2^-P."""
    a = ln_int_lp(mlo)[0]
    b = ln_int_lp(mhi)[1]
    if ex >= 0:
        a += ex * LN2_LP[0]
        b += ex * LN2_LP[1]
    else:
        a += ex * LN2_LP[1]
        b += ex * LN2_LP[0]
    return a >> SH, cdiv(b, 1 << SH)


def ln_fix(v):
    """ln(v / 2^P) for an integer v >= 1, at scale 2^-P."""
    return ln_pos(v, v, -P)


LN2 = (LN2_LP[0] >> SH, cdiv(LN2_LP[1], 1 << SH))


# ---- pi and exp(2 pi i c) ----

def atan_inv(m, bits):
    """Exact rationals lo <= atan(1/m) <= hi with hi - lo < 2^-bits (alternating series, decreasing terms)."""
    S = Fraction(0)
    i = 0
    m2 = m * m
    pw = m
    while True:
        a = Fraction(1, (2 * i + 1) * pw)
        S2 = S + a if i % 2 == 0 else S - a
        if i >= 1 and a < Fraction(1, 2 ** bits):
            return min(S, S2), max(S, S2)
        S = S2
        i += 1
        pw *= m2


def pi_machin(bits):
    a = atan_inv(5, bits)
    b = atan_inv(239, bits)
    return 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]


PI_FR = pi_machin(P + 200)
_PIQ = PH + 32
PI_LO_I = (PI_FR[0].numerator << _PIQ) // PI_FR[0].denominator     # floor(pi 2^_PIQ)
PI_HI_I = cdiv(PI_FR[1].numerator << _PIQ, PI_FR[1].denominator)
EM = cdiv(PI_FR[1].numerator << 16, PI_FR[1].denominator)          # ceil(pi 2^16): e = EM 2^-(L+16) >= pi 2^-L


def unit_point(num, den):
    """exp(2 pi i num/den) for 0 <= num/den <= 1/2, as a complex ball (x, y, r) at scale 2^-P."""
    xl = (2 * PI_LO_I * num) // (den << 32)              # 2 pi num/den at scale 2^-PH, outward
    xh = cdiv(2 * PI_HI_I * num, den << 32)
    xm = (xl + xh) >> 1
    xw = max(xh - xm, xm - xl)
    one = 1 << PH
    c = s = 0
    t, E, errsum, n = one, 0, 0, 0
    while True:
        r = n & 3
        if r == 0:
            c += t
        elif r == 1:
            s += t
        elif r == 2:
            c -= t
        else:
            s -= t
        errsum += E
        n += 1
        tn = (t * xm) // (n << PH)
        E = cdiv(E * xm, n << PH) + 1
        if tn == 0 and n > 4 and (E << 1) < one:
            trunc = E + cdiv((t + E) * xm, n << PH) + 1      # bounds the true next term x^n/n! (Lagrange)
            break
        t = tn
    rad = 2 * (errsum + trunc) + xw + 2
    d = PH - P
    return c >> d, s >> d, cdiv(rad, 1 << d) + 2


# ----------------------------------------------------------------------------------------------------------------
# Integer polynomials (ascending coefficient lists inside the program)
# ----------------------------------------------------------------------------------------------------------------

def compose_T(asc):
    """The ascending coefficients of F(z(1 - z)) for F given by ascending integer coefficients."""
    res = [0]
    for c in reversed(asc):                 # Horner: res = res * (z - z^2) + c
        new = [0] * (len(res) + 2)
        for i, x in enumerate(res):
            if x:
                new[i + 1] += x
                new[i + 2] -= x
        new[0] += c
        while len(new) > 1 and new[-1] == 0:
            new.pop()
        res = new
    return res


def p_deriv(a):
    return [i * a[i] for i in range(1, len(a))]


def fp_trim(a, p):
    a = [c % p for c in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def fp_divmod(a, b, p):
    a = list(a)
    inv = pow(b[-1], p - 2, p)
    q = [0] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and a:
        f = a[-1] * inv % p
        s = len(a) - len(b)
        q[s] = f
        for i in range(len(b)):
            a[s + i] = (a[s + i] - f * b[i]) % p
        while a and a[-1] == 0:
            a.pop()
    return q, a


def fp_gcd(a, b, p):
    a, b = fp_trim(a, p), fp_trim(b, p)
    while b:
        a, b = b, fp_divmod(a, b, p)[1]
    if a:
        inv = pow(a[-1], p - 2, p)
        a = [c * inv % p for c in a]
    return a


PRIMES = [2305843009213693951, 4611686018427387847, 9223372036854775783, 1000000007, 998244353, 2147483647]


def coprime_modp(a, b):
    """(True, p) when gcd(a, b) mod p is constant for a prime p dividing neither leading coefficient: then a, b are
    coprime over Q. (False, None) when every prime tried gives a non-constant gcd."""
    for p in PRIMES:
        if a[-1] % p == 0 or b[-1] % p == 0:
            continue
        if len(fp_gcd(a, b, p)) == 1:
            return True, p
    return False, None


def gcd_Q(a, b):
    """Exact monic gcd over Q of two integer polynomials (ascending), by Euclid over the rationals."""
    A = [Fraction(c) for c in a]
    B = [Fraction(c) for c in b]
    while B and any(B):
        while B and B[-1] == 0:
            B.pop()
        R = list(A)
        while len(R) >= len(B) and R:
            f, s = R[-1] / B[-1], len(R) - len(B)
            for i in range(len(B)):
                R[s + i] -= f * B[i]
            R.pop()
            while R and R[-1] == 0:
                R.pop()
        A, B = B, R
    return [c / A[-1] for c in A]


def squarefree_modp(a):
    for p in PRIMES:
        if a[-1] % p == 0:
            continue
        if len(fp_gcd(a, p_deriv(a), p)) == 1:
            return True, p
    return False, None


# ----------------------------------------------------------------------------------------------------------------
# Roots: our own Aberth iteration (floats, then fixed point at 2^-(Mb + 64) and 2^-(Mb + 300)), then Gershgorin
# ----------------------------------------------------------------------------------------------------------------

def float_aberth(a, maxit=600, z=None):
    """Approximate roots in floats (starting points only)."""
    n = len(a) - 1
    lc = float(a[-1])
    c = [x / lc for x in a]
    if z is None:
        nz = 0
        while nz < n and a[nz] == 0:
            nz += 1
        r0 = abs(c[nz]) ** (1.0 / (n - nz)) if nz < n else 1.0
        r0 = min(max(r0, 0.5), 2.0)
        z = [r0 * cmath.exp(1j * (2 * math.pi * k / n + 0.4)) for k in range(n)]
    z = list(z)
    for it in range(maxit):
        prev = list(z)
        try:
            mx = _aberth_sweep(c, z, n)
        except (OverflowError, ZeroDivisionError):
            return prev
        if not all(math.isfinite(x.real) and math.isfinite(x.imag) for x in z):
            return prev
        if mx < 1e-12:
            break
    return z


def _aberth_sweep(c, z, n):
    mx = 0.0
    for k in range(n):
        zk = z[k]
        p = 0j
        dp = 0j
        for co in reversed(c):
            dp = dp * zk + p
            p = p * zk + co
        if p == 0:
            continue
        if dp == 0:
            z[k] = zk + 1e-8
            continue
        N = p / dp
        s = 0j
        for j in range(n):
            if j != k:
                dz = zk - z[j]
                if dz != 0:
                    s += 1 / dz
        den = 1 - N * s
        w = N / den if den != 0 else N
        z[k] = zk - w
        mx = max(mx, abs(w) / max(1.0, abs(zk)))
    return mx


def _to_fixed(z, prec):
    return int(round(z.real * 2.0 ** 52)) << (prec - 52), int(round(z.imag * 2.0 ** 52)) << (prec - 52)


def _to_float(zf, prec):
    return complex(zf[0] / 2.0 ** prec, zf[1] / 2.0 ** prec)


def fx_aberth(a, zs, prec, tol_ulps, maxit=400):
    """Aberth iteration in fixed point at scale 2^-prec from integer pairs zs (scale 2^-prec). Returns
    (approximations, sweeps) once two consecutive sweeps move no point by more than tol_ulps, else (zs, None)."""
    n = len(a) - 1
    one = 1 << prec
    zs = list(zs)
    ash = [c << prec for c in a]
    done_rounds = 0
    for it in range(maxit):
        mx = 0
        for k in range(n):
            x, y = zs[k]
            px, py = ash[n], 0
            dx, dy = 0, 0
            for j in range(n - 1, -1, -1):
                dx, dy = ((dx * x - dy * y) >> prec) + px, ((dx * y + dy * x) >> prec) + py
                px, py = ((px * x - py * y) >> prec) + ash[j], (px * y + py * x) >> prec
            if px == 0 and py == 0:
                continue
            dd = dx * dx + dy * dy
            if dd == 0:
                zs[k] = (x + (one >> 30), y)
                mx = one
                continue
            Nx = ((px * dx + py * dy) << prec) // dd
            Ny = ((py * dx - px * dy) << prec) // dd
            sx = sy = 0
            for j in range(n):
                if j == k:
                    continue
                ex, ey = x - zs[j][0], y - zs[j][1]
                ee = ex * ex + ey * ey
                if ee == 0:
                    continue
                sx += (ex << (2 * prec)) // ee
                sy -= (ey << (2 * prec)) // ee
            tx = one - ((Nx * sx - Ny * sy) >> prec)
            ty = -((Nx * sy + Ny * sx) >> prec)
            tt = tx * tx + ty * ty
            if tt == 0:
                wx, wy = Nx, Ny
            else:
                wx = ((Nx * tx + Ny * ty) << prec) // tt
                wy = ((Ny * tx - Nx * ty) << prec) // tt
            zs[k] = (x - wx, y - wy)
            mx = max(mx, abs(wx) + abs(wy))
        if mx <= tol_ulps:
            done_rounds += 1
            if done_rounds >= 2:
                return zs, it + 1
        else:
            done_rounds = 0
    return zs, None


def magnitude_bits(a, z):
    """Bits of a bound on the Horner partial sums of a at the points z (floats; chooses precision only)."""
    R = max([1.0] + [abs(x) for x in z]) * 1.25
    lr = math.log2(R)
    return int(max(abs(c).bit_length() + j * lr for j, c in enumerate(a) if c)) + (len(a)).bit_length() + 4


def approx_roots(a, final_bits, start=None):
    """Our own Aberth iteration: floats for a start, then fixed point with Mb + 64 bits until every point moves less
    than 2^-40, then with Mb + final_bits bits until every point moves less than 2^-(final_bits - 40); Mb bounds the
    size of Horner's partial sums, so the attainable accuracy is about 2^-(final_bits). Returns (zs, prec, sweeps)."""
    z0 = float_aberth(a, z=start)
    Mb = magnitude_bits(a, z0)
    p1 = Mb + 64
    zs, its1 = fx_aberth(a, [_to_fixed(z, p1) for z in z0], p1, 1 << (p1 - 40), maxit=600)
    if its1 is None:
        return None, None, None
    p2 = Mb + final_bits
    zs = [(x << (p2 - p1), y << (p2 - p1)) for x, y in zs]
    zs, its2 = fx_aberth(a, zs, p2, 1 << (p2 - final_bits + 40), maxit=30)
    if its2 is None:
        return None, None, None
    return zs, p2, its1 + its2


def gershgorin(a, zs, prec):
    """Certify all roots of the integer polynomial a (ascending) from distinct approximations zs (integer pairs at
    scale 2^-prec). Returns a list of discs (cx, cy, rad) at scale 2^-prec, pairwise disjoint, each containing
    exactly one root; raises ValueError otherwise."""
    n = len(a) - 1
    if len(zs) != n:
        raise ValueError("%d approximations for a polynomial of degree %d" % (len(zs), n))
    if len(set(zs)) != n:
        raise ValueError("approximations not distinct")
    one = 1 << prec
    lc = a[-1]
    ash = [c << prec for c in a]
    discs = []
    for i in range(n):
        x, y = zs[i]
        mz = abs(x) + abs(y)                         # >= |z| 2^prec
        px, py, pr = ash[n], 0, 0                    # Horner in balls: z exact, coefficients exact
        for j in range(n - 1, -1, -1):
            nx = px * x - py * y
            ny = px * y + py * x
            px = (nx >> prec) + ash[j]
            py = ny >> prec
            pr = cdiv(pr * mz, one) + 2
        qx, qy, qr = lc << prec, 0, 0                # lc prod_{j != i} (z_i - z_j)
        for j in range(n):
            if j == i:
                continue
            dx, dy = x - zs[j][0], y - zs[j][1]
            md = abs(dx) + abs(dy)
            nx = qx * dx - qy * dy
            ny = qx * dy + qy * dx
            qx, qy = nx >> prec, ny >> prec
            qr = cdiv(qr * md, one) + 2
        mq = isqrt(qx * qx + qy * qy)                # |q| >= mq (scale 2^-prec)
        if mq <= qr:
            raise ValueError("root %d: the product of differences is not separated from 0" % i)
        qq = qx * qx + qy * qy
        wx = ((px * qx + py * qy) << prec) // qq
        wy = ((py * qx - px * qy) << prec) // qq
        mw = isqrt(wx * wx + wy * wy) + 1
        wr = cdiv((pr << prec) + (mw + 2) * qr, mq - qr) + 2
        rad = wr + (n - 1) * (mw + wr)
        discs.append((x - wx, y - wy, rad))
    for i in range(n):
        xi, yi, ri = discs[i]
        for j in range(i):
            xj, yj, rj = discs[j]
            dx, dy = xi - xj, yi - yj
            if dx * dx + dy * dy <= (ri + rj) ** 2:
                raise ValueError("discs %d and %d overlap" % (j, i))
    return discs


def certified_roots(f_asc):
    """(discs at scale 2^-P, stats) for all roots of f(z(1 - z)), f an integer polynomial (ascending). The roots
    w of f are approximated first; each gives the starting points (1 +- sqrt(1 - 4w))/2 for f o T."""
    t0 = time.time()
    a = compose_T(f_asc)
    ws, pw, its0 = approx_roots(f_asc, 96)
    zs = None
    if ws is not None:
        start = []
        for wf in ws:
            w = _to_float(wf, pw)
            r = cmath.sqrt(1 - 4 * w)
            start += [(1 + r) / 2, (1 - r) / 2]
        zs, prec, its = approx_roots(a, ROOT_BITS, start)
    if zs is None:
        zs, prec, its = approx_roots(a, ROOT_BITS)
    if zs is None:
        raise ValueError("Aberth iteration did not converge")
    discs = gershgorin(a, zs, prec)
    d = prec - P
    out = [(cx >> d, cy >> d, cdiv(r, 1 << d) + 2) for cx, cy, r in discs]
    maxrad = max(r for _, _, r in discs)
    return out, {"deg": len(a) - 1, "iters": its, "prec": prec, "maxrad_log2": maxrad.bit_length() - prec,
                 "time": time.time() - t0}


def _roots_job(args):
    idx, a = args
    try:
        return idx, certified_roots(a), None
    except Exception as exc:          # reported as a failed condition
        return idx, None, str(exc)


# ----------------------------------------------------------------------------------------------------------------
# Reading a certificate (as data) and the exact conditions H1-H3
# ----------------------------------------------------------------------------------------------------------------

RAT = re.compile(r"([0-9]+(/[0-9]+)?|[0-9]+\.[0-9]+)$")       # "a/b" or an exact decimal "x.y"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _intlist(v):
    if not (isinstance(v, list) and v and all(isinstance(c, int) and not isinstance(c, bool) for c in v)):
        raise ValueError("a coefficient array is not a list of integers")
    return list(v)


def load_cert(path):
    """{'den': [(name, ascending coeffs, m)], 'num': [(name, ascending coeffs, q)], 'raw': dict} from a JSON file in
    the [Pin26] format (also used by [Epi26], PR #280 and PR #301): descending integer coefficient arrays,
    'multiplicity' integers, 'q_rational' strings."""
    with open(path) as fh:
        d = json.load(fh)
    den, num = [], []
    for i, f in enumerate(d["denominator"]):
        m = f["multiplicity"]
        if not isinstance(m, int) or isinstance(m, bool):
            raise ValueError("multiplicity is not an integer")
        den.append((f.get("name", "R%d" % (i + 1)), list(reversed(_intlist(f["coeffs"]))), m))
    for i, f in enumerate(d["numerator"]):
        q = f["q_rational"]
        if not (isinstance(q, str) and RAT.match(q)):
            raise ValueError("exponent %r is not a non-negative exact rational" % (q,))
        num.append((f.get("name", "S%d" % (i + 1)), list(reversed(_intlist(f["coeffs"]))), Fraction(q)))
    return {"den": den, "num": num, "raw": d}


def check_hypotheses(rep, cert):
    den, num = cert["den"], cert["num"]
    ok = True
    for name, a, m in den:
        ok &= rep.check("H1 %s: degree >= 1, leading coefficient +-1, multiplicity >= 1" % name,
                        len(a) >= 2 and a[-1] in (1, -1) and m >= 1, "deg %d, lc %d, m %d" % (len(a) - 1, a[-1], m))
    for name, a, q in num:
        ok &= rep.check("H1 %s: degree >= 1, nonzero leading coefficient, q >= 0" % name,
                        len(a) >= 2 and a[-1] != 0 and q >= 0, "deg %d, lc %d, q = %s" % (len(a) - 1, a[-1], q))
    D = sum(m * (len(a) - 1) for _, a, m in den)
    degA = sum(q * (len(a) - 1) for _, a, q in num)
    ok &= rep.check("H2 deg A < D", degA < D, "D = %d, deg A = %s = %.12f D" % (D, degA, float(degA / D)))
    bad, nprs, used = [], 0, set()
    for ni, (nn, s, q) in enumerate(num):
        if q == 0:
            continue
        for di, (dn, r, m) in enumerate(den):
            nprs += 1
            cop, p = coprime_modp(s, r)
            if cop:
                used.add(p)
            else:
                g = gcd_Q(s, r)
                bad.append("%s/%s (exact gcd of degree %d)" % (nn, dn, len(g) - 1))
    ok &= rep.check("H3 gcd(S_i, R_j) = 1 for all %d pairs with q_i > 0" % nprs, not bad,
                    ("modulo prime(s) %s" % sorted(used)) if not bad else "; ".join(bad))
    return ok, D, degA


# ----------------------------------------------------------------------------------------------------------------
# E1, E2: certified roots and the Mahler measure of Q o T
# ----------------------------------------------------------------------------------------------------------------

def all_roots(cert, workers):
    """Certified root discs (scale 2^-P) for every factor's composition with T, in the order den + num."""
    polys = [a for _, a, _ in cert["den"]] + [a for _, a, _ in cert["num"]]
    jobs = list(enumerate(polys))
    if workers <= 1:
        res = [_roots_job(j) for j in jobs]
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers) as pool:
            res = list(pool.imap_unordered(_roots_job, sorted(jobs, key=lambda j: -len(j[1])), chunksize=1))
    res.sort(key=lambda r: r[0])
    return [(out[0] if out else None, out[1] if out else None, err) for _, out, err in res]


def root_abs(x, y, r):
    """(lo, hi) at scale 2^-P with lo <= |z| <= hi for the disc (x, y, r)."""
    s = isqrt(x * x + y * y)
    return max(0, s - r), s + 1 + r


def logM_composed(a_comp, discs):
    """ln M(F o T) = ln|lc| + sum ln max(1, |r|) as an interval at scale 2^-P."""
    lo_m, hi_m, ex = 1, 1, 0
    for x, y, r in discs:
        lo, hi = root_abs(x, y, r)
        lo_m *= max(lo, ONE)
        hi_m *= max(hi, ONE)
        ex -= P
        if lo_m.bit_length() > 600:
            sh = lo_m.bit_length() - 400
            lo_m >>= sh
            hi_m = cdiv(hi_m, 1 << sh)
            ex += sh
    L = ln_pos(lo_m, hi_m, ex) if discs else (0, 0)
    c = ln_pos(abs(a_comp[-1]), abs(a_comp[-1]), 0)
    return L[0] + c[0], L[1] + c[1]


def vieta_check(a_comp, discs):
    """Sum of the root discs contains -a_{n-1}/a_n, and sum ln|r| contains ln|a_0/a_n| (a_0 != 0): self-tests."""
    n = len(a_comp) - 1
    sx = sum(x for x, _, _ in discs)
    sy = sum(y for _, y, _ in discs)
    sr = sum(r for _, _, r in discs)
    tx = Fraction(-a_comp[n - 1], a_comp[n]) * ONE
    ok1 = abs(sx - tx) <= sr + n and abs(sy) <= sr + n
    if a_comp[0] == 0:
        return ok1
    lo_m, hi_m, ex = 1, 1, 0
    for x, y, r in discs:
        lo, hi = root_abs(x, y, r)
        if lo == 0:
            return ok1
        lo_m *= lo
        hi_m *= hi
        ex -= P
        if lo_m.bit_length() > 600:
            sh = lo_m.bit_length() - 400
            lo_m >>= sh
            hi_m = cdiv(hi_m, 1 << sh)
            ex += sh
    L = ln_pos(lo_m, hi_m, ex)
    t = Fraction(abs(a_comp[0]), abs(a_comp[n]))
    T = ln_pos(t.numerator, t.numerator, 0)
    U = ln_pos(t.denominator, t.denominator, 0)
    return ok1 and L[0] <= T[1] - U[0] and T[0] - U[1] <= L[1]


# ----------------------------------------------------------------------------------------------------------------
# E3: the integral of max(0, g) on [0, 1/2]
# ----------------------------------------------------------------------------------------------------------------

L0 = 9                          # first level: panels of width 2^-9
LMAX = 44                       # depth limit
KMAX = 40                       # largest per-root truncation order
NJ = 44                         # terms of the j_k series
R_TOL = ONE >> 16               # a positive panel is accepted when its series remainder is below this (in g units)
CROSS_TOL = ONE >> 40           # an undecided panel is accepted when 2h (upper - lower) is below this

_STIR = {}


def _c_coef(k, m):
    """Coefficient of x^m in (e^x - 1)^k, i.e. k! S2(m, k) / m!, exactly."""
    if (k, m) not in _STIR:
        S = [[0] * (k + 1) for _ in range(m + 1)]
        S[0][0] = 1
        for i in range(1, m + 1):
            for j in range(1, min(i, k) + 1):
                S[i][j] = j * S[i - 1][j] + S[i - 1][j - 1]
        for i in range(k, m + 1):
            _STIR[(k, i)] = Fraction(math.factorial(k) * S[i][k], math.factorial(i))
    return _STIR[(k, m)]


def j_table(L, kmax=KMAX):
    """j_k = (1/(2h e^k)) int_{-h}^{h} (e^{2 pi i s} - 1)^k ds for k = 0..kmax, h = 2^-(L+1), e = EM 2^-(L+16), as
    intervals at scale 2^-P. With theta = pi 2^-L: j_k = (theta/e)^k sum_{n >= 0, k+n even} c_{k,k+n} (-1)^{(k+n)/2}
    theta^n / (k+n+1), and c_{k,k+n} <= k^n / n! bounds the tail."""
    PJ = P + 64
    one = 1 << PJ
    tl = (PI_FR[0].numerator << PJ) // (PI_FR[0].denominator << L)
    th = cdiv(PI_FR[1].numerator << PJ, PI_FR[1].denominator << L)
    rl = (PI_FR[0].numerator << (PJ + 16)) // (PI_FR[0].denominator * EM)     # theta / e
    rh = cdiv(PI_FR[1].numerator << (PJ + 16), PI_FR[1].denominator * EM)
    assert rh <= one
    tpow = [(one, one)]
    for n in range(1, NJ + 2):
        a, b = tpow[-1]
        tpow.append(((a * tl) >> PJ, cdiv(b * th, one)))
    out = [(ONE, ONE)]
    rk = (one, one)
    thf = Fraction(th, one)
    for k in range(1, kmax + 1):
        rk = ((rk[0] * rl) >> PJ, cdiv(rk[1] * rh, one))
        slo = shi = 0
        for n in range(0, NJ + 1):
            if (k + n) % 2:
                continue
            c = _c_coef(k, k + n) / (k + n + 1)
            if ((k + n) // 2) % 2:
                c = -c
            lo, hi = ivq(tpow[n], c.numerator, c.denominator)
            slo += lo
            shi += hi
        x = k * thf
        assert x < 1
        tail = x ** (NJ + 1) / math.factorial(NJ + 1) / (1 - x)
        tail_i = cdiv(tail.numerator << PJ, tail.denominator) + NJ + 2
        slo -= tail_i
        shi += tail_i
        c4 = (slo * rk[0], slo * rk[1], shi * rk[0], shi * rk[1])
        lo, hi = min(c4) >> PJ, cdiv(max(c4), one)
        out.append((lo >> 64, cdiv(hi, 1 << 64)))
    return out


# worker globals (set before the pool forks)
G_ROOTS = []        # per factor: list of discs (x, y, r) at scale 2^-P
G_W = []            # per factor: (p, q) with weight p/q (negative for denominators)
G_AW = []           # per factor: |weight| as a float (chooses truncation orders only)
G_C = (0, 0)        # sum_i q_i ln|lc S_i|
G_JT = {}           # level -> j table
G_MODE = "pos"      # "pos": int max(0, g);  "mean": int g (self-test, every panel must resolve)
G_TAU = 1e-12       # per-root remainder target (float; chooses truncation orders only)


def setup_integrand(factors, roots, C):
    """factors: list of (p, q) weights; roots: list of disc lists; C: interval. Factors of weight 0 are dropped."""
    global G_ROOTS, G_W, G_AW, G_C, G_JT, G_TAU
    G_ROOTS, G_W, G_AW = [], [], []
    for (p, q), rts in zip(factors, roots):
        if p == 0:
            continue
        G_ROOTS.append(rts)
        G_W.append((p, q))
        G_AW.append(abs(p) / q)
    G_C = C
    for L in range(L0, LMAX + 1):
        if L not in G_JT:
            G_JT[L] = j_table(L)
    nroots = sum(len(r) for r in G_ROOTS)
    G_TAU = (R_TOL / ONE) / (8.0 * nroots)


def _choose_K(aw, rho):
    """Smallest K with aw rho^(K+1) / ((K+1)(1 - rho)) <= G_TAU (floats; the remainder is then bounded exactly)."""
    if rho <= 0.0:
        return 1
    K = 1
    val = aw * rho * rho / (2.0 * (1.0 - rho))
    while val > G_TAU and K < KMAX:
        K += 1
        val *= rho * K / (K + 1)
    return K


def panel(L, k):
    """The panel [k 2^-L, (k+1) 2^-L]. Returns (kind, lo, hi): kind 'neg' / 'pos' / 'crude' with [lo, hi] bounding
    the panel MEAN of max(0, g) (or of g in mode 'mean') at scale 2^-P; or ('split', 0, 0); or ('fail', 0, 0)."""
    ux, uy, ur = unit_point(2 * k + 1, 1 << (L + 1))
    e = EM << (P - L - 16)                  # e >= 2 pi h = pi 2^-L, scale 2^-P
    e4 = e << 2
    Glo, Ghi = G_C
    V0 = 0
    near = []
    far = []
    for F in range(len(G_ROOTS)):
        p, q = G_W[F]
        plo = phi = 1
        ex = 0
        T = 0
        fl = []
        for (x, y, rr) in G_ROOTS[F]:
            dx = ux - x
            dy = uy - y
            D2 = dx * dx + dy * dy
            s = isqrt(D2)
            dlt = ur + rr
            dl = s - dlt                        # <= |u_c - r| (true points)
            if dl >= e4:                        # far: e / |u_c - r| <= 1/4
                plo *= dl
                phi *= s + 1 + dlt
                ex -= P
                if plo.bit_length() > 512:
                    sh = plo.bit_length() - 256
                    plo >>= sh
                    phi = cdiv(phi, 1 << sh)
                    ex += sh
                T += cdiv(e << P, dl)
                fl.append((dx, dy, D2, dl, dlt))
            else:
                near.append((F, dl, s + 1 + dlt))
        far.append(fl)
        if fl:
            g0, g1 = ivq(ln_pos(plo, phi, ex), p, q)
            Glo += g0
            Ghi += g1
            V0 += cdiv(abs(p) * T, q)
    V0 = cdiv(4 * V0, 3)                    # |ln|1 + b eta|| <= -ln(1 - |b|) <= (4/3)|b| for |b| <= 1/4
    nh = nl = 0
    nh_inf = nl_inf = False
    for (F, dl, du) in near:
        p, q = G_W[F]
        up = ln_fix(du + e)[1]                  # ln|u - r| <= ln(du + e) on the panel
        lo = ln_fix(dl - e)[0] if dl - e > 0 else None
        if p > 0:
            nh += cdiv(up * p, q)
            if lo is None:
                nl_inf = True
            else:
                nl += (lo * p) // q
        else:
            if lo is None:
                nh_inf = True
            else:
                nh += cdiv(lo * p, q)
            nl += (up * p) // q
    if G_MODE == "pos" and not nh_inf and Ghi + V0 + nh <= 0:
        return "neg", 0, 0
    if nh_inf or (G_MODE == "mean" and near):
        return ("split", 0, 0) if L < LMAX else ("fail", 0, 0)
    # the series
    Sx = [0] * (KMAX + 1)
    Sy = [0] * (KMAX + 1)
    RAD = REM = 0
    KM = 1
    for F in range(len(G_ROOTS)):
        fl = far[F]
        if not fl:
            continue
        p, q = G_W[F]
        aw = G_AW[F]
        sx = [0] * (KMAX + 1)
        sy = [0] * (KMAX + 1)
        rad = rem = 0
        kmF = 1
        for (dx, dy, D2, dl, dlt) in fl:
            nx = ux * dx + uy * dy
            ny = uy * dx - ux * dy
            bx = (e * nx) // D2                 # b = e u_c / (u_c - r) at scale 2^-P (rounded down)
            by = (e * ny) // D2
            rb = cdiv(e * (ur + cdiv(dlt << (P + 1), dl)), dl) + 2      # |b_true - b| (ulps)
            if rb > (ONE >> 3):
                return ("split", 0, 0) if L < LMAX else ("fail", 0, 0)
            rho = cdiv(e << P, dl)              # >= |b_true|
            K = _choose_K(aw, rho / 3.402823669209385e38)
            px, py = bx, by
            sx[1] += px
            sy[1] += py
            for j in range(2, K + 1):
                px, py = (px * bx - py * by) >> P, (px * by + py * bx) >> P
                sx[j] += px
                sy[j] += py
            rad += 2 * rb + 4                   # |b_true^j - computed power| <= 2 rb + 4 ulps for every j
            rp = rho
            for j in range(K):
                rp = cdiv(rp * rho, ONE)
            rem += cdiv(4 * rp, 3 * (K + 1))    # sum_{j > K} |b eta|^j / j <= rho^(K+1) / ((K+1)(1 - rho))
            if K > kmF:
                kmF = K
        for j in range(1, kmF + 1):
            Sx[j] += (sx[j] * p) // q
            Sy[j] += (sy[j] * p) // q
        RAD += cdiv(rad * abs(p), q) + 2
        REM += cdiv(rem * abs(p), q) + 1
        if kmF > KM:
            KM = kmF
    # variation bound: |Re sum_j (-1)^{j+1} S_j eta^j / j| <= |Im S_1| + |S_1| e/2 + sum_{j>=2} |S_j|/j
    V = abs(Sy[1]) + RAD + cdiv((isqrt(Sx[1] ** 2 + Sy[1] ** 2) + 1 + RAD) * e, 2 * ONE) + REM
    for j in range(2, KM + 1):
        V += cdiv(isqrt(Sx[j] ** 2 + Sy[j] ** 2) + 1 + RAD, j)
    sup = Ghi + V + nh
    mu_lo = mu_hi = None
    if not near:
        JT = G_JT[L]
        mu_lo, mu_hi = Glo - REM, Ghi + REM
        for j in range(1, KM + 1):
            lo, hi = ivmul((Sx[j] - RAD, Sx[j] + RAD), JT[j])
            lo, hi = lo // j, cdiv(hi, j)
            if j % 2 == 0:
                lo, hi = -hi, -lo
            mu_lo += lo
            mu_hi += hi
    if G_MODE == "mean":
        return ("pos", mu_lo, mu_hi) if (REM <= R_TOL or L >= LMAX) else ("split", 0, 0)
    if sup <= 0:
        return "neg", 0, 0
    if not near and Glo - V > 0:
        if REM <= R_TOL or L >= LMAX:
            return "pos", mu_lo, mu_hi
        return "split", 0, 0
    lo = max(0, mu_lo) if mu_lo is not None else 0
    if L >= LMAX or ((sup - lo) >> L) <= CROSS_TOL:
        return "crude", lo, sup
    return "split", 0, 0


def _integrate_task(k0):
    stack = [(L0, k0)]
    lo_sum = hi_sum = 0
    counts = {"neg": 0, "pos": 0, "crude": 0, "split": 0}
    maxL, fails, crude_w = L0, [], 0
    while stack:
        L, k = stack.pop()
        kind, lo, hi = panel(L, k)
        if kind == "fail":
            fails.append((L, k))
            continue
        counts[kind] += 1
        maxL = max(maxL, L)
        if kind == "split":
            stack.append((L + 1, 2 * k + 1))
            stack.append((L + 1, 2 * k))
            continue
        if G_MODE == "pos":
            lo = max(0, lo)
        lo_sum += lo << (LS - L)
        hi_sum += hi << (LS - L)
        if kind == "crude":
            crude_w += (hi - lo) << (LS - L)
    return k0, lo_sum, hi_sum, counts, maxL, fails, crude_w


def integrate(workers, mode="pos"):
    """[lo, hi] at scale 2^-(P + LS) for int_0^{1/2} max(0, g) (mode 'pos') or int_0^{1/2} g (mode 'mean')."""
    global G_MODE
    G_MODE = mode
    tasks = list(range(1 << (L0 - 1)))
    t0 = time.time()
    if workers <= 1:
        res = [_integrate_task(k) for k in tasks]
    else:
        import multiprocessing as mp
        order = sorted(tasks, key=lambda i: (i * 7919) % len(tasks))
        with mp.get_context("fork").Pool(workers) as pool:
            res = list(pool.imap_unordered(_integrate_task, order, chunksize=1))
    lo = sum(r[1] for r in res)
    hi = sum(r[2] for r in res)
    counts = {key: sum(r[3][key] for r in res) for key in ("neg", "pos", "crude", "split")}
    fails = [f for r in res for f in r[5]]
    return {"lo": lo, "hi": hi, "counts": counts, "maxL": max(r[4] for r in res), "fails": fails,
            "crude": sum(r[6] for r in res), "time": time.time() - t0}


# ----------------------------------------------------------------------------------------------------------------
# Independent evaluation of g at a point from the exact coefficients (self-test of the root-based model)
# ----------------------------------------------------------------------------------------------------------------

def ball_mul(a, b):
    ax, ay, ar = a
    bx, by, br = b
    x = (ax * bx - ay * by) >> P
    y = (ax * by + ay * bx) >> P
    ma = abs(ax) + abs(ay)
    mb = abs(bx) + abs(by)
    return x, y, cdiv(ma * br + mb * ar + ar * br, ONE) + 2


def g_horner(num, den, factors_asc, weights):
    """g at t = num/den from Horner's rule on the exact coefficients of each F at w = T(u) (no roots used; the
    leading coefficients are inside F)."""
    u = unit_point(num, den)
    u2 = ball_mul(u, u)
    w = (u[0] - u2[0], u[1] - u2[1], u[2] + u2[2])           # T(u) = u - u^2
    lo = hi = 0
    for a, (p, q) in zip(factors_asc, weights):
        if p == 0:
            continue
        acc = (a[-1] << P, 0, 0)
        for c in reversed(a[:-1]):
            m = ball_mul(acc, w)
            acc = (m[0] + (c << P), m[1], m[2])
        s = isqrt(acc[0] ** 2 + acc[1] ** 2)
        if s - acc[2] <= 0:
            return None
        L = ln_pos(s - acc[2], s + 1 + acc[2], -P)
        g0, g1 = ivq(L, p, q)
        lo += g0
        hi += g1
    return lo, hi


def g_roots_point(num, den):
    """g at t = num/den from the certified roots (the integrator's representation), all roots counted."""
    ux, uy, ur = unit_point(num, den)
    lo, hi = G_C
    for F in range(len(G_ROOTS)):
        p, q = G_W[F]
        plo = phi = 1
        ex = 0
        for (x, y, rr) in G_ROOTS[F]:
            dx, dy = ux - x, uy - y
            s = isqrt(dx * dx + dy * dy)
            dl = s - ur - rr
            if dl <= 0:
                return None
            plo *= dl
            phi *= s + 1 + ur + rr
            ex -= P
        g0, g1 = ivq(ln_pos(plo, phi, ex), p, q)
        lo += g0
        hi += g1
    return lo, hi


# ----------------------------------------------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------------------------------------------

class Report:
    def __init__(self, title):
        self.fails = []
        self.t0 = time.time()
        print(title, flush=True)

    def info(self, key, val):
        print("  %s: %s" % (key, val), flush=True)

    def check(self, desc, ok, detail=""):
        print("  check %s: %s%s" % (desc, "PASS" if ok else "FAIL", (" (" + detail + ")") if detail else ""),
              flush=True)
        if not ok:
            self.fails.append(desc)
        return ok

    def done(self):
        print("  run time: %.1f s" % (time.time() - self.t0))
        print("  RESULT: %s" % ("ALL CHECKS PASS" if not self.fails else "FAILED: " + "; ".join(self.fails)))
        return not self.fails


def dec_floor(fr, digits):
    q = (fr.numerator * 10 ** digits) // fr.denominator
    s = str(abs(q)).rjust(digits + 1, "0")
    return ("-" if q < 0 else "") + s[:-digits] + "." + s[-digits:]


def dec_ceil(fr, digits):
    q = -((-fr.numerator * 10 ** digits) // fr.denominator)
    s = str(abs(q)).rjust(digits + 1, "0")
    return ("-" if q < 0 else "") + s[:-digits] + "." + s[-digits:]


def encl(lo, hi, digits=13):
    return "[%s, %s]" % (dec_floor(lo, digits), dec_ceil(hi, digits))


# ----------------------------------------------------------------------------------------------------------------
# One certificate end to end
# ----------------------------------------------------------------------------------------------------------------

def evaluate(rep, cert, workers, tag="", verbose=True):
    """H1-H3, E1-E3 for a loaded certificate; returns a dict with D, the enclosures (Fractions) and timings, or
    None when a condition fails."""
    ok, D, degA = check_hypotheses(rep, cert) if verbose else _quiet_hyp(cert)
    if not ok:
        return None
    t0 = time.time()
    res = all_roots(cert, workers)
    names = [n for n, _, _ in cert["den"]] + [n for n, _, _ in cert["num"]]
    polys = [a for _, a, _ in cert["den"]] + [a for _, a, _ in cert["num"]]
    bad = [(names[i], err) for i, (d, st, err) in enumerate(res) if err]
    nroots = sum(len(d) for d, st, err in res if d)
    vieta = all(vieta_check(compose_T(polys[i]), d) for i, (d, st, err) in enumerate(res) if d)
    if verbose:
        rep.check("E1 every root of every F o T enclosed (%d factors, %d discs, pairwise disjoint per factor)"
                  % (len(polys), nroots), not bad,
                  ("largest disc radius 2^%d; %.1f s" % (max(st["maxrad_log2"] for d, st, err in res if st),
                                                       time.time() - t0)) if not bad else str(bad))
        rep.check("E1 self-test: root sums match Vieta and sum ln|r| encloses ln|a_0/a_n|", vieta)
    if bad or not vieta:
        return None
    discs = [d for d, st, err in res]
    nd = len(cert["den"])
    lm = [0, 0]
    for j, (name, a, m) in enumerate(cert["den"]):
        L = logM_composed(compose_T(a), discs[j])
        lm[0] += m * L[0]
        lm[1] += m * L[1]
    logM = (Fraction(lm[0], ONE), Fraction(lm[1], ONE))
    if verbose:
        rep.info("E2 log M(Q o T)", encl(*logM, digits=10) + "  (/D: " + encl(logM[0] / D, logM[1] / D) + ")")
    weights = [(-m, 1) for _, _, m in cert["den"]] + [(q.numerator, q.denominator) for _, _, q in cert["num"]]
    C = [0, 0]
    for name, a, q in cert["num"]:
        if q:
            c0, c1 = ivq(ln_pos(abs(a[-1]), abs(a[-1]), 0), q.numerator, q.denominator)
            C[0] += c0
            C[1] += c1
    setup_integrand(weights, discs, tuple(C))
    if verbose:
        rng = random.Random(82)
        agree, cmp = True, 0
        for _ in range(12):
            L = rng.randrange(10, 30)
            kk = rng.randrange(1 << (L - 1))
            a1 = g_roots_point(2 * kk + 1, 1 << (L + 1))
            a2 = g_horner(2 * kk + 1, 1 << (L + 1), polys, weights)
            if a1 and a2:
                cmp += 1
                if a1[1] < a2[0] or a2[1] < a1[0]:
                    agree = False
        rep.check("E3 self-test: g from the certified roots overlaps g from Horner on the exact coefficients "
                  "(%d of 12 random points compared)" % cmp, agree and cmp >= 10)
    it = integrate(workers, "pos")
    if it["fails"]:
        rep.check("E3 every panel resolved", False, "%d panels unresolved at depth %d" % (len(it["fails"]), LMAX))
        return None
    I = (Fraction(2 * it["lo"], 1 << (P + LS)), Fraction(2 * it["hi"], 1 << (P + LS)))
    if verbose:
        c = it["counts"]
        rep.info("E3 int_0^1 max(0, g) dt", "%s  (%d panels: %d negative, %d positive, %d crude, %d split; depth "
                 "%d..%d; crude width %.3e; %.1f s)" % (encl(*I, digits=8), sum(c.values()) - c["split"], c["neg"],
                                                       c["pos"], c["crude"], c["split"], L0, it["maxL"],
                                                       2 * it["crude"] / 2.0 ** (P + LS), it["time"]))
    B = ((logM[0] + I[0]) / D, (logM[1] + I[1]) / D)
    if verbose:
        jensen_selftests(rep, cert, discs, workers)
    return {"D": D, "degA": degA, "logM": logM, "I": I, "B": B, "it": it, "discs": discs, "weights": weights}


def jensen_selftests(rep, cert, discs, workers):
    """int_0^1 g dt = sum_F w_F ln M(F o T) (Jensen), with the left side from the panel series over every panel and
    the right side from the root moduli: for the denominators, and for the numerators with no root within 2^-12 of
    the unit circle (the others make g unbounded below)."""
    nd = len(cert["den"])
    for label, idx in (("denominators", list(range(nd))),
                       ("numerators away from the circle", [nd + i for i, (n, a, q) in enumerate(cert["num"])
                                                            if q and all(abs(root_abs(*r)[0] - ONE) > (ONE >> 12)
                                                                         and abs(root_abs(*r)[1] - ONE) > (ONE >> 12)
                                                                         for r in discs[nd + i])])):
        if not idx:
            continue
        ws, ds, J = [], [], [Fraction(0), Fraction(0)]
        C = [0, 0]
        for i in idx:
            if i < nd:
                n, a, m = cert["den"][i]
                w = Fraction(-m)
            else:
                n, a, w = cert["num"][i - nd]
                c0, c1 = ivq(ln_pos(abs(a[-1]), abs(a[-1]), 0), w.numerator, w.denominator)
                C[0] += c0
                C[1] += c1
            L = logM_composed(compose_T(a), discs[i])
            lo, hi = ivq(L, w.numerator, w.denominator)
            J[0] += Fraction(lo, ONE)
            J[1] += Fraction(hi, ONE)
            ws.append((w.numerator, w.denominator))
            ds.append(discs[i])
        setup_integrand(ws, ds, tuple(C))
        it = integrate(workers, "mean")
        I = (Fraction(2 * it["lo"], 1 << (P + LS)), Fraction(2 * it["hi"], 1 << (P + LS)))
        rep.check("E3 self-test (Jensen): int g over all panels equals sum w ln M(F o T), %d %s" % (len(idx), label),
                  not it["fails"] and I[0] <= J[1] and J[0] <= I[1],
                  "%s vs %s; %d panels, %.1f s" % (encl(*I, digits=6), encl(*J, digits=6), it["counts"]["pos"],
                                                  it["time"]))


def _quiet_hyp(cert):
    class Q:
        fails = []

        def check(self, d, ok, detail=""):
            if not ok:
                self.fails.append(d)
            return ok

        def info(self, *a):
            pass
    return check_hypotheses(Q(), cert)


# ----------------------------------------------------------------------------------------------------------------
# Self-tests of the arithmetic and of the integrator against closed forms
# ----------------------------------------------------------------------------------------------------------------

def selftests(rep, workers):
    from decimal import Decimal, getcontext
    getcontext().prec = 110
    ok = True
    for n in (2, 3, 10, 12345678901234567890123, 2 ** 200 + 12345):
        lo, hi = ln_pos(n, n, 0)
        d = Decimal(n).ln() * Decimal(2 ** P)
        ok &= Decimal(lo) <= d <= Decimal(hi)
    rep.check("S1 ln n against decimal's correctly rounded ln (5 values)", ok)
    PI100 = Decimal("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679")
    ok = (Decimal(PI_FR[0].numerator) / Decimal(PI_FR[0].denominator) <= PI100 + Decimal(10) ** -99 and
          Decimal(PI_FR[1].numerator) / Decimal(PI_FR[1].denominator) >= PI100)
    rep.check("S1 pi (Machin) brackets 100 digits of pi", ok)

    def dcs(a):
        c = Decimal(0)
        t = Decimal(1)
        k = 0
        while abs(t) > Decimal(10) ** -100:
            c += t
            k += 2
            t = -t * a * a / (k * (k - 1))
        s = Decimal(0)
        t = a
        k = 1
        while abs(t) > Decimal(10) ** -100:
            s += t
            k += 2
            t = -t * a * a / (k * (k - 1))
        return c, s
    ok = True
    for num, den in ((0, 2), (1, 4), (1, 2), (1, 6), (5, 1024), (12345, 2 ** 20), (1, 2 ** 40), (2 ** 30 - 1, 2 ** 31)):
        x, y, r = unit_point(num, den)
        c, s = dcs(2 * PI100 * num / den)
        ex = Decimal(x) / Decimal(2 ** P) - c
        ey = Decimal(y) / Decimal(2 ** P) - s
        ok &= (ex * ex + ey * ey).sqrt() <= Decimal(r) / Decimal(2 ** P)
    rep.check("S1 exp(2 pi i c) balls contain the decimal cos/sin (8 points)", ok)
    discs, st = certified_roots([-1, 0, 0, 0, 0, 0, 0, 0, 1])                  # f = x^8 - 1 composed with T
    a = compose_T([-1, 0, 0, 0, 0, 0, 0, 0, 1])
    ok = len(discs) == 16
    for kk in range(8):                                     # each root z of z(1-z) = w, w^8 = 1, inside one disc
        w = cmath.exp(2j * math.pi * kk / 8)
        for sgn in (1, -1):
            z = (1 + sgn * cmath.sqrt(1 - 4 * w)) / 2
            inside = [i for i, (x, y, r) in enumerate(discs) if abs(complex(x, y) / 2.0 ** P - z) < 1e-12]
            ok &= len(inside) == 1
    rep.check("S2 roots of (z(1-z))^8 - 1: 16 disjoint discs, one per known root", ok, "max radius 2^%d" %
              st["maxrad_log2"])
    a = compose_T([269, -97, 30, -7, 1])
    zs, prec, its = approx_roots(a, ROOT_BITS)
    try:
        gershgorin(a, zs[:-1], prec)
        ok = False
    except ValueError:
        ok = True
    rep.check("S2 red: a missing root approximation is refused", ok)
    a = compose_T([1, -2, 1])                                # (T - 1)^2: two double roots exp(+-i pi/3)
    zs, prec, its = approx_roots(a, ROOT_BITS)
    try:
        ok = zs is None
        if zs is not None:
            gershgorin(a, zs, prec)
    except ValueError:
        ok = True
    rep.check("S2 red: (z(1-z) - 1)^2, with two double roots, is refused", ok)
    # S3: the integrator against closed forms
    f = [-9, 1]                                             # S = x - 9: g = ln|S(T(u))| > 0 on the circle
    d1, _ = certified_roots(f)
    LM = logM_composed(compose_T(f), d1)
    setup_integrand([(1, 1)], [d1], (0, 0))
    it = integrate(workers, "pos")
    I = (Fraction(2 * it["lo"], 1 << (P + LS)), Fraction(2 * it["hi"], 1 << (P + LS)))
    J = (Fraction(LM[0], ONE), Fraction(LM[1], ONE))
    rep.check("S3 int max(0, ln|S(T(u))|) for S = x - 9 equals ln M(S o T) (Jensen)", I[0] <= J[1] and J[0] <= I[1],
              "%s vs %s" % (encl(*I, digits=15), encl(*J, digits=15)))
    f = [1, -19, 178, -1084, 4786, -16221, 43677, -95554, 172556, -260367, 331775, -360911, 339335, -279897,
         206152, -138288, 86189, -50741, 28457, -15180, 7646, -3618, 1612, -679, 269, -97, 30, -7, 1]   # Doche Q1
    d1, _ = certified_roots(f)
    LM = logM_composed(compose_T(f), d1)
    setup_integrand([(-3, 1)], [d1], (0, 0))
    it = integrate(workers, "mean")
    I = (Fraction(2 * it["lo"], 1 << (P + LS)), Fraction(2 * it["hi"], 1 << (P + LS)))
    J = (Fraction(-3 * LM[1], ONE), Fraction(-3 * LM[0], ONE))
    rep.check("S3 int g for g = -3 ln|Q1(T(u))| equals -3 ln M(Q1 o T) (Jensen; every panel by the series)",
              I[0] <= J[1] and J[0] <= I[1], "%s vs %s, %d panels" % (encl(*I, digits=15), encl(*J, digits=15),
                                                                     it["counts"]["pos"]))
    # a crossing case against a float quadrature (sanity, not a decision)
    f2 = [-3, 1]
    d2, _ = certified_roots(f2)
    setup_integrand([(-3, 1), (5, 2)], [d1, d2], (0, 0))
    it = integrate(workers, "pos")
    I = (Fraction(2 * it["lo"], 1 << (P + LS)), Fraction(2 * it["hi"], 1 << (P + LS)))
    N = 200000
    tot = 0.0
    for i in range(N):
        u = cmath.exp(2j * math.pi * (i + 0.5) / N)
        w = u - u * u
        q1 = 0
        for c in reversed(f):
            q1 = q1 * w + c
        tot += max(0.0, 2.5 * math.log(abs(w - 3)) - 3 * math.log(abs(q1)))
    tot /= N
    rep.check("S3 a sign-changing g (5/2 ln|T - 3| - 3 ln|Q1 o T|) against a 200000-point midpoint sum (floats)",
              abs(float(I[0]) - tot) < 1e-6 * max(1, tot), "%s vs %.10f; %d crude panels" % (encl(*I, digits=12), tot,
                                                                                       it["counts"]["crude"]))


# ----------------------------------------------------------------------------------------------------------------
# Rows
# ----------------------------------------------------------------------------------------------------------------

ROWS = {
    "pin": {"title": "[Pin26] C_82 < 0.2534449* (82a page; PR #242 asks to remove the asterisk)",
            "file": "pin/certificate.json",
            "sha": "b33626f367b6f52922746c0fee7bb008449cd0bb18efcb983ab20ff073350fce", "sha_src": "82a page",
            "printed": "0.2534449", "D": 999999992,
            "stated": [("README: verify.py certified upper bound", "0.25344468733362658", "upper"),
                       ("README: independent_check.py certified upper bound", "0.25344468733741928", "upper"),
                       ("README: referee's evaluator, below", "0.253444687331849", "upper"),
                       ("82a page: certified upper bound", "0.2534446873", "digits"),
                       ],
            "encl": ("PR #242 [T2026c]", "0.253444687248", "0.253444687515")},
    "epi": {"title": "[Epi26] C_82 <= 0.2535744749068227* (82a page; PR #242 asks to remove the asterisk)",
            "file": "epi/construction.json",
            "sha": "e85e26fa7122a2e68b19b2fa01e27e78bbd802450bce2a0d45717f09b930e559",
            "sha_src": "result.json artifact_sha256", "printed": "0.2535744749068227", "D": None,
            "stated": [("82a page: [Pin26]'s Arb re-evaluation 'encloses it at'", "0.25356893", "digits")],
            "encl": ("PR #242 [T2026c]", "0.253568930187", "0.253568930424")},
    "pr280": {"title": "PR #280 [Pat26] C_82 < 0.2534350402* (open)",
              "file": "pr280/certificate.json",
              "sha": "557da3e3cfdaee57eee84247a4daf9c470f88b0513c300dcb127dadfee775fe0",
              "sha_src": "as fetched at the pinned commit 32898ac (no hash is stated)",
              "printed": "0.2534350402", "D": 1000000000004,
              "stated": [("PR body: vcheck.py tol 1e-10", "0.2534350401297", "upper"),
                         ("PR body: vcheck.py tol 1e-11, 192-bit", "0.2534350401282", "upper"),
                         ("PR body: Pinho's verify.py", "0.2534350401467", "upper"),
                         ("PR body: Pinho's independent_check.py", "0.2534350401334", "upper"),
                         ("PR body: certified upper bound", "0.25343504013", "digits")],
              "encl": None},
    "pr280b": {"title": "PR #280 [Pat26] updated 4 Nov: C_82 < 0.2524177045* (open; tao/certificates/c82/v2 at 3a7a84f)",
               "file": "pr280b/certificate.json",
               "sha": "1125ddf42e9a544c6fc566aa2ebdfb59ec7c35177030c5d50767427e52f841d1",
               "sha_src": "SHA256SUMS of the folder", "printed": "0.2524177045", "D": None,
               "stated": [("README: vcheck.py tol 1e-11, 192-bit", "0.25241770444734", "upper"),
                          ("README: Pinho's verify_exact.py (unmodified), 256-bit", "0.25241770444732", "upper"),
                          ("README: certified value", "0.25241770445", "digits")],
               "encl": None},
    "pr280c": {"title": "PR #280 [Pat26] updated 5 Nov: C_82 < 0.2520591386* (open; tao/certificates/c82/v3 at 88c39d5)",
               "file": "pr280c/certificate.json",
               "sha": "00125f6984ad2646a79d8782ee38f5161bd24c17d4899e4af74ff3e529ba1dbe",
               "sha_src": "SHA256SUMS of the folder", "printed": "0.2520591386", "D": None,
               "stated": [("README: vcheck.py (Arb, 192 bits, --tol 1e-11) enclosure 0.2520591385260044869 +/- 1.67e-26, midpoint", "0.2520591385260044869", "noted")],
               "encl": None},
    "pr301": {"title": "PR #301 [Pin26b] C_82 < 0.253091834* (open; degree44-denominators-v2 at 6d4e80c)",
              "file": "pr301/certificate.json",
              "sha": "adfca0754f4d8b3a64398df22a506932743de7b0e49df542961530cb26c67ac2",
              "sha_src": "README and SHA256SUMS of the folder", "printed": "0.253091834", "D": 1000000000064,
              "stated": [("README: verify_exact.py rigorous upper endpoint", "0.25309183355151", "upper"),
                         ("README: the review/ integrator 'gives'", "0.25309183355146", "digits")],
              "encl": None},
    "pr301a": {"title": "PR #301 first filing, C_82 < 0.253123275 (degree44-denominators at 6d4e80c; superseded)",
               "file": "pr301a/certificate.json",
               "sha": "6a24d45c18012bf7cc3f0ea05f7bcec685a43fd2f7e942f1ca9ce69fd558f056",
               "sha_src": "SHA256SUMS of the folder", "printed": "0.253123275", "D": None,
               "stated": [], "encl": None},
    "pr330": {"title": "PR #330 [Pin26c] C_82 < 0.252404806* (open; rpinho/math-certificates c82/degree104-cvp at 5461d49)",
              "file": "pr330/certificate.json",
              "sha": "27a4e61121b0c882c249fe3f06ed04d3c6b16067b7d1a0f6d2d6a80d40b3d3f9",
              "sha_src": "SHA256SUMS of the folder", "printed": "0.252404806", "D": None,
              "stated": [], "encl": None},
    "doche": {"title": "positive control: Doche 2001b, C_82 <= 0.25444 (data as in [Pin26] controls/doche2001.json)",
              "file": "doche/doche2001.json",
              "sha": "835a10b9ac5454798475a3c171796f87296ff14ed9d979d7bf14ab88d1c37fc2", "sha_src": "as fetched",
              "printed": "0.25444", "D": None, "stated": [], "encl": None},
}


def modified(cert, den=None, num=None):
    return {"den": den if den is not None else list(cert["den"]),
            "num": num if num is not None else list(cert["num"]), "raw": cert["raw"]}


def controls(rep, cert, base, printed, workers, full, root=None):
    B = base["B"]
    lowered = Fraction(int(B[0] * 10 ** 13) - 1, 10 ** 13)          # 1e-13 (or less) below our lower end
    rep.check("R1 red: the bound %s, just below our enclosure, is not implied (B_lo > it)" % dec_floor(lowered, 13),
              B[0] > lowered)
    n0, a0, m0 = cert["den"][0]
    c2 = modified(cert, num=list(cert["num"]) + [("copy_of_" + n0, a0, Fraction(1, 1000))])
    ok, _, _ = _quiet_hyp(c2)
    rep.check("R2 red: %s added as a numerator (q = 1/1000) is caught by H3 (exact gcd of degree %d)" %
              (n0, len(gcd_Q(a0, a0)) - 1), not ok)
    c3 = modified(cert, den=[(n0, [2 * c for c in a0], m0)] + list(cert["den"][1:]))
    rep.check("R3 red: %s doubled (not monic) is refused by H1" % n0, not _quiet_hyp(c3)[0])
    c4 = modified(cert, num=[(n, a, q * Fraction(D_OVER(cert), 1)) for n, a, q in cert["num"]])
    rep.check("R4 red: every q_i scaled so that deg A = D + 1 is refused by H2", not _quiet_hyp(c4)[0])
    if not full:
        return
    # R5: a perturbed numerator coefficient changes B
    i5 = max(range(len(cert["num"])), key=lambda i: cert["num"][i][2] * (len(cert["num"][i][1]) - 1))
    n5, a5, q5 = cert["num"][i5]
    a5p = [a5[0] + 1] + a5[1:]
    num5 = list(cert["num"])
    num5[i5] = (n5 + "+1", a5p, q5)
    sub = Report("    [R5] %s with constant term + 1" % n5)
    r5 = evaluate(sub, modified(cert, num=num5), workers, verbose=False)
    if r5 is None:
        rep.check("R5 red: perturbed factor %s (constant term + 1)" % n5, True, "refused by H1-H3 or E1")
    else:
        rep.check("R5 red: perturbed factor %s (constant term + 1) changes B" % n5,
                  r5["B"][0] > B[1] or r5["B"][1] < B[0], "B becomes %s (%s the printed %s)" %
                  (encl(*r5["B"]), "above" if r5["B"][0] > printed else "still below", dec_floor(printed, 10)))
    # R6: a perturbed multiplicity changes B
    den6 = [(n0, a0, m0 + 1)] + list(cert["den"][1:])
    r6 = evaluate(Report("    [R6] m(%s) + 1" % n0), modified(cert, den=den6), workers, verbose=False)
    rep.check("R6 red: m(%s) + 1 changes B" % n0, r6 is not None and (r6["B"][0] > B[1] or r6["B"][1] < B[0]),
              "B becomes %s" % (encl(*r6["B"], digits=17) if r6 else "?") + ", was " + encl(*B, digits=17))
    # S4: the same certificate on a different tiling (first level 2^-11, a 16x tighter remainder target)
    global L0, R_TOL
    saved = (L0, R_TOL)
    L0, R_TOL = 11, ONE >> 20
    try:
        r8 = evaluate(Report("    [S4] first level 11, remainder target 2^-20"), cert, workers, verbose=False)
    finally:
        L0, R_TOL = saved
    rep.check("S4 the same certificate on another tiling (first level 2^-11, remainder target 2^-20): the "
              "enclosures overlap", r8 is not None and r8["B"][0] <= B[1] and B[0] <= r8["B"][1],
              "B = %s, %d panels" % (encl(*r8["B"], digits=17), sum(r8["it"]["counts"].values())) if r8 else "?")
    # R7: [Pin26]'s own red-control file (the same certificate with m(Q1) doubled), read as data
    mpath = os.path.join(root, "pinm", "mutated_v3_should_fail.json")
    if os.path.exists(mpath):
        h = sha256(mpath)
        r7 = evaluate(Report("    [R7] %s" % mpath), load_cert(mpath), workers, verbose=False)
        rep.check("R7 red: [Pin26] controls/mutated_v3_should_fail.json (sha %s..., m(Q1) doubled) fails the target"
                  " %s" % (h[:12], dec_floor(printed, 7)), r7 is not None and r7["B"][0] > printed,
                  "B = %s" % (encl(*r7["B"]) if r7 else "?"))


def D_OVER(cert):
    D = sum(m * (len(a) - 1) for _, a, m in cert["den"])
    degA = sum(q * (len(a) - 1) for _, a, q in cert["num"])
    return (D + 1) / degA


def run_row(key, root, do_controls, workers):
    row = ROWS[key]
    rep = Report(row["title"])
    path = os.path.join(root, row["file"])
    h = sha256(path)
    rep.check("SHA-256 of %s matches %s" % (row["file"], row["sha_src"]), h == row["sha"], h)
    cert = load_cert(path)
    rep.info("factors", "%d denominators (%s), %d numerators; %d roots after composing with T" % (
        len(cert["den"]), ", ".join("%s^%d" % (n, m) for n, _, m in cert["den"]), len(cert["num"]),
        2 * sum(len(a) - 1 for _, a, _ in cert["den"] + cert["num"])))
    res = evaluate(rep, cert, workers)
    if res is None:
        return rep.done()
    D = res["D"]
    if row["D"] is not None:
        rep.check("D = %d as stated" % row["D"], D == row["D"])
    printed = Fraction(row["printed"])
    B = res["B"]
    rep.info("B", "%s  (width %.2e)" % (encl(*B, digits=16), float(B[1] - B[0])))
    rep.check("E4 B < %s, the printed upper bound (B_hi = %s; margin %.3e)" % (row["printed"], dec_ceil(B[1], 16),
                                                                              float(printed - B[1])), B[1] < printed)
    for what, v, kind in row["stated"]:
        fv = Fraction(v)
        if kind == "upper":
            rep.check("consistent with the stated upper bound: %s %s (our B_lo is not above it)" % (what, v),
                      B[0] <= fv)
        elif kind == "noted":   # a stated value compared and reported, not a condition of the bound
            gap = fv - B[1] if fv > B[1] else (B[0] - fv if fv < B[0] else Fraction(0))
            rep.info("stated value", "%s %s: %s our enclosure by %.2e%s" % (
                what, v, "above" if fv > B[1] else "below" if fv < B[0] else "inside", float(gap),
                "" if gap == 0 else " (not a condition of the bound; reported in the replay)"))
        else:                   # displayed digits: rounded or truncated, so B meets [v - ulp/2, v + ulp)
            ulp = Fraction(1, 10 ** (len(v.split(".")[1])))
            rep.check("consistent with the displayed digits: %s %s (rounded or truncated)" % (what, v),
                      B[0] < fv + ulp and B[1] >= fv - ulp / 2)
    if row["encl"]:
        who, a, b = row["encl"]
        rep.check("inside the enclosure %s [%s, %s]" % (who, a, b), Fraction(a) <= B[0] and B[1] <= Fraction(b))
    if do_controls:
        controls(rep, cert, res, printed, workers, full=(key == "pin"), root=root)
    return rep.done()


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    workers = max(1, (os.cpu_count() or 2) - 4)
    for a in sys.argv:
        if a.startswith("--workers="):
            workers = int(a.split("=")[1])
    if len(args) == 1 and args[0] == "selftest":
        rep = Report("self-tests of the arithmetic, root certification and integrator")
        selftests(rep, workers)
        sys.exit(0 if rep.done() else 1)
    if len(args) != 2 or args[0] not in ROWS:
        print(__doc__)
        sys.exit(2)
    sys.exit(0 if run_row(args[0], args[1], "--controls" in sys.argv, workers) else 1)


if __name__ == "__main__":
    main()
