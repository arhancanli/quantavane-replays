"""Independent check of three C_42 upper bounds that rest on Bakeas's Theorem A (Turan's pure power-sum constant,
constants/42a.md): 0.68637* (Bakeas, PR #302), 0.6863607950* (Pinho, PR #305) and 0.686360751* (Vatopoulos, PR #327).
Written from the mathematics of Bakeas's NOTES.md (Theorem A and Lemma 3) and the three profile files: no checker,
solver or notebook of any of the three repositories (rigorous_limit.py, independent_check.py, verify_volterra*.py,
check_bakeas_*.py, rigorous_check.py, spectral_witness.py, ...) was read or run. We also read, as notes only, Pinho's
README.md, BAKEAS_ADAPTATION.md and VOLTERRA_VERIFIER_MATH.md and Vatopoulos's NOTES_CERTIFIER.md (and the headings
and opening claims of his README.md, NOTES.md and MATH.md). Two devices below also appear in NOTES_CERTIFIER.md, the
U/V identity of its Lemma C1 and the inequality sqrt(X) >= (3/2) X u - (1/2) X^2 u^3; we re-derived both.

  ./fetch.sh && python3 -I thmA_check.py [CERTS_ROOT] [--rows=bak,pin,vat] [--selftest] [--controls] [--write]

The statement checked (Theorem A). With T > 0, knots 0 = u_0 < ... < u_K = 1/2, phases psi_j, phi the piecewise linear
interpolant, g = T e^(i phi), gamma = g(0) = T e^(i psi_0), a = Re gamma, d = gamma - g on (0, 1/2] and d = gamma on
(1/2, 1]: if kappa = sup |d| < 1 and a >= 0, then limsup R_n <= max(T, F), where rho is the bounded solution of
  t rho(t) = d(t) + int_0^t d(u) rho(t - u) du      (0 < t <= 1),
  p(t) = t^-gamma + int_0^t (t - v)^-gamma rho(v) dv,   F = |p(1)| / int_0^(1/2) |p(t)| dt / (1 - t).
We certify an upper bound for F, then compare it with the claim C >= T.

How (all our own; see the README section for the derivations):
 1. Hypotheses in exact rationals and decimal intervals: 0 < T <= C < 1, the knots, a = T cos psi_0 > 0, and
    kappa = max(T, 2 T sin(E/2)) < 1 with E = max |psi_j - psi_0| < 3.
 2. A trial function rho~ (a discontinuous polynomial of degree P on each cell of the uniform grid h = 1/N, which
    contains every knot): floating-point collocation, then rounded to exact dyadic coefficients. Nothing below depends
    on how it was found.
 3. Its residual R = d + d*rho~ - t rho~, divided by t, bounded on every cell: d is replaced by its Taylor polynomial
    of degree Q on each cell (remainder bounded), the convolution of the two piecewise polynomials is computed exactly
    as a two-dimensional convolution of integer arrays (one product of decimal integers, Kronecker substitution, in
    the decimal module's exact arithmetic), and the sup of each residual polynomial is bounded by the sum of the
    moduli of its Chebyshev coefficients, computed exactly.
 4. From the residual to e = rho - rho~: |e(t)| <= r(t) + (1/t) int_0^t |d(u)| |e(t - u)| du, so any bounded E >= 0 with
    E >= r + (1/t) int |d| E bounds |e| (the operator has norm <= kappa < 1). We build such an E, constant on blocks
    of cells, by a forward recursion, and integrate it against the two kernels of F.
 5. p~(1) by the binomial series of (c - y)^-gamma on every cell (tails bounded), and the denominator of F from below:
    on each cell p~ = h^-gamma (A(xi) + S), S = h x^(1-gamma) J(x), where A is a polynomial (from a second exact
    two-dimensional convolution for the far cells, explicit series for the rest) and J holds the jumps of rho~ at the
    left node; then |A| >= (3/2) X u - (1/2) X^2 u^3 with X = |A|^2 and any u >= 0 (equality at u = 1/|A|, second
    order), integrated exactly, and |A + S| >= |A| + Re(S conj A) u - eta |S| with eta >= sup |u^2 X - 1|; at the cells
    that start at a knot, the first inequality applied to X = |A + S|^2 in generalized monomials
    x^(k + m (1 - gamma) + m' (1 - conj gamma)).
 6. F <= (|p~(1)| + EN) / (Den~ - ED) < C. A diagnostic lower bound for F (not part of the claim) comes from upper
    bounds of the same integrals; the red controls use it to show that their claims are false, not only unproved.
Python 3.9+, standard library only."""
import cmath
import decimal
import json
import math
import os
import sys
import time
from decimal import Decimal
from fractions import Fraction
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)            # for c42_check.py next to this file (python3 -I leaves the script's folder out)
from c42_check import I, Z, q, iln, iexp, isincos, LO, HI, NEAR  # noqa: E402

ZERO, ONE = Decimal(0), Decimal(1)
KINK = Decimal("1e-12")             # a cell whose crude jump bound would cost more is integrated with the jump kept
EXACT = decimal.Context(prec=decimal.MAX_PREC, Emax=decimal.MAX_EMAX, Emin=decimal.MIN_EMIN,
                        traps=[decimal.Inexact, decimal.Rounded, decimal.InvalidOperation, decimal.Overflow])


# ---- arithmetic helpers ----------------------------------------------------------------------------------------------
def _pi():
    """pi = 16 arctan(1/5) - 4 arctan(1/239) (Machin), between consecutive partial sums (as in kos_check.py)."""
    def atan_inv(n, terms):
        s = sum(Fraction((-1) ** k, (2 * k + 1) * n ** (2 * k + 1)) for k in range(terms))
        nxt = Fraction(1, (2 * terms + 1) * n ** (2 * terms + 1))
        return (s, s + nxt) if terms % 2 == 0 else (s - nxt, s)
    a, b = atan_inv(5, 60), atan_inv(239, 20)
    return I(q(16 * a[0] - 4 * b[1]).lo, q(16 * a[1] - 4 * b[0]).hi)


TWO_PI = _pi() * I(2)


def isincos_any(x):
    """(sin x, cos x) for any interval x: reduce by a multiple of the rigorous 2 pi, then c42_check.isincos."""
    k = round(float(x.lo) / (2 * math.pi))
    return isincos(x - TWO_PI * I(k) if k else x)


def up(x):
    """A decimal upper bound for the rational x."""
    return q(x).hi


def dn(x):
    return q(x).lo


def zmul(x, y):
    return Z(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re)


def zscale(x, r):
    """x times the rational r."""
    rr = q(r)
    return Z(x.re * rr, x.im * rr)


def zabs_hi(z):
    return NEAR.next_plus(NEAR.sqrt(z.abs2_hi()))


def zinv_z(z):
    """1 / z for a rectangle z away from 0."""
    d = z.re * z.re + z.im * z.im
    return Z(z.re / d, -z.im / d)


def fix(x, W):
    """(n, e): n = round(2^W mid x) and e >= |t - n 2^-W| for every t in the interval x."""
    S = Decimal(1 << W)
    lo, hi = LO.multiply(x.lo, S), HI.multiply(x.hi, S)
    n = int(NEAR.divide(NEAR.add(lo, hi), 2).to_integral_value(rounding=decimal.ROUND_HALF_EVEN))
    e = max(HI.subtract(hi, Decimal(n)), HI.subtract(Decimal(n), lo), ZERO)
    return n, HI.divide(e, S)


def fixz(z, W):
    nr, er = fix(z.re, W)
    ni, ei = fix(z.im, W)
    return nr, ni, HI.add(er, ei)


def dint(n, W=0):
    """An upper bound for |n| 2^-W (n an integer)."""
    return HI.divide(Decimal(abs(n)), Decimal(1 << W))


def dsum(xs):
    s = ZERO
    for x in xs:
        s = HI.add(s, x)
    return s


def frac_lo(x):
    return LO.divide(Decimal(x.numerator), Decimal(x.denominator))


def frac_hi(x):
    return HI.divide(Decimal(x.numerator), Decimal(x.denominator))


# ---- the profile and the hypotheses of Theorem A ---------------------------------------------------------------------
def load_profile(path):
    d = json.load(open(path))
    T = Fraction(d["T"])
    ph = [Fraction(s) for s in d["phases"]]
    if "breaks" in d:
        kn = [Fraction(s) for s in d["breaks"]]
    else:
        assert d.get("knots") == "j/48, j = 0..24", "unknown knot format"
        kn = [Fraction(j, 48) for j in range(25)]
    C = Fraction(d["C"]) if "C" in d else T
    return dict(T=T, phases=ph, knots=kn, C=C)


def hypotheses(prof, N, out=print):
    """The hypotheses of Theorem A (and the grid), rigorously; returns (ok, info)."""
    T, ph, kn, C = prof["T"], prof["phases"], prof["knots"], prof["C"]
    msgs = []
    ok = True
    if not (0 < T <= C < 1):
        ok = False
        msgs.append("0 < T <= C < 1 fails")
    if len(ph) != len(kn) or kn[0] != 0 or kn[-1] != Fraction(1, 2) or any(
            kn[j + 1] <= kn[j] for j in range(len(kn) - 1)):
        return False, dict(msgs=msgs + ["knots must increase from 0 to 1/2, one phase per knot"])
    if N % 2 or any((u * N).denominator != 1 for u in kn):
        return False, dict(msgs=msgs + ["grid N = %d does not contain every knot" % N])
    s0, c0 = isincos_any(q(ph[0]))
    a, b = q(T) * c0, q(T) * s0
    if not a.lo > 0:
        ok = False
        msgs.append("a = Re gamma = T cos psi_0 > 0 not certified: [%s, %s]" % (fmt(a.lo, 6), fmt(a.hi, 6)))
    E = max(abs(p - ph[0]) for p in ph)
    if not E < 3:
        return False, dict(msgs=msgs + ["max |psi_j - psi_0| >= 3"])
    sE, _ = isincos(q(E / 2))
    dmax = q(2 * T) * sE                      # sup |d| on (0, 1/2]: |d| = 2T |sin((phi - psi_0)/2)|, phi between knots
    kappa_hi = max(up(T), dmax.hi)
    if not kappa_hi < 1:
        ok = False
        msgs.append("kappa < 1 not certified")
    slopes = [(ph[j + 1] - ph[j]) / (kn[j + 1] - kn[j]) for j in range(len(kn) - 1)]
    return ok, dict(a=a, b=b, E=E, dmax=dmax, kappa_hi=kappa_hi, slopes=slopes, msgs=msgs)


def fmt(x, digits=12):
    with decimal.localcontext(NEAR):
        return format(+Decimal(x), ".%dg" % digits) if not isinstance(x, str) else x


# ---- d on the grid: its Taylor polynomials of degree Q, rounded, with bounds ---------------------------------------
def dtilde(prof, info, N, Q, W):
    """E[i] = [(re, im) for a <= Q]: 2^W times the rounded Taylor coefficients of delta(i + theta) = d(h (i + theta)),
    theta in [0, 1]; eta[i] >= |delta - delta~| on cell i, and for cell 0 |delta - delta~| <= eta0 theta."""
    T, ph, kn = prof["T"], prof["phases"], prof["knots"]
    H, K = N // 2, len(kn) - 1
    M = [int(u * N) for u in kn]
    gam = Z(q(T) * isincos_any(q(ph[0]))[1], q(T) * isincos_any(q(ph[0]))[0])
    gr, gi, ge = fixz(gam, W)
    E, eta = [None] * N, [None] * N
    eta0 = None
    for j in range(K):
        b = info["slopes"][j] / N                      # phase increment per cell (exact)
        s, c = isincos_any(q(ph[j]))
        z = Z(c, s)                                    # e^(i phi) at the left end of the current cell
        sb, cb = isincos_any(q(b))
        step = Z(cb, sb)
        coef = [T * b ** k / factorial(k) for k in range(Q + 1)]     # T b^k / k!
        tail = up(T * abs(b) ** (Q + 1) / factorial(Q + 1))
        for i in range(M[j], M[j + 1]):
            row, err = [], tail
            for k in range(Q + 1):
                # c_k = -T (i b)^k / k! e^(i phi_i) = -coef_k i^k z, plus gamma for k = 0
                m = k % 4
                zr, zi = (z.re, z.im) if m == 0 else (-z.im, z.re) if m == 1 else (-z.re, -z.im) if m == 2 else (
                    z.im, -z.re)
                cr, ci = -(zr * q(coef[k])), -(zi * q(coef[k]))
                if k == 0:
                    cr, ci = cr + gam.re, ci + gam.im
                if i == 0 and k == 0:
                    row.append((0, 0))                 # gamma - T e^(i psi_0) = 0 exactly
                    continue
                nr, er = fix(cr, W)
                ni, ei = fix(ci, W)
                row.append((nr, ni))
                err = HI.add(err, HI.add(er, ei))
            E[i], eta[i] = row, err
            if i == 0:
                eta0 = err                             # every term of cell 0 carries a factor theta^k, k >= 1
            z = zmul(z, step)
    for i in range(H, N):
        E[i] = [(gr, gi)] + [(0, 0)] * Q
        eta[i] = ge
    return E, eta, eta0, M


# ---- the trial function: floating-point collocation (only finds rho~; nothing below depends on it) ---------------
def witness(prof, N, P, W, Lmax=14):
    """rho~ on cell m is sum_k r_mk theta^k; returns R[m] = [(re, im)] with r = (re + i im) 2^-W exactly."""
    T, ph, kn = prof["T"], prof["phases"], prof["knots"]
    K = len(kn) - 1
    M = [int(u * N) for u in kn]
    b = [float((ph[j + 1] - ph[j]) / (kn[j + 1] - kn[j])) / N for j in range(K)]
    G = [float(T) * cmath.exp(1j * float(ph[j])) for j in range(K)]
    gam = G[0]
    Lj = [M[j + 1] - M[j] for j in range(K)]
    s = [(1 - math.cos((2 * c + 1) * math.pi / (2 * (P + 1)))) / 2 for c in range(P + 1)]
    pts = s + [1.0]
    fa = [math.factorial(i) for i in range(P + Lmax + 3)]
    Ztab = [[[fa[k] / fa[l + k + 1] * x ** (l + k + 1) for x in pts] for k in range(P + 1)] for l in range(Lmax + 1)]
    nterm = []
    for j in range(K):
        bb, l, t = abs(b[j]), 0, 1.0
        while l < Lmax and t > 1e-19:
            l += 1
            t = t * bb / (l + 1)
        nterm.append(l + 1)
    ib = [[(1j * b[j]) ** l for l in range(Lmax + 1)] for j in range(K)]
    eib = [[cmath.exp(1j * b[j] * x) for x in pts] for j in range(K)]
    eibL = [cmath.exp(1j * b[j] * Lj[j]) for j in range(K)]
    Bm = [[-gam * s[c] ** (k + 1) / (k + 1) + G[0] * sum(ib[0][l] * Ztab[l][k][c] for l in range(nterm[0]))
           for k in range(P + 1)] for c in range(P + 1)]
    V = [[s[c] ** k for k in range(P + 1)] for c in range(P + 1)]
    acc = [None] * N
    Eg, Ehat = 0j, [0j] * K
    R = []
    jn = 0
    scale = float(1 << W)
    for n in range(N):
        while jn < K - 1 and n >= M[jn + 1]:
            jn += 1
        a_n = acc[n] or [0j] * (P + 1)
        acc[n] = None
        rhs = []
        for c in range(P + 1):
            dl = gam - G[jn] * cmath.exp(1j * b[jn] * (n + s[c] - M[jn])) if n < M[K] else gam
            rhs.append(dl * N + gam * Eg - G[0] * eib[0][c] * Ehat[0] + a_n[c])
        A = [[(n + s[c]) * V[c][k] + Bm[c][k] for k in range(P + 1)] for c in range(P + 1)]
        r = _solve(A, rhs)
        # round to the exact dyadic trial; continue the march with the rounded values
        rr = [(int(round(z.real * scale)), int(round(z.imag * scale))) for z in r]
        R.append(rr)
        r = [complex(x, y) / scale for x, y in rr]
        Zv = [[sum(r[k] * Ztab[l][k][c] for k in range(P + 1)) for c in range(P + 2)] for l in range(Lmax + 1)]
        Eg += Zv[0][P + 1]
        for j in range(K):
            vals = [eib[j][c] * Ehat[j] + sum(ib[j][l] * Zv[l][c] for l in range(nterm[j])) for c in range(P + 2)]
            Ehat[j] = vals[P + 1]
            if j >= 1 and n + M[j] < N:
                tgt = acc[n + M[j]] = acc[n + M[j]] or [0j] * (P + 1)
                for c in range(P + 1):
                    tgt[c] -= G[j] * vals[c]
            if n + M[j + 1] < N:
                tgt = acc[n + M[j + 1]] = acc[n + M[j + 1]] or [0j] * (P + 1)
                f = G[j] * eibL[j]
                for c in range(P + 1):
                    tgt[c] += f * vals[c]
    return R


def _solve(A, y):
    n = len(y)
    A = [row[:] + [y[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(A[i][c]))
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] / piv
            if f:
                Ai, Ac = A[i], A[c]
                for k in range(c, n + 1):
                    Ai[k] -= f * Ac[k]
    x = [0j] * n
    for i in range(n - 1, -1, -1):
        x[i] = (A[i][n] - sum(A[i][k] * x[k] for k in range(i + 1, n))) / A[i][i]
    return x


# ---- exact two-dimensional convolution of integer arrays (Kronecker substitution in decimal integers) -------------
def _pack(A, S, D, nrows):
    zero = "0" * D
    pos, neg = [], []
    for i in range(nrows - 1, -1, -1):
        row = A[i]
        for a in range(S - 1, -1, -1):
            v = row[a] if a < len(row) else 0
            if v > 0:
                pos.append(str(v).zfill(D))
                neg.append(zero)
            elif v < 0:
                neg.append(str(-v).zfill(D))
                pos.append(zero)
            else:
                pos.append(zero)
                neg.append(zero)
    return EXACT.subtract(Decimal("".join(pos)), Decimal("".join(neg)))


def conv2d(X, Y, S, NZ):
    """Z[n][c] = sum over i + j = n, a + b = c of X[i][a] Y[j][b], for n < NZ and c < S, exactly; every row of X and
    Y has length <= S and the caller guarantees a + b < S for every nonzero pair. The arrays are packed into decimal
    integers with D digits per slot, multiplied in exact decimal arithmetic (any rounding would raise), and unpacked;
    D is chosen so that every slot of the product is below half of 10^D in modulus."""
    NX, NY = len(X), len(Y)
    mx = max((abs(v) for row in X for v in row), default=0)
    my = max((abs(v) for row in Y for v in row), default=0)
    bound = mx * my * min(NX, NY) * S
    D = max(len(str(2 * bound)), len(str(mx)), len(str(my))) + 1
    assert len(str(mx)) < D and len(str(my)) < D
    prod = EXACT.multiply(_pack(X, S, D, NX), _pack(Y, S, D, NY))
    npos = (NX + NY - 1) * S
    half = 5 * 10 ** (D - 1)
    s = str(EXACT.add(prod, Decimal(("5" + "0" * (D - 1)) * npos)))
    assert len(s) <= npos * D
    s = s.zfill(npos * D)
    L = len(s)
    out = []
    for n in range(NZ):
        base = L - n * S * D
        out.append([int(s[base - (c + 1) * D:base - c * D]) - half for c in range(S)])
    return out


def conv2d_complex(X, Y, S, NZ):
    """The same for Gaussian-integer arrays given as rows of (re, im) pairs (three real products, Gauss's trick)."""
    def part(A, k, n):
        return [[v[k] for v in row] + [0] * (n - len(row)) for row in A]
    Xr, Xi, Yr, Yi = part(X, 0, S), part(X, 1, S), part(Y, 0, S), part(Y, 1, S)
    P1 = conv2d(Xr, Yr, S, NZ)
    P2 = conv2d(Xi, Yi, S, NZ)
    P3 = conv2d([[u + v for u, v in zip(r1, r2)] for r1, r2 in zip(Xr, Xi)],
                [[u + v for u, v in zip(r1, r2)] for r1, r2 in zip(Yr, Yi)], S, NZ)
    return [[(p1 - p2, p3 - p1 - p2) for p1, p2, p3 in zip(a, b, c)] for a, b, c in zip(P1, P2, P3)]


def reversed_coeffs(row):
    """Coefficients of f(1 - s) from those of f(theta) (exact integers)."""
    n = len(row)
    out = []
    for k in range(n):
        sr = si = 0
        for a in range(k, n):
            cb = comb(a, k)
            sr += cb * row[a][0]
            si += cb * row[a][1]
        sg = -1 if k % 2 else 1
        out.append((sg * sr, sg * si))
    return out


# ---- Chebyshev bound for the sup of a polynomial on [0, 1] ---------------------------------------------------------
_CHEB = {}


def cheb01_matrix(D):
    """2^(2D) c_km with theta^k = sum_m c_km T_m(2 theta - 1): c_km = 2^(1-2k) C(2k, k-m) (m >= 1), 2^(-2k) C(2k, k)."""
    if D not in _CHEB:
        _CHEB[D] = [[(2 if m else 1) * comb(2 * k, k - m) * 4 ** (D - k) if m <= k else 0 for m in range(D + 1)]
                    for k in range(D + 1)]
    return _CHEB[D]


def sup01(qr, qi):
    """An integer upper bound for 2^(2D) sup_{[0,1]} |sum_k (qr_k + i qi_k) theta^k|, D = len - 1: the sum of the moduli
    of the Chebyshev coefficients, computed exactly."""
    D = len(qr) - 1
    Cm = cheb01_matrix(D)
    tot = 0
    for m in range(D + 1):
        cr = sum(qr[k] * Cm[k][m] for k in range(m, D + 1) if qr[k])
        ci = sum(qi[k] * Cm[k][m] for k in range(m, D + 1) if qi[k])
        tot += math.isqrt(cr * cr + ci * ci) + 1
    return tot


# ---- the residual on every cell -----------------------------------------------------------------------------------
def residual_bounds(N, P, Q, E, R, eta, eta0, WE, WR, out=print):
    """r[n] >= sup |R(t)/t| over cell n (decimals), from the exact residual polynomial of rho~ against delta~ plus the
    bounds for delta - delta~."""
    S = P + Q + 1                      # c = a + b <= P + Q
    DR = P + Q + 1                     # degree of the residual polynomial
    FD = factorial(DR)
    t0 = time.time()
    X = [[(factorial(a) * v[0], factorial(a) * v[1]) for a, v in enumerate(row)] for row in E]
    Y = [[(factorial(b) * v[0], factorial(b) * v[1]) for b, v in enumerate(row)] for row in R]
    U = conv2d_complex(X, Y, S, N)
    Xc = [[(factorial(a) * v[0], factorial(a) * v[1]) for a, v in enumerate(reversed_coeffs(row))] for row in E]
    Yc = [[(factorial(b) * v[0], factorial(b) * v[1]) for b, v in enumerate(reversed_coeffs(row))] for row in R]
    V = conv2d_complex(Xc, Yc, S, N - 1)
    t1 = time.time()
    # V-part: (1 - theta)^(c+1) / (c+1)! expanded, scaled by FD
    VM = [[(FD // factorial(c + 1)) * comb(c + 1, k) * (-1) ** k for k in range(DR + 1)] for c in range(S)]
    UM = [FD // factorial(c + 1) for c in range(S)]
    sE = N * (1 << WR) * FD
    sR = (1 << WE) * FD
    scale = Decimal(1 << (WE + WR)) * Decimal(FD)          # Q~ times scale has integer coefficients
    rnorm = max(dsum(dint(v[0], WR) + dint(v[1], WR) for v in row) for row in R)   # sup |rho~| <= sum |r_k|
    etamax = max(eta)
    glob = HI.multiply(etamax, rnorm)
    r = []
    worst = (ZERO, -1)
    for n in range(N):
        qr, qi = [0] * (DR + 1), [0] * (DR + 1)
        for k in range(Q + 1):
            qr[k] += sE * E[n][k][0]
            qi[k] += sE * E[n][k][1]
        Un = U[n]
        for c in range(S):
            qr[c + 1] += UM[c] * Un[c][0]
            qi[c + 1] += UM[c] * Un[c][1]
        if n >= 1:
            Vn = V[n - 1]
            for c in range(S):
                vr, vi = Vn[c]
                if vr or vi:
                    row = VM[c]
                    for k in range(c + 2):
                        qr[k] += row[k] * vr
                        qi[k] += row[k] * vi
        Rn = R[n]
        for k in range(P + 1):
            qr[k] -= sR * n * Rn[k][0]
            qi[k] -= sR * n * Rn[k][1]
            qr[k + 1] -= sR * Rn[k][0]
            qi[k + 1] -= sR * Rn[k][1]
        if n == 0:
            assert qr[0] == 0 and qi[0] == 0
            sup = sup01(qr[1:], qi[1:])
            val = HI.divide(Decimal(sup), HI.multiply(scale, Decimal(4 ** (DR - 1))))
            rn = HI.add(HI.add(val, HI.multiply(Decimal(N), eta0)), glob)
        else:
            sup = sup01(qr, qi)
            val = HI.divide(Decimal(sup), HI.multiply(scale, Decimal(4 ** DR)))
            rn = HI.add(HI.divide(HI.add(val, HI.multiply(Decimal(N), eta[n])), Decimal(n)), glob)
        r.append(rn)
        if rn > worst[0]:
            worst = (rn, n)
    out("  residual: two exact 2-D convolutions in %.1f s, cell bounds in %.1f s; max |R/t| <= %s at cell %d, "
        "delta-error term %s" % (t1 - t0, time.time() - t1, fmt(worst[0], 4), worst[1], fmt(glob, 3)))
    return r


# ---- from the residual to e = rho - rho~: a block majorant ----------------------------------------------------------
def majorant(prof, info, N, r, g, out=print):
    """E_k >= |e| on block k = [k w, (k+1) w), w = g h, from the comparison principle: |e| <= r + K|e| with
    (K f)(t) = (1/t) int_0^t |d(u)| f(t - u) du (norm <= kappa < 1), so any bounded E >= 0 with E >= r + K E bounds |e|.
    On block k >= 1: (K E)(t) <= (1/k) sum_{i=0..k} D_i max(E_{k-i}, E_{k-i-1}), D_i >= sup |d| on block i; on block 0,
    |d(u)| <= M0 u gives (K E)(t) <= E_0 M0 w / 2. Returns (E, EN, ED): EN >= int_0^1 (1-v)^-a |e|,
    ED >= 2/(1-a) int_0^(1/2) (1/2 - v)^(1-a) |e| >= int_0^(1/2) |p - p~| dt / (1 - t)."""
    T, ph, kn = prof["T"], prof["phases"], prof["knots"]
    H = N // 2
    assert H % g == 0
    KB, w = N // g, Fraction(g, N)
    rB = [max(r[k * g:(k + 1) * g]) for k in range(KB)]
    slopes = info["slopes"]

    def phi(u):
        for j in range(len(kn) - 1):
            if u <= kn[j + 1]:
                return ph[j] + slopes[j] * (u - kn[j])
        return ph[-1]
    D = []
    for k in range(KB):
        if k < KB // 2:
            lo, hi = k * w, (k + 1) * w
            pts = [lo, hi] + [u for u in kn if lo < u < hi]
            D.append(up(T * max(abs(phi(u) - ph[0]) for u in pts)))
        else:
            D.append(up(T))
    lam0 = max(abs(slopes[j]) for j in range(len(kn) - 1) if kn[j] < w)
    m0 = up(T * lam0 * w / 2)
    assert m0 < 1 and D[0] < 1
    E = [HI.divide(rB[0], LO.subtract(ONE, m0))]
    for k in range(1, KB):
        s = HI.multiply(D[0], E[k - 1])
        for i in range(1, k + 1):
            l = k - i
            e = E[l] if l == 0 or E[l] >= E[l - 1] else E[l - 1]
            s = HI.add(s, HI.multiply(D[i], e))
        den = LO.subtract(ONE, HI.divide(D[0], Decimal(k)))
        E.append(HI.divide(HI.add(rB[k], HI.divide(s, Decimal(k))), den))
    A = info["a"].hi                              # (1-v)^-a <= (1-v)^-A and (1/2-v)^(1-a) <= (1/2-v)^(1-A)

    def pw(x, e):                                 # x^e for rationals x > 0, as an interval
        return iexp(q(e) * iln(q(x)))
    EN = ED = ZERO
    one_A, two_A = LO.subtract(ONE, A), LO.subtract(Decimal(2), A)
    for k in range(KB):
        x0, x1 = 1 - k * w, 1 - (k + 1) * w
        hi0 = pw(x0, 0).hi if x0 == 1 else iexp(I(one_A) * iln(q(x0))).hi
        lo1 = iexp(I(one_A) * iln(q(x1))).lo if x1 > 0 else ZERO
        EN = HI.add(EN, HI.multiply(E[k], HI.divide(HI.subtract(hi0, lo1), one_A)))
        if k < KB // 2:
            y0, y1 = Fraction(1, 2) - k * w, Fraction(1, 2) - (k + 1) * w
            hi0 = iexp(I(two_A) * iln(q(y0))).hi
            lo1 = iexp(I(two_A) * iln(q(y1))).lo if y1 > 0 else ZERO
            ED = HI.add(ED, HI.multiply(E[k], HI.divide(HI.subtract(hi0, lo1), two_A)))
    ED = HI.multiply(ED, HI.divide(Decimal(2), one_A))
    out("  majorant: %d blocks of %d cells; max E = %s (block %d); EN <= %s, ED <= %s"
        % (KB, g, fmt(max(E), 3), E.index(max(E)), fmt(EN, 3), fmt(ED, 3)))
    return E, EN, ED


# ---- constants in gamma ------------------------------------------------------------------------------------------------
class Consts:
    def __init__(self, info, Kmax, P, DA, SF):
        a, b = info["a"], info["b"]
        self.a, self.b = a, b
        self.gam = Z(a, b)
        # k^-gamma = k^-a (cos(b ln k) - i sin(b ln k)), from the primes, multiplicatively
        spf = list(range(Kmax + 1))
        for p in range(2, int(Kmax ** 0.5) + 1):
            if spf[p] == p:
                for m in range(p * p, Kmax + 1, p):
                    if spf[m] == m:
                        spf[m] = p
        tbl = [None, Z(I(1), I(0))]
        for k in range(2, Kmax + 1):
            p = spf[k]
            if p == k:
                lk = iln(q(k))
                mod = iexp(-(a * lk))
                s, c = isincos_any(b * lk)
                tbl.append(Z(mod * c, -(mod * s)))
            else:
                tbl.append(zmul(tbl[p], tbl[k // p]))
        self.tbl = tbl
        l2 = iln(q(2))
        mod = iexp(a * l2)
        s, c = isincos_any(b * l2)
        self.two_g = Z(mod * c, mod * s)                       # 2^gamma
        one = Z(I(1), I(0))
        # B_k = B(k+1, 1-gamma) = k! / prod_{i=1..k+1} (i - gamma)
        self.B = []
        cur = zinv_z(Z(I(1) - a, -b))
        for k in range(P + 1):
            if k:
                cur = zmul(zscale(cur, k), zinv_z(Z(I(k + 1) - a, -b)))
            self.B.append(cur)
        # binom(-gamma, i)
        self.bin = [one]
        for i in range(1, DA + 2):
            self.bin.append(zscale(zmul(self.bin[-1], Z(-a - I(i - 1), -b)), Fraction(1, i)))
        # falling factorial (-gamma)(-gamma-1)...(-gamma-s+1)
        self.fall = [one]
        for s_ in range(1, SF + 1):
            self.fall.append(zmul(self.fall[-1], Z(-a - I(s_ - 1), -b)))
        # (gamma)_j / j!
        self.poch = [one]
        for j in range(1, 400):
            self.poch.append(zscale(zmul(self.poch[-1], Z(a + I(j - 1), b)), Fraction(1, j)))


# ---- the numerator p~(1) ---------------------------------------------------------------------------------------------
def numerator(N, P, R, WR, cs, out=print):
    """p~(1) = 1 + h^(1-gamma) sum_{c=1..N} int_0^1 (c - y)^-gamma rho~_{N-c}(y) dy (tau units). For c >= 2,
    (c - y)^-gamma = c^-gamma sum_j (gamma)_j/j! (y/c)^j with |(gamma)_j/j!| <= 1 (tail bounded); for c = 1 the
    moments are B(k+1, 1-gamma)."""
    a, b = cs.a, cs.b
    t0 = time.time()
    lcm = [0] * 400
    for j in range(400):
        L = 1
        for k in range(j + 1, j + P + 2):
            L = L * k // math.gcd(L, k)
        lcm[j] = L
    two_W = 1 << WR
    tot = Z(I(0), I(0))
    tail = ZERO
    a_lo = a.lo
    for c in range(2, N + 1):
        row = R[N - c]
        l1 = dsum(HI.divide(HI.add(dint(v[0], WR), dint(v[1], WR)), Decimal(k + 1)) for k, v in enumerate(row))
        # terms: c^-(J+1)/(1 - 1/c) l1 < 1e-24 l1
        J = int(math.ceil(24 * math.log(10) / math.log(c)))
        S = Z(I(0), I(0))
        for j in range(J + 1):
            L = lcm[j]
            nr = sum(v[0] * (L // (j + k + 1)) for k, v in enumerate(row))
            ni = sum(v[1] * (L // (j + k + 1)) for k, v in enumerate(row))
            den = L * two_W * c ** j
            xr, xi = q(Fraction(nr, den)), q(Fraction(ni, den))
            pj = cs.poch[j]
            S = S + Z(pj.re * xr - pj.im * xi, pj.re * xi + pj.im * xr)
        tot = tot + zmul(cs.tbl[c], S)
        t = HI.multiply(HI.divide(Decimal(1), LO.subtract(ONE, HI.divide(ONE, Decimal(c)))),
                        HI.multiply(l1, iexp(-(I(a_lo) * iln(q(c)))).hi))
        tail = HI.add(tail, HI.multiply(t, HI.divide(ONE, LO.power(Decimal(c), J + 1))))
    row = R[N - 1]
    for k, v in enumerate(row):
        x = Z(q(Fraction(v[0], two_W)), q(Fraction(v[1], two_W)))
        tot = tot + zmul(x, cs.B[k])
    tot = tot.widen(tail)
    lN = iln(q(N))
    mod = iexp(a * lN)
    s, c_ = isincos_any(b * lN)
    hg = Z(mod * c_ * q(Fraction(1, N)), mod * s * q(Fraction(1, N)))     # h^(1-gamma) = N^gamma / N
    p1 = zmul(hg, tot)
    p1 = Z(p1.re + I(1), p1.im)
    hi = zabs_hi(p1)
    re_lo = p1.re.lo if p1.re.lo > 0 else (p1.re.hi.copy_abs() if p1.re.hi < 0 else ZERO)
    im_lo = p1.im.lo if p1.im.lo > 0 else (p1.im.hi.copy_abs() if p1.im.hi < 0 else ZERO)
    lo = NEAR.next_minus(NEAR.sqrt(LO.add(LO.multiply(re_lo, re_lo), LO.multiply(im_lo, im_lo))))
    out("  p~(1) = %s + %s i (series tails %s), |p~(1)| in [%s, %s]   (%.1f s)"
        % (fmt(p1.re.lo, 16), fmt(p1.im.lo, 16), fmt(tail, 3), fmt(lo, 16), fmt(hi, 16), time.time() - t0))
    return lo, hi


# ---- the denominator ------------------------------------------------------------------------------------------------
def _horner(c, x):
    v = 0j
    for z in reversed(c):
        v = v * x + z
    return v


def _rsolve(A, y):
    """Real linear solve (Gaussian elimination, partial pivoting), floats."""
    n = len(y)
    A = [row[:] + [y[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(A[i][c]))
        A[c], A[p] = A[p], A[c]
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            for k in range(c, n + 1):
                A[i][k] -= f * A[c][k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (A[i][n] - sum(A[i][k] * x[k] for k in range(i + 1, n))) / A[i][i]
    return x


def _ufit(f, pts, du, WU):
    """Integer coefficients (scale 2^WU) of a polynomial u that interpolates the floats 1/|f(x)| at the points."""
    vals = [1.0 / abs(f(x)) for x in pts]
    u = _rsolve([[x ** k for k in range(du + 1)] for x in pts], vals)
    return [int(round(v * (1 << WU))) for v in u]


def _pmul(p1, p2, deg=None):
    """Product of two polynomials with Gaussian-integer coefficients (lists of (re, im))."""
    n1, n2 = len(p1), len(p2)
    out_r, out_i = [0] * (n1 + n2 - 1), [0] * (n1 + n2 - 1)
    for i, (ar, ai) in enumerate(p1):
        if not (ar or ai):
            continue
        for j, (br, bi) in enumerate(p2):
            out_r[i + j] += ar * br - ai * bi
            out_i[i + j] += ar * bi + ai * br
    return list(zip(out_r, out_i))


def _rmul(p1, p2):
    """Product of two real integer polynomials."""
    out = [0] * (len(p1) + len(p2) - 1)
    for i, a in enumerate(p1):
        if a:
            for j, b in enumerate(p2):
                out[i + j] += a * b
    return out


_CHEBC = {}


def supc(p):
    """An integer upper bound for 4^D sup_{|xi| <= 1/2} |sum_k p_k xi^k| (real integer coefficients), D = len - 1: the sum
    of the moduli of the Chebyshev coefficients on [-1/2, 1/2] (xi^k = 2^-k t^k, t^k = sum_m 2^(1-k) C(k, (k-m)/2) T_m)."""
    D = len(p) - 1
    if D not in _CHEBC:
        _CHEBC[D] = [[(2 if m else 1) * comb(k, (k - m) // 2) * 2 ** (2 * D - 2 * k) if m <= k and (k - m) % 2 == 0
                      else 0 for m in range(D + 1)] for k in range(D + 1)]
    Cm = _CHEBC[D]
    return sum(abs(sum(p[k] * Cm[k][m] for k in range(m, D + 1, 2))) for m in range(D + 1)) + D + 1


def reg_cell(A, WA, n, N, du=6, WU=40):
    """Bounds for int_{-1/2}^{1/2} |A(xi)| w(xi) dxi, w = 1/(1 - h(n + 1/2 + xi)), A = sum (re + i im) 2^-WA xi^k:
    |A| >= g = (3/2) X u - (1/2) X^2 u^3 with X = |A|^2 (an exact real polynomial) and u >= 0; the integral of g w is
    computed exactly (w by its geometric series in xi, tail bounded). Upper bound: |A| <= g + eta^2 (1+eta)(3+eta) /
    (2 u_min) with eta >= sup |u^2 X - 1| < 1. Returns (lower, upper) as decimals."""
    D = len(A) - 1
    fl = [complex(x, y) / float(1 << WA) for x, y in A]
    pts = [math.cos((2 * c + 1) * math.pi / (2 * (du + 1))) / 2 for c in range(du + 1)]
    U = _ufit(lambda x: _horner(fl, x), pts, du, WU)
    umin_num = (U[0] << du) - sum(abs(U[k]) << (du - k) for k in range(1, du + 1))
    assert umin_num > 0, "u >= 0 not certified"
    Xp = [0] * (2 * D + 1)
    for i, (ar, ai) in enumerate(A):
        Xp[2 * i] += ar * ar + ai * ai
        for j in range(i + 1, D + 1):
            br, bi = A[j]
            Xp[i + j] += 2 * (ar * br + ai * bi)
    X2 = [0] * (4 * D + 1)
    for i, x in enumerate(Xp):
        if x:
            X2[2 * i] += x * x
            for j in range(i + 1, 2 * D + 1):
                X2[i + j] += 2 * x * Xp[j]
    Y1 = _rmul(Xp, U)
    U3 = _rmul(_rmul(U, U), U)
    Y2 = _rmul(X2, U3)
    sc1 = 3 << (2 * WA + 2 * WU)
    G2 = [-y for y in Y2]
    for k, y in enumerate(Y1):
        G2[k] += sc1 * y
    S = 1 << (4 * WA + 3 * WU + 1)               # g = G2 / S
    nG = len(G2)
    Lw = _nterms(Fraction(1, 2 * N - 2 * n - 1))      # lambda/2 = 1/(2N - 2n - 1)
    # moments of xi^j on [-1/2, 1/2]: [j even] / ((j+1) 2^j)
    Lmu = 1
    for j in range(1, nG + Lw + 1):
        Lmu = Lmu * j // math.gcd(Lmu, j)
    Lmu <<= nG + Lw
    Il = []
    for l in range(Lw + 1):
        s = 0
        for k, gk in enumerate(G2):
            j = k + l
            if gk and j % 2 == 0:
                s += gk * (Lmu // ((j + 1) << j))
        Il.append(s)
    cn = Fraction(2 * N - 2 * n - 1, 2 * N)
    lam = Fraction(2, 2 * N - 2 * n - 1)
    val = sum(Fraction(Il[l]) * lam ** l for l in range(Lw + 1)) / (cn * Lmu * S)
    supg = HI.divide(Decimal(sum(abs(gk) >> k for k, gk in enumerate(G2)) + nG), Decimal(S))
    tailw = HI.multiply(supg, up(1 / cn * (lam / 2) ** (Lw + 1) / (1 - lam / 2)))
    lower = LO.subtract(frac_lo(val), tailw)
    # upper bound
    UY = _rmul(U, Y1)
    UY[0] -= 1 << (2 * WA + 2 * WU)
    eta = HI.divide(Decimal(sum(abs(c) >> k for k, c in enumerate(UY)) + len(UY)), Decimal(1 << (2 * WA + 2 * WU)))
    if not eta < ONE:                              # Chebyshev coefficients on [-1/2, 1/2] instead
        eta = HI.divide(Decimal(supc(UY)), Decimal((1 << (2 * WA + 2 * WU)) * 4 ** (len(UY) - 1)))
    assert eta < ONE, "eta >= 1"
    umin = LO.divide(Decimal(umin_num), Decimal(1 << (du + WU)))
    corr = HI.divide(HI.multiply(HI.multiply(HI.multiply(eta, eta), HI.add(ONE, eta)), HI.add(Decimal(3), eta)),
                     LO.multiply(Decimal(2), umin))
    supw = up(Fraction(N, N - n - 1))
    upper = HI.add(HI.add(frac_hi(val), tailw), HI.multiply(corr, supw))
    umax = HI.divide(Decimal((U[0] << du) + sum(abs(U[k]) << (du - k) for k in range(1, du + 1))),
                     Decimal(1 << (du + WU)))
    return lower, upper, supw, U, eta, umax


def first_order(A, WA, U, WU, J, WJ, n, N, inv1, W1, e1max):
    """An interval for Re int_0^1 S(x) conj(A(x - 1/2)) u(x - 1/2) w(x) dx, S = h x^(1 - gamma) J(x),
    w = 1/(1 - h(n + x)): A, u in xi = x - 1/2 (integers at 2^-WA, 2^-WU), J in x (Gaussian integers at 2^-WJ),
    inv1[j] = 2^W1 times a rounded 1/(j + 2 - gamma) with error <= e1max. The polynomial parts are exact; the moments
    are int_0^1 x^(k + 1 - gamma) x^l dx = 1/(k + l + 2 - gamma), w = (N/(N-n)) sum_l (x/(N-n))^l (tail bounded)."""
    D, du = len(A) - 1, len(U) - 1
    Dt = D + du
    Br, Bi = [0] * (Dt + 1), [0] * (Dt + 1)
    for i, (ar, ai) in enumerate(A):
        for j, u in enumerate(U):
            Br[i + j] += ar * u
            Bi[i + j] -= ai * u
    Bx = []
    for k in range(Dt + 1):                           # p(x - 1/2), times 2^Dt
        sr = si = 0
        for i in range(k, Dt + 1):
            cb = comb(i, k) * (-1) ** (i - k) << (Dt - (i - k))
            sr += Br[i] * cb
            si += Bi[i] * cb
        Bx.append((sr, si))
    c = _pmul(J, Bx)
    Nn = N - n
    Lw = _nterms(Fraction(1, Nn))
    assert len(c) + Lw <= len(inv1)
    tot = 0
    absc = 0
    for k, (cr, ci) in enumerate(c):
        if not (cr or ci):
            continue
        mr = mi = 0
        for l in range(Lw + 1):
            f = Nn ** (Lw - l)
            vr, vi = inv1[k + l]
            mr += f * vr
            mi += f * vi
        tot += cr * mr - ci * mi
        absc += abs(cr) + abs(ci)
    scale = Decimal(1 << (WJ + WA + WU + Dt))
    den = HI.multiply(HI.multiply(scale, Decimal(1 << W1)), Decimal(Nn ** Lw))
    lo, hi = LO.divide(Decimal(tot), den), HI.divide(Decimal(tot), den)
    # errors: the rounded moments, and the tail of the series for w
    ac = HI.divide(Decimal(absc), scale)
    err = HI.multiply(ac, HI.multiply(e1max, up(Fraction(Nn, Nn - 1))))
    err = HI.add(err, HI.multiply(ac, up(Fraction(1, Nn ** (Lw + 1)) * Fraction(Nn, Nn - 1))))
    f = q(Fraction(1, Nn))                           # h N / (N - n)
    return LO.multiply(LO.subtract(lo, err), f.lo if lo >= err else f.hi), HI.multiply(HI.add(hi, err),
                                                                                   f.hi if hi + err >= 0 else f.lo)


def second_order(J, WJ, U, WU, n, N, cs):
    """An upper bound for int_0^1 |S(x)|^2 u(x - 1/2) w(x) dx, S = h x^(1 - gamma) J(x): |S|^2 = h^2 x^(2 - 2a) |J|^2,
    |J(x)|^2 = J(x) conj(J)(x) for real x, u in xi = x - 1/2 (moments 1/(k + l + 3 - 2a), interval arithmetic)."""
    a = cs.a
    du = len(U) - 1
    Ux = []
    for k in range(du + 1):                          # u(x - 1/2), times 2^du
        Ux.append(sum(U[i] * comb(i, k) * (-1) ** (i - k) << (du - (i - k)) for i in range(k, du + 1)))
    JJ = _pmul(J, [(x, -y) for x, y in J])
    c = _rmul([x for x, _ in JJ], Ux)               # imaginary parts of |J|^2 cancel exactly
    assert all(y == 0 for _, y in JJ)
    Nn = N - n
    Lw = _nterms(Fraction(1, Nn))
    tot = I(0)
    absc = 0
    for k, ck in enumerate(c):
        if ck:
            mom = I(0)
            for l in range(Lw + 1):
                mom = mom + q(Fraction(1, Nn ** l)) / (I(k + l + 3) - a * I(2))
            tot = tot + I(Decimal(ck)) * mom
            absc += abs(ck)
    sc = Decimal(1 << (2 * WJ + WU + du))
    tail = HI.divide(HI.multiply(Decimal(absc), up(Fraction(1, Nn ** (Lw + 1)) * Fraction(Nn, Nn - 1))), sc)
    val = HI.add(HI.divide(tot.hi, sc), tail)
    return HI.multiply(val, up(Fraction(N, Nn) * Fraction(1, N * N)))


def gen_cell(F, Sf, omega, n, N, cs, du=12, WU=40):
    """Bounds for int_0^1 x^omega |f(x)| w(x) dx with f = sum over types t = (m1, m2) of x^(m1 e1 + m2 e2) F_t(x) / Sf,
    e1 = 1 - gamma, e2 = 1 - conj(gamma), F_t polynomials with Gaussian-integer coefficients, w = 1/(1 - h(n + x)),
    omega = 0 or -a. Same inequality as reg_cell, now with generalized monomials x^(k + m1 e1 + m2 e2), integrated by
    their exact moments (interval arithmetic); |x^(k + E)| <= 1 on [0, 1] bounds every tail."""
    a, b = cs.a, cs.b
    gam_f = complex(float(NEAR.divide(NEAR.add(a.lo, a.hi), 2)), float(NEAR.divide(NEAR.add(b.lo, b.hi), 2)))
    e1f = 1 - gam_f
    sf = float(Sf)
    Ff = {t: [complex(x, y) / sf for x, y in p] for t, p in F.items()}

    def f_at(x):
        v = 0j
        for (m1, m2), p in Ff.items():
            v += cmath.exp((m1 * e1f + m2 * e1f.conjugate()) * math.log(x)) * _horner(p, x)
        return v
    pts = [(1 - math.cos((2 * c + 1) * math.pi / (2 * (du + 1)))) / 2 for c in range(du + 1)]
    U = _ufit(f_at, pts, du, WU)
    Cm = cheb01_matrix(du)
    ch = [sum(U[k] * Cm[k][m] for k in range(m, du + 1)) for m in range(du + 1)]
    umin_num = ch[0] - sum(abs(c) for c in ch[1:])
    assert umin_num > 0, "u >= 0 not certified"
    Fc = {(t[1], t[0]): [(x, -y) for x, y in p] for t, p in F.items()}

    def gmul(P1, P2):
        out = {}
        for t1, p1 in P1.items():
            for t2, p2 in P2.items():
                t = (t1[0] + t2[0], t1[1] + t2[1])
                pr = _pmul(p1, p2)
                o = out.setdefault(t, [])
                if len(o) < len(pr):
                    o.extend([(0, 0)] * (len(pr) - len(o)))
                for i, (x, y) in enumerate(pr):
                    o[i] = (o[i][0] + x, o[i][1] + y)
        return out
    Uc = {(0, 0): [(u, 0) for u in U]}
    X = gmul(F, Fc)
    Y1 = gmul(X, Uc)
    Y2 = gmul(gmul(X, X), gmul(gmul(Uc, Uc), Uc))
    sc1 = 3 * Sf * Sf << (2 * WU)
    G2 = {t: [(-x, -y) for x, y in p] for t, p in Y2.items()}
    for t, p in Y1.items():
        g = G2.setdefault(t, [])
        if len(g) < len(p):
            g.extend([(0, 0)] * (len(p) - len(g)))
        for k, (x, y) in enumerate(p):
            g[k] = (g[k][0] + sc1 * x, g[k][1] + sc1 * y)
    S = 2 * Sf ** 4 << (3 * WU)                   # g = G2 / S
    w0 = q(Fraction(N, N - n))
    mu = Fraction(1, N - n)
    Lw = _nterms(mu)
    tot = Z(I(0), I(0))
    supg = ZERO
    for t, p in G2.items():                       # sup |g| <= sum over types of sup |G2_t| / S (|x^E| <= 1)
        supg = HI.add(supg, HI.divide(Decimal(sup01([x for x, _ in p], [y for _, y in p])),
                                      HI.multiply(Decimal(S), Decimal(4 ** (len(p) - 1)))))
    for (m1, m2), p in G2.items():
        Er = I(m1 + m2) * (I(1) - a) + omega
        Ei = -(I(m1 - m2) * b)
        invs = [zinv_z(Z(Er + I(j + 1), Ei)) for j in range(len(p) + Lw + 1)]
        for k, (x, y) in enumerate(p):
            if not (x or y):
                continue
            mom = Z(I(0), I(0))
            ml = Fraction(1)
            for l in range(Lw + 1):
                mom = mom + zscale(invs[k + l], ml)
                ml *= mu
            gk = Z(I(Decimal(x)), I(Decimal(y)))
            tot = tot + zmul(gk, mom)
    tot = Z(tot.re * w0, tot.im * w0)
    Sd = Decimal(S)
    val_lo, val_hi = LO.divide(tot.re.lo, Sd), HI.divide(tot.re.hi, Sd)
    om_lo = omega.lo if isinstance(omega, I) else Decimal(omega)
    int_xw = HI.divide(ONE, LO.add(ONE, om_lo))          # int_0^1 x^omega dx
    tailw = HI.multiply(HI.multiply(supg, int_xw), up(Fraction(N, N - n) * mu ** (Lw + 1) / (1 - mu)))
    lower = LO.subtract(val_lo, tailw)
    UY = gmul(Uc, Y1)
    one = Sf * Sf << (2 * WU)
    p00 = UY.setdefault((0, 0), [(0, 0)])
    p00[0] = (p00[0][0] - one, p00[0][1])
    eta = ZERO
    for t, p in UY.items():                       # sup |u^2 X - 1|
        eta = HI.add(eta, HI.divide(Decimal(sup01([x for x, _ in p], [y for _, y in p])),
                                    HI.multiply(Decimal(one), Decimal(4 ** (len(p) - 1)))))
    assert eta < ONE, "eta >= 1"

    umin = LO.divide(Decimal(umin_num), Decimal((1 << (2 * du)) << WU))
    corr = HI.divide(HI.multiply(HI.multiply(HI.multiply(eta, eta), HI.add(ONE, eta)), HI.add(Decimal(3), eta)),
                     LO.multiply(Decimal(2), umin))
    supw = up(Fraction(N, N - n - 1))
    upper = HI.add(HI.add(val_hi, tailw), HI.multiply(HI.multiply(corr, supw), int_xw))
    return lower, upper, HI.multiply(supw, int_xw)


def _nterms(r):
    """The least L with r^(L+1) / (1 - r) < 10^-26, for a rational 0 < r < 1 (length of the series for w)."""
    L = 1
    while r ** (L + 1) / (1 - r) >= Fraction(1, 10 ** 26):
        L += 1
    return L


def _rdiv(num, den):
    """round(num / den) for integers, den > 0."""
    return (2 * num + den) // (2 * den)


def denominator(prof, info, N, P, R, WR, cs, M, DA=24, SF=40, out=print):
    """Lower (and, for the record, upper) bounds for int_0^(1/2) |p~(t)| dt / (1 - t)."""
    t0 = time.time()
    H = N // 2
    a, b = cs.a, cs.b
    W = WA = WG = WB = 64
    WK = WM = 64
    two = lambda e: Decimal(1 << e)                                     # noqa: E731
    # ---- far field (source cells m = n - c with c >= 2), one exact 2-D convolution ----
    Fs = [fixz(cs.fall[s], W) for s in range(SF + 1)]
    Tc = [None, None] + [fixz(cs.tbl[c], W) for c in range(2, H)]
    eT = max(t[2] for t in Tc[2:]) if H > 2 else ZERO
    Kar = [[(0, 0)] * (SF + 1) for _ in range(min(2, H))]
    eK = ZERO
    for s in range(SF + 1):
        e = HI.divide(HI.add(HI.add(HI.multiply(Decimal(factorial(s)), eT), Fs[s][2]), HI.multiply(eT, Fs[s][2])),
                      Decimal(16 ** s))
        eK = max(eK, e)
    eK = HI.add(eK, HI.divide(ONE, two(WK)))
    for c in range(2, H):
        tr, ti, _ = Tc[c]
        row = []
        for s in range(SF + 1):
            fr, fi, _ = Fs[s]
            den = (8 * c) ** s << (2 * W - WK)
            row.append((_rdiv(fr * tr - fi * ti, den), _rdiv(fr * ti + fi * tr, den)))
        Kar.append(row)
    # moments: M[m][SF - j] = round(2^WM (-1)^j 8^j / j! int_0^1 (y - 1/2)^j rho~_m(y) dy)
    gjk = []
    for j in range(SF + 1):
        row = []
        for k in range(P + 1):
            v = Fraction(0)
            for i in range(k + 1):
                if (j + i) % 2 == 0:
                    v += comb(k, i) * Fraction(1, 2 ** (k - i)) * Fraction(1, (j + i + 1) * 2 ** (j + i))
            row.append(v * (-1) ** j * Fraction(8 ** j, factorial(j)) * Fraction(1 << WM, 1 << WR))
        dj = 1
        for v in row:
            dj = dj * v.denominator // math.gcd(dj, v.denominator)
        gjk.append(([int(v * dj) for v in row], dj))
    Mar = []
    for m in range(H):
        Rm = R[m]
        row = [None] * (SF + 1)
        for j in range(SF + 1):
            nums, dj = gjk[j]
            row[SF - j] = (_rdiv(sum(nm * v[0] for nm, v in zip(nums, Rm)), dj),
                           _rdiv(sum(nm * v[1] for nm, v in zip(nums, Rm)), dj))
        Mar.append(row)
    eM = HI.divide(ONE, two(WM))
    SKsum = HI.divide(Decimal(sum(abs(x) + abs(y) for row in Kar for x, y in row)), two(WK))
    SMsum = HI.divide(Decimal(sum(abs(x) + abs(y) for row in Mar for x, y in row)), two(WM))
    npairs = Decimal(H * (SF + 1))
    e4 = Decimal("54.6")                       # sum_q 8^q / q! 2^-q <= e^4 < 54.6
    far_round = HI.multiply(e4, HI.add(HI.add(HI.multiply(eK, SMsum), HI.multiply(eM, SKsum)),
                                       HI.multiply(HI.multiply(eK, eM), npairs)))
    t1 = time.time()
    Prod = conv2d_complex(Kar, Mar, 2 * SF + 1, H) if H > 2 else None
    t2 = time.time()
    # ---- c = 1 extended: G_k(xi) = B_k (3/2 + xi)^(k+1-gamma), coefficients of xi^i for i <= DA ----
    th = zmul(cs.tbl[3], cs.two_g)                            # (3/2)^-gamma
    Gk, eGk = [], []
    for k in range(P + 1):
        s_ = Z(I(k + 1) - a, -b)
        base = zscale(zmul(cs.B[k], th), Fraction(3, 2) ** (k + 1))
        bi = Z(I(1), I(0))
        row, err = [], ZERO
        for i in range(DA + 1):
            if i:
                bi = zscale(zmul(bi, Z(s_.re - I(i - 1), s_.im)), Fraction(2, 3 * i))
            nr, ni, e = fixz(zmul(base, bi), WG)
            row.append((nr, ni))
            err = HI.add(err, HI.multiply(e, up(Fraction(1, 2 ** i))))
        # tail: the terms beta_i = binom(s, i) (2/3)^i 2^-i with i > DA; beta_(i+1)/beta_i = |s - i| / (i + 1) / 3
        # <= 1/3 there (|s - i| <= i - Re s + |Im s| <= i + 1 since i >= Re s = k + 1 - a and |Im s| < 1), so the
        # tail is at most 3/2 |binom(s, DA+1) (2/3)^(DA+1)| 2^-(DA+1), times |B_k (3/2)^(k+1-gamma)|
        assert DA + 1 >= k + 1
        bi = zscale(zmul(bi, Z(s_.re - I(DA), s_.im)), Fraction(2, 3 * (DA + 1)))
        tail = HI.multiply(HI.multiply(zabs_hi(bi), zabs_hi(base)), up(Fraction(3, 2) / 2 ** (DA + 1)))
        Gk.append(row)
        eGk.append(HI.add(err, tail))
    # ---- t^-gamma part: (n + 1/2)^-gamma binom(-gamma, q) (n + 1/2)^-q ----
    Bq = [fixz(zscale(cs.bin[i], Fraction(2 ** i)), W) for i in range(DA + 1)]
    Bint = [fixz(cs.B[k], WB) for k in range(P + 1)]
    Babs = [zabs_hi(cs.B[k]) for k in range(P + 1)]
    W1 = 80
    inv1, e1max = [], ZERO
    for j in range(P + DA + 6 + 2 + _nterms(Fraction(2, N)) + 2):
        vr, vi, e = fixz(zinv_z(Z(I(j + 2) - a, -b)), W1)
        inv1.append((vr, vi))
        e1max = max(e1max, e)
    l1 = [dsum(HI.divide(HI.add(dint(v[0], WR), dint(v[1], WR)), Decimal(k + 1)) for k, v in enumerate(R[m]))
          for m in range(H)]
    l1max = max(l1)
    t3 = time.time()
    knots = set(M[1:-1])
    hfac_lo = iexp((I(1) - a) * iln(q(Fraction(1, N)))).lo
    hfac_hi = iexp((I(1) - a) * iln(q(Fraction(1, N)))).hi
    low_sum = up_sum = ZERO
    crude_loss = ZERO
    nknot = 0
    worst_err = ZERO
    # cell 0: p~ = h^-gamma x^-gamma q0(x), q0 = 1 + h sum_k r_0k B_k x^(k+1)
    Sf0 = (1 << (WR + WB)) * N
    q0 = [(Sf0, 0)] + [(0, 0)] * (P + 1)
    err0 = ZERO
    for k, (rr, ri) in enumerate(R[0]):
        br, bi_, eb = Bint[k]
        q0[k + 1] = (rr * br - ri * bi_, rr * bi_ + ri * br)
        err0 = HI.add(err0, HI.multiply(HI.add(dint(rr, WR), dint(ri, WR)), eb))
    err0 = HI.divide(err0, Decimal(N))
    lo0, hi0, wint = gen_cell({(0, 0): q0}, Sf0, -a, 0, N, cs)
    low_sum = LO.subtract(lo0, HI.multiply(err0, wint))
    up_sum = HI.add(hi0, HI.multiply(err0, wint))
    for n in range(1, H):
        # A(xi) at scale 2^-WA
        tr, ti, et = fixz(zmul(cs.tbl[2 * n + 1], cs.two_g), W)
        Ar, Ai = [0] * (DA + 1), [0] * (DA + 1)
        err = ZERO
        tg_abs = HI.add(HI.add(dint(tr, W), dint(ti, W)), et)
        for i in range(DA + 1):
            br, bi_, eb = Bq[i]
            den = (2 * n + 1) ** i << (2 * W - WA)
            Ar[i] += _rdiv(tr * br - ti * bi_, den)
            Ai[i] += _rdiv(tr * bi_ + ti * br, den)
            e = HI.add(HI.multiply(tg_abs, eb), HI.multiply(up(Fraction(2 ** i)), et))
            err = HI.add(err, HI.add(HI.divide(e, Decimal((2 * n + 1) ** i)), HI.divide(ONE, two(WA))))
        err = HI.add(err, HI.multiply(tg_abs, up(Fraction(1, (2 * n + 1) ** (DA + 1)) / (1 - Fraction(1, 2 * n + 1)))))
        # far field
        if n >= 2:
            row = Prod[n]
            for qq in range(SF + 1):
                pr, pi = row[SF + qq]
                den = factorial(qq) * N << (WK + WM - WA)
                vr, vi = _rdiv(pr * 8 ** qq, den), _rdiv(pi * 8 ** qq, den)
                if qq <= DA:
                    Ar[qq] += vr
                    Ai[qq] += vi
                    err = HI.add(err, HI.divide(ONE, two(WA)))
                else:
                    err = HI.add(err, HI.divide(HI.add(dint(vr, WA + qq), dint(vi, WA + qq)), ONE))
                    err = HI.add(err, HI.divide(ONE, two(WA + qq - 1)))
            # tail of the s-series, and the rounding of K and M
            # sum_c c^-(SF+1)/(1 - 1/c) int_0^1 |rho~_(n-c)|, local norms for c <= 12, the largest beyond
            # (sum_{c >= 13} c^-(SF+1) / (1 - 1/c) <= 2 * 13^-SF)
            ft = ZERO
            for c in range(2, min(n, 12) + 1):
                t = HI.divide(l1[n - c], HI.multiply(LO.power(Decimal(c), SF + 1),
                                                     LO.subtract(ONE, HI.divide(ONE, Decimal(c)))))
                ft = HI.add(ft, t)
            if n > 12:
                ft = HI.add(ft, HI.multiply(l1max, HI.divide(Decimal(2), LO.power(Decimal(13), SF))))
            err = HI.add(err, HI.divide(HI.add(ft, far_round), Decimal(N)))
        # c = 1 extended
        Rp = R[n - 1]
        cr, ci = [0] * (DA + 1), [0] * (DA + 1)
        e1 = ZERO
        for k, (rr, ri) in enumerate(Rp):
            if not (rr or ri):
                continue
            g = Gk[k]
            for i in range(DA + 1):
                gr, gi = g[i]
                cr[i] += rr * gr - ri * gi
                ci[i] += rr * gi + ri * gr
            e1 = HI.add(e1, HI.multiply(HI.add(dint(rr, WR), dint(ri, WR)), eGk[k]))
        den = N << (WR + WG - WA)
        for i in range(DA + 1):
            Ar[i] += _rdiv(cr[i], den)
            Ai[i] += _rdiv(ci[i], den)
            err = HI.add(err, HI.divide(ONE, two(WA)))
        err = HI.add(err, HI.divide(e1, Decimal(N)))
        # singular part: jumps of rho~ at node n
        Rn = R[n]
        Dl = []
        sig = ZERO
        for k in range(P + 1):
            sr = Rn[k][0] - sum(Rp[j][0] * comb(j, k) for j in range(k, P + 1))
            si = Rn[k][1] - sum(Rp[j][1] * comb(j, k) for j in range(k, P + 1))
            Dl.append((sr, si))
            sig = HI.add(sig, HI.multiply(HI.add(dint(sr, WR), dint(si, WR)), Babs[k]))
        sig = HI.divide(sig, Decimal(N))
        worst_err = max(worst_err, err)
        A = list(zip(Ar, Ai))
        supw = up(Fraction(N, N - n - 1))
        if n in knots or HI.multiply(HI.multiply(sig, supw), hfac_hi) > KINK:
            nknot += 1
            # f = A(x' - 1/2) + h x'^(1-gamma) J(x'), J = sum_k Delta_k B_k x'^k
            Ax = []
            for jj in range(DA + 1):
                sr = si = 0
                for i in range(jj, DA + 1):
                    cbn = comb(i, jj) * (-1) ** (i - jj) << (DA - (i - jj))
                    sr += Ar[i] * cbn
                    si += Ai[i] * cbn
                Ax.append((sr, si))
            J = []
            eJ = ZERO
            for k, (dr, di) in enumerate(Dl):
                br, bi_, eb = Bint[k]
                J.append((dr * br - di * bi_, dr * bi_ + di * br))
                eJ = HI.add(eJ, HI.multiply(HI.add(dint(dr, WR), dint(di, WR)), eb))
            Sf = (1 << (WA + DA + WR + WB)) * N
            F = {(0, 0): [((x << (WR + WB)) * N, (y << (WR + WB)) * N) for x, y in Ax],
                 (1, 0): [(x << (WA + DA), y << (WA + DA)) for x, y in J]}
            lo_c, hi_c, wint = gen_cell(F, Sf, ZERO, n, N, cs)
            # a second upper bound (diagnostic side only), first order in S:
            # |A + S| <= |A| + Re(S conj A) u + eta |S| + |S|^2 u / (2 (1 - eta))
            _, hi2, supw2, U2, eta2, _ = reg_cell(A, WA, n, N)
            _, fo_hi = first_order(A, WA, U2, 40, J, WR + WB, n, N, inv1, W1, e1max)
            so = HI.divide(second_order(J, WR + WB, U2, 40, n, N, cs), LO.multiply(Decimal(2), LO.subtract(ONE, eta2)))
            hi_alt = HI.add(HI.add(HI.add(hi2, fo_hi), HI.multiply(HI.multiply(eta2, sig), supw2)), so)
            hi_c = min(hi_c, hi_alt)
            wint = max(wint, supw2)
            e_all = HI.add(err, HI.divide(eJ, Decimal(N)))
            low_sum = LO.add(low_sum, LO.subtract(lo_c, HI.multiply(e_all, wint)))
            up_sum = HI.add(up_sum, HI.add(hi_c, HI.multiply(e_all, wint)))
        else:
            # |A + S| >= |A| + Re(S conj A)/|A| >= g + Re(S conj A) u - eta |S|   (|1 - u|A|| <= eta), and
            # |A + S| <= |A| + Re(S conj A)/|A| + |S|^2/(2|A|) <= ... + eta |S| + |S|^2 u_max / (2 (1 - eta))
            lo_c, hi_c, supw, U, eta, umax = reg_cell(A, WA, n, N)
            J = []
            for k, (dr, di) in enumerate(Dl):
                br, bi_, eb = Bint[k]
                J.append((dr * br - di * bi_, dr * bi_ + di * br))
            fo_lo, fo_hi = first_order(A, WA, U, 40, J, WR + WB, n, N, inv1, W1, e1max)
            low_sum = LO.add(low_sum, LO.subtract(LO.add(lo_c, fo_lo),
                                                  HI.multiply(HI.add(err, HI.multiply(eta, sig)), supw)))
            sq = HI.divide(HI.multiply(HI.multiply(sig, sig), umax), LO.multiply(Decimal(2), LO.subtract(ONE, eta)))
            up_sum = HI.add(up_sum, HI.add(HI.add(hi_c, fo_hi), HI.multiply(HI.add(HI.add(err, HI.multiply(eta, sig)),
                                                                                    sq), supw)))
            crude_loss = HI.add(crude_loss, HI.multiply(HI.multiply(eta, sig), supw))
    den_lo = LO.multiply(low_sum, hfac_lo)
    den_hi = HI.multiply(up_sum, hfac_hi)
    out("  denominator: far-field convolution %.1f s, cells %.1f s (%d with kinks treated exactly); largest error in A "
        "%s; jump terms elsewhere exact to first order, remainder %s in total" % (t2 - t1, time.time() - t3, nknot, fmt(worst_err, 3),
                                                             fmt(HI.multiply(crude_loss, hfac_hi), 3)))
    out("  int_0^(1/2) |p~| dt/(1-t) in [%s, %s]" % (fmt(den_lo, 16), fmt(den_hi, 16)))
    return den_lo, den_hi


# ---- one row --------------------------------------------------------------------------------------------------------
def certify(prof, N, P, Q=None, g=None, out=print, label=""):
    """Returns a dict with the certified bounds, or the refusal."""
    t0 = time.time()
    ok, info = hypotheses(prof, N)
    T, C = prof["T"], prof["C"]
    out("%s: T = %s, C = %s, %d knots, grid N = %d, degree P = %d" % (label, prof["T"], prof["C"], len(prof["knots"]),
                                                                       N, P))
    for m in info.get("msgs", []):
        out("  " + m)
    if not ok:
        out("  REFUSED: a hypothesis of Theorem A fails or is not certified")
        return dict(verdict="REFUSED", msgs=info.get("msgs", []))
    out("  hypotheses: a = Re gamma in [%s, %s] > 0; max |psi_j - psi_0| = %s < pi/3; sup |d| on (0,1/2] <= %s; "
        "kappa = max(T, that) <= %s < 1" % (fmt(info["a"].lo, 10), fmt(info["a"].hi, 10), fmt(Decimal(
            info["E"].numerator) / Decimal(info["E"].denominator), 6), fmt(info["dmax"].hi, 6), fmt(info["kappa_hi"], 10)))
    bmax = max(abs(x) for x in info["slopes"]) / N
    if Q is None:
        Q = 4
        while T * bmax ** (Q + 1) / factorial(Q + 1) > Fraction(1, 2 ** 80):
            Q += 1
    if g is None:
        g = 1
        while N // (2 * g) > 1024 and (N // 2) % (2 * g) == 0:
            g *= 2
    WE = WR = 64
    E, eta, eta0, M = dtilde(prof, info, N, Q, WE)
    t1 = time.time()
    R = witness(prof, N, P, WR)
    out("  trial rho~: floating-point collocation, %.1f s; delta~ of degree Q = %d, |delta - delta~| <= %s"
        % (time.time() - t1, Q, fmt(max(eta), 3)))
    r = residual_bounds(N, P, Q, E, R, eta, eta0, WE, WR, out=out)
    _, EN, ED = majorant(prof, info, N, r, g, out=out)
    cs = Consts(info, N + 1, P, 26, 40)
    p_lo, p_hi = numerator(N, P, R, WR, cs, out=out)
    d_lo, d_hi = denominator(prof, info, N, P, R, WR, cs, M, out=out)
    num_hi = HI.add(p_hi, EN)
    den_lo = LO.subtract(d_lo, ED)
    assert den_lo > 0
    F_hi = HI.divide(num_hi, den_lo)
    F_lo = LO.divide(LO.subtract(p_lo, EN), HI.add(d_hi, ED))
    Cd = Decimal(C.numerator) / Decimal(C.denominator)      # C is a terminating decimal: exact
    passed = F_hi < Cd
    out("  F <= (|p~(1)| + EN) / (Den~ - ED) = %s / %s = %s" % (fmt(num_hi, 16), fmt(den_lo, 16), fmt(F_hi, 16)))
    out("  F >= %s (diagnostic lower bound); enclosure width %s" % (fmt(F_lo, 16), fmt(HI.subtract(F_hi, F_lo), 3)))
    out("  claim C = %s: C - F_upper = %s  ->  %s   (%.0f s)" % (C, fmt(HI.subtract(Cd, F_hi) if passed else
                                                                     LO.subtract(Cd, F_hi), 4),
                                                                 "PASS" if passed else "FAIL", time.time() - t0))
    return dict(verdict="PASS" if passed else "FAIL", F_hi=F_hi, F_lo=F_lo, C=C, T=T, EN=EN, ED=ED, p_hi=p_hi,
                d_lo=d_lo, maxr=max(r))


# ---- cross-checks -----------------------------------------------------------------------------------------------------
def harcos_closed(T, psi):
    """For a constant phase, rho = 0 on (0, 1/2] and rho = gamma/t on (1/2, 1], so
    F = |1 + gamma int_0^(1/2) t^-gamma/(1-t) dt| / int_0^(1/2) t^-a/(1-t) dt, with
    int_0^(1/2) t^-g/(1-t) dt = 2^g sum_k 2^-(k+1) / (k+1-g) (tail after K terms <= 2^(a-K) / (K+1-a))."""
    s, c = isincos_any(q(psi))
    a, b = q(T) * c, q(T) * s
    l2 = iln(q(2))
    mod = iexp(a * l2)
    s2, c2 = isincos_any(b * l2)
    two_g = Z(mod * c2, mod * s2)
    K = 170
    Sg, Sa = Z(I(0), I(0)), I(0)
    for k in range(K):
        w = q(Fraction(1, 2 ** (k + 1)))
        z = zinv_z(Z(I(k + 1) - a, -b))
        Sg = Sg + Z(z.re * w, z.im * w)
        Sa = Sa + w / (I(k + 1) - a)
    tail = HI.divide(HI.multiply(Decimal(2), Decimal(2) ** -K), LO.subtract(Decimal(K), Decimal(1)))
    Ig = zmul(two_g, Sg.widen(tail))
    Ia = (mod * Sa).widen(HI.multiply(Decimal(2), tail))
    num = zmul(Z(a, b), Ig)
    num = Z(num.re + I(1), num.im)
    nhi = zabs_hi(num)
    re_lo = num.re.lo if num.re.lo > 0 else (num.re.hi.copy_abs() if num.re.hi < 0 else ZERO)
    im_lo = num.im.lo if num.im.lo > 0 else (num.im.hi.copy_abs() if num.im.hi < 0 else ZERO)
    nlo = NEAR.next_minus(NEAR.sqrt(LO.add(LO.multiply(re_lo, re_lo), LO.multiply(im_lo, im_lo))))
    return LO.divide(nlo, Ia.hi), HI.divide(nhi, Ia.lo)


def finite_n(prof, n):
    """Floating point, not part of the proof: F_n = |P_n| / D_n of the exact construction for n (Lemma 2 with
    s_k = T e^(i phi(k/n)), k <= n/2), from m alpha_m = -sum_k s_k alpha_(m-k)."""
    T, ph, kn = float(prof["T"]), prof["phases"], prof["knots"]
    K = n // 2
    s = [0j]
    j = 0
    for k in range(1, K + 1):
        u = Fraction(k, n)
        while u > kn[j + 1]:
            j += 1
        phi = ph[j] + (ph[j + 1] - ph[j]) / (kn[j + 1] - kn[j]) * (u - kn[j])
        s.append(T * cmath.exp(1j * float(phi)))
    al = [1 + 0j]
    for m in range(1, n + 1):
        tot = 0j
        for k in range(1, min(m, K) + 1):
            tot += s[k] * al[m - k]
        al.append(-tot / m)
    P, acc = [], 0j
    for x in al:
        acc += x
        P.append(acc)
    D = sum(abs(P[n - k]) / k for k in range(K + 1, n + 1))
    return abs(P[n]) / D


# ---- the three rows, the self-tests and the red controls ---------------------------------------------------------
ROWS = {
    "bak": dict(key="[Bak26] Bakeas, PR #302", path=os.path.join("bakeas", "c42_profile_0.68637.json"),
                sha="dcea4edbec1d5935f2f04d4d1b5b8379e2d2d7d78072a1f089c7abbe11ded395", N=384, P=10,
                out="thmA_302_bak2026.txt",
                theirs="their two Arb checkers (N = 122880): F <= 0.6863660 and F <= 0.6863615; their float solve "
                       "F = 0.68636008"),
    "pin": dict(key="[Pin26b] Pinho, PR #305", path=os.path.join("pinho_vqk", "C42_compact_best.json"),
                sha="ac8d052d567069c6f2b32583791515628fde0fdfe3a1c15e8b2967d45e898765", N=9216, P=10,
                out="thmA_305_pin2026b.txt",
                theirs="their two Arb verifiers (N = 14155776, piecewise linear): F <= 0.6863607742 and "
                       "F <= 0.6863607773; Bakeas's float solve F = 0.686360744"),
    "vat": dict(key="[Vat26] Vatopoulos, PR #327", path=os.path.join("vatopoulos", "c42_profile_0.686360751.json"),
                sha="d015a8295d809fa26c9f18888ece53d77c9fea71712420313e75911065855652", N=18432, P=14,
                out="thmA_327_vat2026.txt",
                theirs="their checker #1 (N = 73728, P = 10): F in [0.6863607494206, 0.6863607494291]; checker #2 "
                       "(P = 12): F <= 0.6863607494249"),
}


class Tee:
    def __init__(self, path):
        self.f = open(path, "w") if path else None

    def __call__(self, s=""):
        print(s, flush=True)
        if self.f:
            self.f.write(s + "\n")
            self.f.flush()


def read_pinned(root, row):
    import hashlib
    raw = open(os.path.join(root, row["path"]), "rb").read()
    h = hashlib.sha256(raw).hexdigest()
    assert h == row["sha"], "SHA-256 mismatch for %s: %s" % (row["path"], h)
    tmp = json.loads(raw.decode("utf-8"))
    prof = dict(T=Fraction(tmp["T"]), phases=[Fraction(s) for s in tmp["phases"]], C=Fraction(tmp.get("C", tmp["T"])))
    if "breaks" in tmp:
        prof["knots"] = [Fraction(s) for s in tmp["breaks"]]
    else:
        assert tmp.get("knots") == "j/48, j = 0..24"
        prof["knots"] = [Fraction(j, 48) for j in range(25)]
    return prof, h


def run_row(root, name, outdir, finite=(1000, 4000, 16000)):
    row = ROWS[name]
    out = Tee(os.path.join(outdir, row["out"]) if outdir else None)
    prof, h = read_pinned(root, row)
    out("%s: C_42 <= %s, from Bakeas's Theorem A and the profile alone" % (row["key"], prof["C"].numerator /
                                                                            prof["C"].denominator))
    out("  profile %s, SHA-256 %s (as pinned)" % (row["path"], h))
    res = certify(prof, row["N"], row["P"], out=out, label="  " + row["key"])
    out("  for comparison only, read after our run: %s" % row["theirs"])
    if finite:
        vals = []
        for n in finite:
            vals.append("n = %d: %.7f" % (n, finite_n(prof, n)))
        out("  floating point, not part of the proof: F_n = |P_n| / D_n of the construction itself, " + ", ".join(vals))
    out("RESULT C42 <= %s: %s" % (prof["C"], res["verdict"]))
    return res


def conv_test(out, trials=6):
    """conv2d_complex against direct summation on random signed Gaussian-integer arrays of mixed sizes."""
    import random
    rng = random.Random(20261107)
    ok = True
    for t in range(trials):
        NX, NY = rng.randint(1, 40), rng.randint(1, 40)
        AX, BY = rng.randint(1, 8), rng.randint(1, 8)
        S = AX + BY - 1 + rng.randint(0, 3)
        big = 10 ** rng.randint(1, 60)
        X = [[(rng.randint(-big, big), rng.randint(-big, big)) for _ in range(AX)] for _ in range(NX)]
        Y = [[(rng.randint(-big, big), rng.randint(-big, big)) for _ in range(BY)] for _ in range(NY)]
        NZ = rng.randint(1, NX + NY - 1)
        Zc = conv2d_complex(X, Y, S, NZ)
        for n in range(NZ):
            for c in range(S):
                zr = zi = 0
                for i in range(max(0, n - NY + 1), min(n, NX - 1) + 1):
                    for a_ in range(AX):
                        b_ = c - a_
                        if 0 <= b_ < BY:
                            xr, xi = X[i][a_]
                            yr, yi = Y[n - i][b_]
                            zr += xr * yr - xi * yi
                            zi += xr * yi + xi * yr
                ok = ok and Zc[n][c] == (zr, zi)
    out("0. The exact two-dimensional convolution (decimal Kronecker products) against direct summation, %d random "
        "signed arrays with entries up to 10^60: %s" % (trials, "identical" if ok else "MISMATCH"))
    return ok


def selftest(root, outdir, quick=False):
    out = Tee(os.path.join(outdir, "thmA_selftest.txt") if outdir else None)
    out("Self-tests of thmA_check.py")
    conv_test(out)
    # 1. Harcos's constant phase, where F has a closed form
    T, psi = Fraction("0.69367634"), Fraction("0.89767188")
    lo, hi = harcos_closed(T, psi)
    out("1. Constant phase psi = %s, T = %s (gamma = 0.43246... + 0.54237... i, Harcos's case): rho = 0 on (0,1/2], "
        "rho = gamma/t on (1/2,1], closed form F in [%s, %s], so max(T, F) = %s, Harcos's 0.693676"
        % (psi, T, fmt(lo, 14), fmt(hi, 14), T))
    prof = dict(T=T, C=Fraction("0.6936764"), phases=[psi, psi], knots=[Fraction(0), Fraction(1, 2)])
    for N in ((16, 64) if not quick else (16,)):
        res = certify(prof, N, 10, out=out, label="   constant phase, N = %d" % N)
        out("   certified F_upper - closed-form F = %s (>= 0: %s); closed form inside [F_lower, F_upper]: %s"
            % (fmt(HI.subtract(res["F_hi"], lo), 3), res["F_hi"] >= hi, res["F_lo"] <= lo and hi <= res["F_hi"]))
    if not quick:
        convergence(root, out)
    return out


def convergence(root, out, runs=(("bak", 96, 10), ("bak", 192, 10), ("bak", 384, 10), ("bak", 768, 10),
                                  ("pin", 4608, 12), ("vat", 18432, 12))):
    out("2. Convergence: the same certificate at other grids N and degrees P (the rows use bak 384/10, pin 9216/10, "
        "vat 18432/14):")
    for name, N, P in runs:
        prof, _ = read_pinned(root, ROWS[name])
        t = time.time()
        res = certify(prof, N, P, out=lambda s: None, label="")
        out("   %s N = %5d, P = %2d: max |R/t| <= %s, EN <= %s, ED <= %s, F in [%s, %s]   (%.0f s)"
            % (ROWS[name]["key"].split()[0], N, P, fmt(res["maxr"], 3), fmt(res["EN"], 2), fmt(res["ED"], 2),
               fmt(res["F_lo"], 13), fmt(res["F_hi"], 13), time.time() - t))


def controls(root, outdir, which=("bak", "pin", "vat")):
    out = Tee(os.path.join(outdir, "thmA_controls.txt") if outdir else None)
    out("Red controls (each must FAIL or be REFUSED)")
    lowerT = dict(bak=Fraction("0.68636"), pin=Fraction("0.68636"), vat=Fraction("0.68636"))
    grids = dict(bak=(384, 10), pin=(4608, 12), vat=(18432, 12))
    perturb = dict(bak=(0, Fraction("0.05")), pin=(3, Fraction("0.02")), vat=(20, Fraction("0.005")))
    for name in which:
        prof, _ = read_pinned(root, ROWS[name])
        N, P = grids[name]
        p2 = dict(prof, T=lowerT[name], C=lowerT[name])
        res = certify(p2, N, P, out=out, label="%s with T = C = %s (below its F)" % (ROWS[name]["key"], lowerT[name]))
        out("   -> %s; F_lower - C = %s (> 0: the claim is false, not only unproved)"
            % (res["verdict"], fmt(LO.subtract(res["F_lo"], Decimal(p2["C"].numerator) / Decimal(p2["C"].denominator)),
                                   3)))
        j, dlt = perturb[name]
        ph = list(prof["phases"])
        ph[j] += dlt
        p3 = dict(prof, phases=ph)
        res = certify(p3, N, P, out=out, label="%s with psi_%d + %s" % (ROWS[name]["key"], j, dlt))
        out("   -> %s; F_lower - C = %s" % (res["verdict"], fmt(LO.subtract(res["F_lo"], Decimal(
            p3["C"].numerator) / Decimal(p3["C"].denominator)), 3)))
    prof, _ = read_pinned(root, ROWS["bak"])
    p4 = dict(prof, phases=[p - Fraction("0.3") for p in prof["phases"]])
    res = certify(p4, 192, 10, out=out, label="[Bak26] with every phase - 0.3 (psi_0 = -1.626, Re gamma < 0)")
    out("   -> %s" % res["verdict"])
    p5 = dict(prof, T=Fraction(1), C=Fraction(1))
    res = certify(p5, 192, 10, out=out, label="[Bak26] with T = C = 1 (kappa = 1)")
    out("   -> %s" % res["verdict"])
    prof, _ = read_pinned(root, ROWS["vat"])
    res = certify(prof, 9216, 10, out=out, label="[Vat26] on the grid N = 9216, which misses the knot 1/18432")
    out("   -> %s" % res["verdict"])


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    root = args[0] if args else os.path.join(HERE, "certs")
    rows = ["bak", "pin", "vat"]
    for f in flags:
        if f.startswith("--rows="):
            rows = [r for r in f.split("=", 1)[1].split(",") if r]
    outdir = os.path.join(HERE, "results") if "--write" in flags else None
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    if "--selftest" in flags:
        selftest(root, outdir)
    for name in rows:
        run_row(root, name, outdir)
    if "--controls" in flags:
        controls(root, outdir)


if __name__ == "__main__":
    main()
