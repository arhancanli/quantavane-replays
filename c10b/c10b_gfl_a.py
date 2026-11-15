#!/usr/bin/env python3
"""Program A for the 10b lower bound [GFL26] (Guo, Fang, Lu, arXiv 2609.07000v1, K_G^C > 1.35584631827168), recorded in
PR #196: the finite conditions of the paper's Theorem 3.1 (weighted criterion) for the exact parameters of its Section 4,
recomputed from the paper's definitions with our own code. Ball arithmetic on exact rationals: every quantity is a
Fraction centre with a Fraction radius, centres rounded to 2^-PREC with the rounding added to the radius.

Definitions used (paper, Sections 2-4):
  e_{a,b}(r) = (a! b!)^(-1/2) sum_{j <= min(a,b)} (-1)^j j! C(a,j) C(b,j) r^(a+b-2j)   (radial factor of h_{a,b});
  phi_l = e_{l+1,l};  d nu = e^(-s) ds with s = r^2;  Lambda = {1..10};
  I_{p,q} = {l in Lambda : p <= l+1, q <= l};
  M^{p,q}_{lm} = sqrt(beta_l beta_m) int e_{l+1-p,l-q} e_{m+1-p,m-q} / W(r^2) d nu   (need I - M^{p,q} > 0);
  D = 2 rho + W(r^2), R = int r^2/D, b_l = int r phi_l / D, G_lm = int phi_l phi_m / D, B = diag(beta),
  S = R - b^T (B^-1 + G)^-1 b  (need S <= 1);  mu_W = int W(r^2) d nu;  ratio = (1 - rho)/(rho + mu_W).
Route (ours): 1/W = Delta/N with Delta = (s + d1)(s + d2), N = Delta W (a cubic); its roots -lambda_i by exact rational
bisection; residues a_i = Delta(-lambda_i)/N'(-lambda_i); for a polynomial P, int P e^-s / (s + lambda) ds by synthetic
division P = (s + lambda) Q + P(-lambda): sum_k Q_k k! + P(-lambda) e^lambda E1(lambda); E1 by its power series with
Euler's gamma (Euler-Maclaurin with Bernoulli numbers) and ln (atanh series). Linear algebra: interval LDL without
pivoting, indices in increasing order (as the paper states), and S by solving (B^-1 + G) x = b by elimination.

  c10b_gfl_a.py [--set NAME VALUE]... [--claim-ratio DEC] [--claim-mu DEC]
--set replaces one parameter (rho, beta1..beta10, alpha, eta, zeta1, zeta2, delta1, delta2) by an exact decimal
(red controls). Standard library only; run with python3 -I.
"""
import sys
from fractions import Fraction as Fr
from math import comb, factorial, isqrt

PREC = 256
ONE = 1 << PREC

PARAMS = {  # Section 4.1 of the paper, every decimal read as an exact rational
    "rho": "0.012039377730497942",
    "beta1": "0.4283744445791419", "beta2": "0.1506077793447579", "beta3": "0.08328691492582467",
    "beta4": "0.05294104124824108", "beta5": "0.03632381026351961", "beta6": "0.026165030774068086",
    "beta7": "0.01945831930582338", "beta8": "0.01490133597851163", "beta9": "0.011682554694871422",
    "beta10": "0.014696581084273385",
    "alpha": "0.16212090047056932029265504", "eta": "0.2995274614735440217489782",
    "zeta1": "0.1705006229073915451327806", "delta1": "0.6122395153199937",
    "zeta2": "1.3819822296622917537966546", "delta2": "6.340370329139007",
}


class B:
    """A ball: centre c (Fraction), radius r (Fraction >= 0)."""
    __slots__ = ("c", "r")

    def __init__(self, c, r=0):
        self.c, self.r = Fr(c), Fr(r)
        self._round()

    def _round(self):
        c = self.c
        if c.denominator > ONE:
            cc = Fr((c.numerator * ONE) // c.denominator, ONE)
            self.r += c - cc
            self.c = cc
        r = self.r
        if r.denominator > ONE:                     # the radius is rounded up to the same grid
            self.r = Fr(-((-r.numerator * ONE) // r.denominator), ONE)

    def __add__(s, o):
        o = o if isinstance(o, B) else B(o)
        return B(s.c + o.c, s.r + o.r)
    __radd__ = __add__

    def __neg__(s):
        return B(-s.c, s.r)

    def __sub__(s, o):
        o = o if isinstance(o, B) else B(o)
        return B(s.c - o.c, s.r + o.r)

    def __rsub__(s, o):
        return B(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, B) else B(o)
        return B(s.c * o.c, abs(s.c) * o.r + abs(o.c) * s.r + s.r * o.r)
    __rmul__ = __mul__

    def inv(s):
        m = abs(s.c) - s.r
        assert m > 0, "inverse of a ball containing 0"
        return B(1 / s.c, s.r / (abs(s.c) * m))

    def __truediv__(s, o):
        o = o if isinstance(o, B) else B(o)
        return s * o.inv()

    def lo(s):
        return s.c - s.r

    def hi(s):
        return s.c + s.r


def b_sqrt(x):
    lo, hi = x.lo(), x.hi()
    assert lo > 0
    a = Fr(isqrt((lo.numerator * ONE * ONE) // lo.denominator), ONE)            # <= sqrt(lo)
    t = isqrt(-((-hi.numerator * ONE * ONE) // hi.denominator))
    bb = Fr(t + 1, ONE)                                                       # >= sqrt(hi)
    return B((a + bb) / 2, (bb - a) / 2)


def atanh_series(t, eps):
    """atanh(t) for a rational |t| <= 1/2: partial sum plus a geometric tail bound."""
    s, k, p = Fr(0), 0, t
    while True:
        term = p / (2 * k + 1)
        tail = abs(term) * Fr(4, 3)
        if tail < eps:
            return B(s, tail)
        s += term
        p *= t * t
        k += 1


EPS = Fr(1, 1 << (PREC - 10))
LN2 = atanh_series(Fr(1, 3), EPS) * 2
LN5 = atanh_series(Fr(1, 3), EPS) * 4 + atanh_series(Fr(1, 9), EPS) * 2   # ln 5 = 2 ln 2 + ln(5/4)


def ln_ball(x):
    """ln of a positive ball: ln(centre) by range reduction, plus radius / (centre - radius)."""
    c = x.c
    e = 0
    while c > Fr(3, 2):
        c /= 2; e += 1
    while c < Fr(3, 4):
        c *= 2; e -= 1
    v = atanh_series((c - 1) / (c + 1), EPS) * 2 + LN2 * e
    return v + B(0, x.r / (x.c - x.r))


def bernoulli(nmax):
    """B_0..B_nmax (B_1 = +1/2 convention irrelevant here: only even indices are used), Akiyama-Tanigawa."""
    out, a = [], [Fr(0)] * (nmax + 1)
    for m in range(nmax + 1):
        a[m] = Fr(1, m + 1)
        for j in range(m, 0, -1):
            a[j - 1] = j * (a[j - 1] - a[j])
        out.append(a[0])
    return out


def euler_gamma():
    n, K = 200, 30
    Bn = bernoulli(2 * K + 2)
    H = sum(Fr(1, i) for i in range(1, n + 1))
    s = H - Fr(1, 2 * n) + sum(Bn[2 * k] / (2 * k * Fr(n) ** (2 * k)) for k in range(1, K + 1))
    rem = abs(Bn[2 * K + 2]) / ((2 * K + 2) * Fr(n) ** (2 * K + 2))
    lnn = LN2 * 3 + LN5 * 2                                                  # 200 = 2^3 5^2
    return B(s, rem) - lnn


GAMMA = euler_gamma()


def exp_ball(x):
    c = x.c
    assert 0 <= c < 40
    s, t, k = Fr(0), Fr(1), 0
    while True:
        if k > 2 * c and t < EPS:
            break
        s += t
        k += 1
        t = t * c / k
    v = B(s, 2 * t)
    bound = Fr(3) ** (int(c + x.r) + 1)                                     # e^(c + r) < 3^(ceil)
    return v + B(0, x.r * bound)


def e1_ball(x):
    """E1 of a positive ball by the series -gamma - ln x + sum_{k>=1} (-1)^(k+1) x^k/(k k!)."""
    c = x.c
    s, k, p = Fr(0), 1, c          # p = c^k / k!
    while True:
        term = p / k
        if k > 2 * c + 2 and term < EPS:
            break
        s += term if k % 2 == 1 else -term
        k += 1
        p = p * c / k
    v = B(s, term) - GAMMA - ln_ball(B(c))
    return v + B(0, x.r / (x.c - x.r))          # |E1'| = e^-x / x <= 1/x


def poly_eval(coefs, x):
    acc = B(0)
    for a in reversed(coefs):
        acc = acc * x + a
    return acc


def poly_mul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def roots_cubic(N, brackets):
    """Roots of the rational cubic N in the given rational brackets (sign change checked), as balls."""
    f = lambda x: sum(a * x ** i for i, a in enumerate(N))
    out = []
    for lo, hi in brackets:
        flo, fhi = f(lo), f(hi)
        assert flo * fhi < 0, "no sign change in a bracket"
        while hi - lo > Fr(1, 1 << (PREC - 20)):
            mid = Fr(((lo + hi) / 2 * ONE).__floor__(), ONE)
            fm = f(mid)
            if fm == 0:
                lo = hi = mid; break
            if (fm < 0) == (flo < 0):
                lo, flo = mid, fm
            else:
                hi = mid
        out.append(B((lo + hi) / 2, (hi - lo) / 2))
    return out


def setup(p):
    d1, d2 = p["delta1"], p["delta2"]
    Delta = poly_mul([d1, 1], [d2, 1])
    s_d2 = poly_mul([0, 1], [d2, 1])
    s_d1 = poly_mul([0, 1], [d1, 1])
    N = [Fr(0)] * 4
    for i, a in enumerate(poly_mul([p["alpha"], p["eta"]], Delta)):
        N[i] += a
    for i, a in enumerate(s_d2):
        N[i] += p["zeta1"] * a
    for i, a in enumerate(s_d1):
        N[i] += p["zeta2"] * a
    ND = [N[i] + 2 * p["rho"] * (Delta[i] if i < 3 else 0) for i in range(4)]
    return Delta, N, ND


def partial_fractions(Delta, N, d1, d2):
    """Roots -lambda_i of N in (-inf,-d2), (-d2,-d1), (-d1,0) and residues a_i: Delta/N = sum a_i/(s + lambda_i)."""
    big = Fr(-1000)
    rts = roots_cubic(N, [(big, -d2), (-d2, -d1), (-d1, Fr(0))])
    dN = [i * N[i] for i in range(1, 4)]
    out = []
    for r in rts:
        lam = -r
        a = poly_eval(Delta, r) / poly_eval(dN, r)
        Jl = exp_ball(lam) * e1_ball(lam)          # int e^-s/(s + lambda) ds = e^lambda E1(lambda)
        out.append((lam, a, Jl))
    return out


def integral(P, pf):
    """int_0^inf P(s) e^-s Delta/N ds for a polynomial P with rational coefficients (ascending)."""
    total = B(0)
    for lam, a, Jl in pf:
        # synthetic division of P by (s + lambda): P = (s + lambda) Q + P(-lambda)
        n = len(P) - 1
        q = [None] * n
        acc = B(P[n])
        for k in range(n - 1, -1, -1):
            q[k] = acc
            acc = acc * (-lam) + P[k]
        val = acc * Jl
        for k in range(n):
            val = val + q[k] * factorial(k)
        total = total + a * val
    return total


def radial(a, b):
    """(a! b!)^(1/2) e_{a,b}(r) as a polynomial in r (ascending coefficients, integers)."""
    out = [0] * (a + b + 1)
    for j in range(min(a, b) + 1):
        out[a + b - 2 * j] += (-1) ** j * factorial(j) * comb(a, j) * comb(b, j)
    return out


def r_to_s(poly_r):
    assert all(c == 0 for c in poly_r[1::2]), "odd powers of r in an even product"
    return [Fr(c) for c in poly_r[0::2]]


def prod_r(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return out


def ldl_pivots(A):
    """Pivots of an interval LDL^T decomposition without pivoting (A symmetric, list of lists of balls)."""
    n = len(A)
    L = [[None] * n for _ in range(n)]
    d = []
    for j in range(n):
        dj = A[j][j]
        for k in range(j):
            dj = dj - L[j][k] * L[j][k] * d[k]
        d.append(dj)
        if dj.lo() <= 0:
            return d, False
        for i in range(j + 1, n):
            v = A[i][j]
            for k in range(j):
                v = v - L[i][k] * L[j][k] * d[k]
            L[i][j] = v / dj
    return d, True


def solve(A, bvec):
    """Gaussian elimination without pivoting on balls (A positive definite)."""
    n = len(A)
    M = [row[:] + [bvec[i]] for i, row in enumerate(A)]
    for j in range(n):
        for i in range(j + 1, n):
            f = M[i][j] / M[j][j]
            for k in range(j, n + 1):
                M[i][k] = M[i][k] - f * M[j][k]
    x = [None] * n
    for i in range(n - 1, -1, -1):
        v = M[i][n]
        for k in range(i + 1, n):
            v = v - M[i][k] * x[k]
        x[i] = v / M[i][i]
    return x


def fmt(x, d=25):
    lo, hi = x.lo(), x.hi()
    f = lambda v, up: (-((-v.numerator * 10 ** d) // v.denominator) if up else (v.numerator * 10 ** d) // v.denominator)
    def s(q):
        neg = q < 0
        t = str(abs(q)).rjust(d + 1, "0")
        return ("-" if neg else "") + t[:-d] + "." + t[-d:]
    return "[%s, %s]" % (s(f(lo, False)), s(f(hi, True)))


def main():
    p = {k: Fr(v) for k, v in PARAMS.items()}
    claim_ratio = claim_mu = None
    a = sys.argv[1:]
    while a:
        if a[0] == "--set":
            p[a[1]] = Fr(a[2]); print("CONTROL %s set to %s" % (a[1], a[2])); a = a[3:]
        elif a[0] == "--claim-ratio":
            claim_ratio = Fr(a[1]); a = a[2:]
        elif a[0] == "--claim-mu":
            claim_mu = Fr(a[1]); a = a[2:]
        else:
            raise SystemExit("unknown option " + a[0])
    Lam = list(range(1, 11))
    beta = {l: p["beta%d" % l] for l in Lam}
    rho = p["rho"]
    assert all(v > 0 for v in p.values()) and 0 < rho < 1 and p["delta1"] < p["delta2"]
    print("Euler's gamma in %s" % fmt(GAMMA, 40))
    # self-test of the radial factors: orthonormal in L2(nu) within each angular frequency
    for a_ in range(0, 7):
        for b_ in range(0, 7):
            for c_ in range(0, 7):
                d_ = c_ - (a_ - b_)
                if d_ < 0:
                    continue
                pr = prod_r(radial(a_, b_), radial(c_, d_))
                val = sum(c * factorial(k) for k, c in enumerate(r_to_s(pr)))
                expect = factorial(a_) * factorial(b_) if (a_, b_) == (c_, d_) else 0
                assert val == expect, (a_, b_, c_, d_)
    print("self-test: the radial factors e_{a,b} (a, b <= 6) are orthonormal in L2(nu), exactly")
    Delta, N, ND = setup(p)
    pfW = partial_fractions(Delta, N, p["delta1"], p["delta2"])
    pfD = partial_fractions(Delta, ND, p["delta1"], p["delta2"])
    for name, pf in (("W", pfW), ("2 rho + W", pfD)):
        print("1/(%s): poles -lambda with lambda in %s" % (name, ", ".join(fmt(l, 12) for l, _, _ in pf)))
    chk = integral([Fr(1)], pfW)
    print("int e^-s / W(s) ds in %s" % fmt(chk))
    # Gram blocks
    minpiv, minblk, npiv, nblk, allok = None, None, 0, 0, True
    for pp in range(0, 12):
        for qq in range(0, 11):
            I = [l for l in Lam if pp <= l + 1 and qq <= l]
            if not I:
                continue
            nblk += 1
            rad = {l: radial(l + 1 - pp, l - qq) for l in I}
            fac = {l: factorial(l + 1 - pp) * factorial(l - qq) for l in I}
            A = []
            for l in I:
                row = []
                for m in I:
                    Pint = integral(r_to_s(prod_r(rad[l], rad[m])), pfW)
                    scale = b_sqrt(B(beta[l] * beta[m] / (fac[l] * fac[m])))
                    v = scale * Pint
                    row.append((B(1) - v) if l == m else -v)
                A.append(row)
            d, ok = ldl_pivots(A)
            npiv += len(d)
            allok &= ok
            for v in d:
                if minpiv is None or v.lo() < minpiv.lo():
                    minpiv, minblk = v, (pp, qq)
    print("Gram blocks: %d nonempty, %d pivots, all positive: %s; least pivot in %s, block (p,q) = %s"
          % (nblk, npiv, allok, fmt(minpiv, 15), minblk))
    # Schur condition
    phi = {l: radial(l + 1, l) for l in Lam}
    nphi = {l: factorial(l + 1) * factorial(l) for l in Lam}
    R = integral([Fr(0), Fr(1)], pfD)
    bv = [integral(r_to_s(prod_r([0, 1], phi[l])), pfD) * b_sqrt(B(Fr(1, nphi[l]))) for l in Lam]
    G = [[integral(r_to_s(prod_r(phi[l], phi[m])), pfD) * b_sqrt(B(Fr(1, nphi[l] * nphi[m]))) for m in Lam] for l in Lam]
    A = [[G[i][j] + (B(1 / beta[Lam[i]]) if i == j else B(0)) for j in range(10)] for i in range(10)]
    dA, okA = ldl_pivots(A)
    x = solve(A, bv)
    S = R
    for i in range(10):
        S = S - bv[i] * x[i]
    print("B^-1 + G positive definite (LDL): %s" % okA)
    print("R in %s" % fmt(R))
    print("S in %s" % fmt(S))
    print("1 - S in %s" % fmt(B(1) - S, 20))
    # weight mean and ratio
    mu = B(p["alpha"] + p["eta"])
    for j in (1, 2):
        dl = B(p["delta%d" % j])
        mu = mu + p["zeta%d" % j] * (B(1) - dl * exp_ball(dl) * e1_ball(dl))
    ratio = (B(1) - rho) / (B(rho) + mu)
    print("mu_W in %s" % fmt(mu))
    print("ratio (1 - rho)/(rho + mu_W) in %s" % fmt(ratio))
    ok = allok and okA and (B(1) - S).lo() > 0
    print("conditions of Theorem 3.1 (all Gram blocks, S < 1): %s" % ("hold" if ok else "FAIL"))
    if claim_mu is not None:
        h = mu.hi() < claim_mu
        print("claimed mu_W < %s: %s" % (claim_mu, "holds" if h else "FAILS")); ok &= h
    if claim_ratio is not None:
        h = ratio.lo() > claim_ratio
        print("claimed ratio > %s: %s" % (claim_ratio, "holds" if h else "FAILS")); ok &= h
    print("VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
