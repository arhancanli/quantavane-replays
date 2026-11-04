#!/usr/bin/env python3
"""Independent replay of the asterisked upper bounds for C_86, the Schur-Siegel-Smyth trace constant
(teorth/optimizationproblems, constant 86a), from the certificates' data files.

Written from the definition and the criterion as the 86a page states them, and from the certificate format notes
(FORMAT.md / PROOF.md of [Kos2026], README of [Pin26]; for PR #298 the README, PROOF.md, METHOD.md and
VERIFIER_NOTES.md of c86/adaptive-cone). The claimants' checkers and audit scripts were neither read nor run.
Standard library only (Python 3.9+); run with python3 -I. No floating-point number decides anything: floats only
choose where to split, where to put a tangent point and where to start a root iteration; every decision is an exact
integer or rational comparison.

  c86_check.py ROW CERTS_ROOT [--controls] [--workers=N]     ROW in: kos pin vid p298

The certificate: edges e_0 < ... < e_N, densities d_i >= 0 on [e_i, e_{i+1}], integer polynomials Q_j (constant
term first), weights b_j >= 0, and `bound`. With mu = sum_i d_i 1_[e_i,e_{i+1}] dx, m = mu(R), T = int x dmu,
S = supp mu (the union of the cells with d_i > 0), U(z) = int log|z - y| dmu(y) and f = U - sum_j b_j log|Q_j|:
  E1 edges strictly increasing, every d_i >= 0 and b_j >= 0, m > 0, min S > 0;
  E2 every Q_j primitive, of degree >= 1, irreducible over Z, no two associate; and all roots real (the page states
     it for both rows; the argument does not need it);
  E3 sum_j b_j deg Q_j < m;
  E4 int log|Q_j| dmu >= 0 for every j (the page says > 0);
  E5 f(x) >= 0 for every x in S (and >= the margin a row claims);
  E6 T/m equals `bound` exactly, and the printed decimal is >= T/m.

How each is decided
- Arithmetic. Exact Fractions; transcendental quantities are integer intervals [lo, hi] at scale 2^-128. Logarithms
  are our own: ln n = e ln 2 + ln(1 + j/256) + 2 atanh(z), z < 1/512, every series term rounded outward and the tail
  bounded; tested against the correctly rounded `decimal` ln.
- E2. Primitivity by gcd. Irreducibility two ways: (a) distinct-degree factorisation modulo primes; the degrees of
  any factorisation over Z must be subset sums of the factor degrees modulo every good prime, and when only 0 and n
  survive Q is irreducible; (b) for monic real-rooted Q, every subset of at most n/2 isolated roots: some elementary
  symmetric function encloses no integer, or the candidate factor is tested by exact division. Real roots: exact
  Sturm counts; every root is isolated and refined by exact sign bisection to width 2^-112.
- E4. int log|Q| dmu = m log|lc Q| + sum_roots U(r), U(x) = sum_k c_k G(x - e_k), G(u) = u log|u| - u,
  c_k = d_k - d_{k-1}: each G(r - e_k) is enclosed at the two ends of r's enclosure (G is monotone there).
- E5. On each support cell [e_i, e_{i+1}], f = S + (the two edge terms c_i G(x - e_i), c_{i+1} G(x - e_{i+1}))
  - (b log|x - r| for roots r within H of the cell). S, the rest, is analytic on the cell and is replaced by its
  Taylor polynomial of degree 28 at the cell midpoint: the coefficients are exact rational sums (edge part) and
  outward-rounded root sums (root part), with explicit bounds for the remainder of S, S' and S''. A piece J of
  the cell is certified by the larger of two lower bounds, each the minimum of a concave function, hence attained
  at an end of J:
    A  chord of S, minus (|J|^2/8) sup_J S'' (when positive), plus the tangent group;
    B  S(t) + S'(t)(x - t) + min(0, inf_J S'')(x - t)^2 / 2 plus the tangent group.
  The tangent group is the convex terms (an edge term convex on the cell; a nearby root outside J), each replaced
  by its tangent at the same point t of J; concave edge terms are replaced by their chords; a root inside J by
  -b log(max distance). Pieces are bisected until every bound reaches the goal. Separate checks: the Taylor model
  agrees with direct termwise sums at random points; f >= 0 at gap and exterior points (the minimum principle);
  the energy I(mu) = int U dmu closed form, and int f dmu = I(mu) - sum_j b_j E4_j >= m min_S f.

PR #298 (row p298) has a different format: node heights h_k, the density linear between the nodes e_k (5,375 of
them, unequal spacing) and zero outside; auxiliaries may have nonreal roots (E2 there: primitive, irreducible by the
sieve, positive leading coefficient, distinct). U(z) = Re sum_k j_k F(z - e_k), j_k the slope jumps,
F(w) = w^2 log(w)/2 - 3w^2/4 (principal logarithm; G2(u) = u^2 log|u|/2 - 3u^2/4 on the real line).
- Roots: real ones by exact Sturm isolation; every root also in a Gershgorin disc of the Weierstrass matrix (radius
  about 2^-305; the code of our C82 replay), the discs off the real axis holding exactly the nonreal roots; Vieta's
  e1, e2 and |e_n| are checked against the enclosures.
- E4 from the closed form at every root: log by our own series, arg by our own arctangent (table of atan(j/256) by
  Euler's series, then the alternating series below 1/512; pi by Machin's formula); a real root costs
  5 rad sum|j_k| (Lipschitz), a nonreal one rad m/(Im c - rad).
- E5 per support cell, x = x0 + r sigma: nodes and roots within 4r are exact on each piece (chords for concave,
  tangents for convex terms; every near node is within distance 1, where G2 is concave), the rest is a Taylor
  polynomial of degree 24 in sigma (power sums of r/(x0 - e_k) from exact integers, table logarithms) with explicit
  remainder bounds for the value and two derivatives; bounds A and B as above, bisection until each piece clears
  the goal (0.9 x the cell's float minimum), and the cells near the global minimum again at (1 - 1e-8) x it.

What rests on mathematics rather than computation (cited, not re-derived): the minimum principle (f >= 0 on S gives
f >= 0 on C), Smith's Prop. 5.15 (E2-E5 make mu/m arithmetic) and Thm 1.5 (arithmetic measures on an interval of
capacity > 1 are limits of conjugate distributions), as [Kos2026] PROOF.md, the [Pin26] README and the #298
PROOF.md state them.
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
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import isqrt
from operator import mul

# ----------------------------------------------------------------------------------------------------------------
# Fixed-point intervals: an integer v stands for v / 2^P; an interval is a pair (lo, hi)
# ----------------------------------------------------------------------------------------------------------------

P = 128
ONE = 1 << P
LP = 168                       # internal scale of the logarithms and the Taylor sums
SH = LP - P
NT = 28                        # Taylor degree of the per-cell model


def cdiv(a, b):
    """ceil(a / b) for b > 0."""
    return -((-a) // b)


def ivq(iv, p, q):
    """The interval iv times the exact rational p/q (q > 0), rounded outward."""
    lo, hi = iv
    if p >= 0:
        return lo * p // q, cdiv(hi * p, q)
    return hi * p // q, cdiv(lo * p, q)


def qiv(fr, scale=P):
    """An exact rational as an interval at the given scale."""
    n, d = fr.numerator, fr.denominator
    return (n << scale) // d, cdiv(n << scale, d)


def ivfl(iv, scale=P):
    return (iv[0] + iv[1]) / 2.0 / float(1 << scale)


def fx(v, scale=P):
    """A fixed-point integer as a float (display only)."""
    return v / float(1 << scale)


def frq(fr, p, q):
    """floor(fr * p / q) for a Fraction fr (helper for rational times integer)."""
    return (fr.numerator * p) // (fr.denominator * q)


# ----------------------------------------------------------------------------------------------------------------
# Rigorous logarithms
# ----------------------------------------------------------------------------------------------------------------

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
    return lo, hi + 2 * t + 2           # remaining terms sum to at most (4/3) t / k


_LN2 = tuple(2 * v for v in _atanh(1, 3))
_TAB = [tuple(2 * v for v in _atanh(j, 512 + j)) for j in range(256)]   # ln(1 + j/256)


def ln_int_lp(n):
    """ln n for an integer n >= 1, as an interval at scale 2^-LP."""
    e = n.bit_length() - 1
    j = ((n << 8) >> e) - 256
    a = (n << 8) - ((256 + j) << e)
    b = (n << 8) + ((256 + j) << e)
    zl, zh = _atanh(a, b)
    return e * _LN2[0] + _TAB[j][0] + 2 * zl, e * _LN2[1] + _TAB[j][1] + 2 * zh


@lru_cache(maxsize=1 << 21)
def ln_int(n):
    lo, hi = ln_int_lp(n)
    return lo >> SH, cdiv(hi, 1 << SH)


def ln_q(a, b):
    """ln(a/b) for positive integers a, b, at scale 2^-P."""
    x = ln_int(a)
    y = ln_int(b)
    return x[0] - y[1], x[1] - y[0]


def ln_q_lp(a, b):
    x = ln_int_lp(a)
    y = ln_int_lp(b)
    return x[0] - y[1], x[1] - y[0]


def G_rat(num, den):
    """G(u) = u ln|u| - u at u = num/den (den > 0), G(0) = 0."""
    if num == 0:
        return 0, 0
    L = ln_q(abs(num), den)
    return ivq((L[0] - ONE, L[1] - ONE), num, den)


# ----------------------------------------------------------------------------------------------------------------
# Reading a certificate (as data)
# ----------------------------------------------------------------------------------------------------------------

RAT = re.compile(r"-?[0-9]+(/[0-9]+)?$")


def rat(v):
    if isinstance(v, bool):
        raise ValueError("boolean where a rational is expected")
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str) and RAT.match(v):
        return Fraction(v)
    raise ValueError("not an exact rational: %r" % (v,))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_cert(path):
    d = json.load(open(path))
    for key in ("edges", "densities", "A", "b", "bound"):
        if key not in d:
            raise ValueError("missing key " + key)
    edges = [rat(x) for x in d["edges"]]
    dens = [rat(x) for x in d["densities"]]
    A = []
    for Q in d["A"]:
        if not (isinstance(Q, list) and Q and all(isinstance(c, int) and not isinstance(c, bool) for c in Q)):
            raise ValueError("a polynomial is not a list of integers: %r" % (Q,))
        A.append(list(Q))
    b = [rat(x) for x in d["b"]]
    return {"edges": edges, "dens": dens, "A": A, "b": b, "bound": rat(d["bound"]), "raw": d}


def pstr(Q):
    """x^2 - 3x + 1 style."""
    out = []
    for i in range(len(Q) - 1, -1, -1):
        c = Q[i]
        if c == 0:
            continue
        mono = "" if i == 0 else ("x" if i == 1 else "x^%d" % i)
        a = abs(c)
        body = (str(a) if (a != 1 or i == 0) else "") + mono
        out.append(("- " if c < 0 else "+ ") + body)
    s = " ".join(out)
    return s[2:] if s.startswith("+ ") else "-" + s[2:]


# ----------------------------------------------------------------------------------------------------------------
# Integer polynomials: evaluation, Sturm sequences, real roots, irreducibility
# ----------------------------------------------------------------------------------------------------------------

def ev_int(Q, p, q):
    """q^d Q(p/q) as an exact integer (q > 0, so it has the sign of Q(p/q))."""
    d = len(Q) - 1
    acc, qq = Q[d], 1
    for m in range(d - 1, -1, -1):
        qq *= q
        acc = acc * p + Q[m] * qq
    return acc


def ev_frac(Q, x):
    v = Fraction(0)
    for c in reversed(Q):
        v = v * x + c
    return v


def p_trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def p_rem(a, b):
    a = list(a)
    db, lb = len(b) - 1, b[-1]
    while a and len(a) - 1 >= db:
        if a[-1] == 0:
            a.pop()
            continue
        f, s = a[-1] / lb, len(a) - 1 - db
        for i in range(db + 1):
            a[s + i] -= f * b[i]
        assert a[-1] == 0
        a.pop()
    return p_trim(a)


def sturm_seq(Q):
    s0 = [Fraction(c) for c in Q]
    s1 = [Fraction(i * Q[i]) for i in range(1, len(Q))]
    seq = [s0, s1]
    while len(seq[-1]) > 1:
        r = p_rem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])
    return seq


def _sgn(v):
    return (v > 0) - (v < 0)


def _var(signs):
    s = [x for x in signs if x]
    return sum(1 for a, b in zip(s, s[1:]) if a != b)


def var_at(seq, x):
    return _var([_sgn(ev_frac(p, x)) for p in seq])


def var_inf(seq, positive):
    return _var([_sgn(p[-1]) * (1 if positive or (len(p) - 1) % 2 == 0 else -1) for p in seq])


def real_root_count(Q):
    """(number of distinct real roots, squarefree?) by Sturm's theorem."""
    seq = sturm_seq(Q)
    squarefree = len(seq[-1]) == 1
    return var_inf(seq, False) - var_inf(seq, True), squarefree


def isolate_roots(Q, bits=112):
    """All real roots of a squarefree integer Q as dyadic enclosures (lo, hi), hi - lo <= 2^-bits, Q(lo)Q(hi) < 0;
    (r, r) for the exact root of a linear Q."""
    d = len(Q) - 1
    if d == 1:
        r = Fraction(-Q[0], Q[1])
        return [(r, r)]
    seq = sturm_seq(Q)
    bound = 1 + max(abs(Fraction(c, Q[-1])) for c in Q[:-1])
    B = 1
    while B < bound:
        B *= 2
    out, stack = [], [(Fraction(-B), Fraction(B))]
    while stack:
        a, b = stack.pop()
        n = var_at(seq, a) - var_at(seq, b)
        if n == 0:
            continue
        sa, sb = _sgn(ev_frac(Q, a)), _sgn(ev_frac(Q, b))
        if n == 1 and sa * sb < 0:
            out.append((a, b))
            continue
        mid = (a + b) / 2
        if ev_frac(Q, mid) == 0:
            raise ValueError("rational root %s of a polynomial of degree %d" % (mid, d))
        stack += [(a, mid), (mid, b)]
    res = []
    for a, b in out:
        K = max(bits + 8, (b - a).denominator.bit_length())
        A, Bn = int(a * 2 ** K), int(b * 2 ** K)          # dyadic, exact
        sA = _sgn(ev_int(Q, A, 1 << K))
        while Bn - A > (1 << (K - bits)):
            M = (A + Bn) // 2
            sM = _sgn(ev_int(Q, M, 1 << K))
            assert sM != 0
            if sM == sA:
                A = M
            else:
                Bn = M
        res.append((Fraction(A, 1 << K), Fraction(Bn, 1 << K)))
    return sorted(res)


# ---- arithmetic modulo p (constant term first) ----

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


def fp_mulmod(a, b, f, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return fp_divmod(fp_trim(r, p), f, p)[1]


def fp_powmod(base, e, f, p):
    result, b = [1], fp_divmod(base, f, p)[1]
    while e:
        if e & 1:
            result = fp_mulmod(result, b, f, p)
        b = fp_mulmod(b, b, f, p)
        e >>= 1
    return result


def fp_gcd(a, b, p):
    a, b = fp_trim(a, p), fp_trim(b, p)
    while b:
        a, b = b, fp_divmod(a, b, p)[1]
    if a:
        inv = pow(a[-1], p - 2, p)
        a = [c * inv % p for c in a]
    return a


def ddf_degrees(Q, p):
    """Degrees of the irreducible factors of Q modulo p (distinct-degree factorisation), or None when p divides the
    leading coefficient or Q is not squarefree modulo p."""
    f = fp_trim(Q, p)
    if len(f) != len(Q):
        return None
    inv = pow(f[-1], p - 2, p)
    f = [c * inv % p for c in f]
    df = fp_trim([i * f[i] for i in range(1, len(f))], p)
    if not df or len(fp_gcd(f, df, p)) > 1:
        return None
    degs, g, h, d = [], f, [0, 1], 0
    while 2 * (d + 1) <= len(g) - 1:
        d += 1
        h = fp_powmod(h, p, g, p)
        hx = list(h) + [0] * max(0, 2 - len(h))
        hx[1] = (hx[1] - 1) % p
        fac = fp_gcd(g, hx, p)
        if len(fac) > 1:
            degs += [d] * ((len(fac) - 1) // d)
            g = fp_divmod(g, fac, p)[0]
            h = fp_divmod(h, g, p)[1]
    if len(g) > 1:
        degs.append(len(g) - 1)
    return sorted(degs)


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def irreducible_sieve(Q, maxp=600):
    """(proved, detail): proved is True when the degree sets modulo primes leave only {0, n}."""
    n = len(Q) - 1
    if n == 1:
        return True, "degree 1"
    possible = set(range(n + 1))
    used = []
    for p in primes_upto(maxp):
        degs = ddf_degrees(Q, p)
        if degs is None:
            continue
        sums = {0}
        for x in degs:
            sums |= {s + x for s in sums}
        possible &= sums
        used.append((p, degs))
        if possible == {0, n}:
            return True, "primes %s" % ", ".join("%d%s" % (p, degs) for p, degs in used[-3:])
    return False, "factor degrees still possible: %s" % sorted(possible - {0, n})


def iv_mul_frac(a, b):
    c = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(c), max(c)


def irreducible_roots(Q, roots):
    """For monic (up to sign) real-rooted Q with all roots isolated: (True, None) if no monic integer factor of
    degree <= n/2 exists, (False, factor) otherwise. Exact rational interval arithmetic."""
    s = 1 if Q[-1] > 0 else -1
    M = [s * c for c in Q]
    assert M[-1] == 1
    n = len(M) - 1
    for k in range(1, n // 2 + 1):
        for sub in combinations(roots, k):
            poly = [(Fraction(1), Fraction(1))]            # constant term first
            for lo, hi in sub:
                new = [(Fraction(0), Fraction(0))] * (len(poly) + 1)
                for i, c in enumerate(poly):
                    up = new[i + 1]
                    new[i + 1] = (up[0] + c[0], up[1] + c[1])
                    rc = iv_mul_frac(c, (-hi, -lo))
                    new[i] = (new[i][0] + rc[0], new[i][1] + rc[1])
                poly = new
            cand = []
            for lo, hi in poly:
                a = -((-lo.numerator) // lo.denominator)    # ceil
                if a <= hi:
                    cand.append(a)
                else:
                    break
            if len(cand) == len(poly):
                rem = p_rem([Fraction(c) for c in M], [Fraction(c) for c in cand])
                if not rem:
                    return False, cand
    return True, None


# ----------------------------------------------------------------------------------------------------------------
# The certificate as data for the checks
# ----------------------------------------------------------------------------------------------------------------

class Cert:
    def __init__(self, raw):
        self.edges, self.dens, self.A, self.b = raw["edges"], raw["dens"], raw["A"], raw["b"]
        self.bound = raw["bound"]
        n = len(self.dens)
        self.n = n
        self.m = sum(d * (self.edges[i + 1] - self.edges[i]) for i, d in enumerate(self.dens))
        self.T = sum(d * (self.edges[i + 1] ** 2 - self.edges[i] ** 2) / 2 for i, d in enumerate(self.dens))
        self.support = [i for i, d in enumerate(self.dens) if d > 0]
        # jumps c_k = d_k - d_{k-1} at edge k (d_{-1} = d_n = 0)
        self.c = [(self.dens[k] if k < n else 0) - (self.dens[k - 1] if k > 0 else 0) for k in range(n + 1)]
        Dc = 1
        for d in self.dens:
            Dc = Dc * d.denominator // math.gcd(Dc, d.denominator)
        self.Dc = Dc
        self.k_lo, self.k_hi = (self.support[0], self.support[-1] + 1) if self.support else (0, 0)
        self.C = [int(self.c[k] * Dc) for k in range(self.k_lo, self.k_hi + 1)]

    def uniform(self):
        H = self.edges[1] - self.edges[0]
        return all(self.edges[i + 1] - self.edges[i] == H for i in range(self.n)), H


# ----------------------------------------------------------------------------------------------------------------
# Shared data for E4 and E5 (module globals, so forked workers inherit them copy-on-write)
# ----------------------------------------------------------------------------------------------------------------

SHD = {}


def prepare(cert, roots_per_Q):
    """Fill SHD: grid, integer jumps, the root list, Taylor tables."""
    ok, H = cert.uniform()
    assert ok, "the E5 engine needs equal cells"
    SHD.clear()
    ln_int.cache_clear()
    GH.cache_clear()
    SHD["H"] = H
    SHD["HN"], SHD["HD"] = H.numerator, H.denominator
    SHD["e0"] = cert.edges[0]
    SHD["LNH"] = ln_q(H.numerator, H.denominator)
    SHD["LNH_LP"] = ln_q_lp(H.numerator, H.denominator)
    Dc = 1
    for d in cert.dens:
        Dc = Dc * d.denominator // math.gcd(Dc, d.denominator)
    SHD["Dc"] = Dc
    first, last = cert.support[0], cert.support[-1]
    k_lo, k_hi = first, last + 1
    SHD["k_lo"], SHD["k_hi"] = k_lo, k_hi
    C = [int(cert.c[k] * Dc) for k in range(k_lo, k_hi + 1)]
    assert all(cert.c[k] == 0 for k in range(len(cert.c)) if not k_lo <= k <= k_hi)
    SHD["C"] = C
    SHD["Cabs"] = [abs(x) for x in C]
    SHD["Csum"] = sum(C)
    SHD["Cpos"] = sum(x for x in C if x > 0)
    SHD["Cneg"] = sum(x for x in C if x < 0)
    SHD["cfr"] = cert.c
    # roots: (b_j, lo, hi, j)
    roots = []
    for j, rs in enumerate(roots_per_Q):
        if cert.b[j] == 0:
            continue
        for lo, hi in rs:
            roots.append((cert.b[j], lo, hi, j))
    SHD["roots"] = roots
    SHD["sum_b_roots"] = sum(r[0] for r in roots)
    lcterm = (0, 0)
    for j, Q in enumerate(cert.A):
        if abs(Q[-1]) != 1 and cert.b[j] != 0:
            L = ln_int_lp(abs(Q[-1]))
            t = ivq(L, -cert.b[j].numerator, cert.b[j].denominator)
            lcterm = (lcterm[0] + t[0], lcterm[1] + t[1])
    SHD["lcterm_LP"] = lcterm
    # Taylor tables, indexed by idx = jmax - j, j = i - k, m = 2j + 1
    L = last - first
    jmin, jmax = -(L + 1), L
    SHD["jmax"] = jmax
    XM, XMa, XL, XML = [], [], [], []
    XV = {n: [] for n in range(2, NT + 1)}
    W = 0
    lncache = {}
    for idx in range(jmax - jmin + 1):
        jj = jmax - idx
        mm = 2 * jj + 1
        if abs(mm) == 1:
            XM.append(0), XMa.append(0), XL.append(0), XML.append(0)
            for n in XV:
                XV[n].append(0)
            continue
        am = abs(mm)
        if am not in lncache:
            lncache[am] = ln_int_lp(am)
        lo, hi = lncache[am]
        W = max(W, hi - lo)
        XM.append(mm), XMa.append(am), XL.append(lo), XML.append(mm * lo)
        for n in XV:
            pw = am ** (n - 1)
            if mm > 0 or (n - 1) % 2 == 0:
                XV[n].append((1 << LP) // pw)
            else:
                XV[n].append(-cdiv(1 << LP, pw))
    SHD.update(XM=XM, XMa=XMa, XL=XL, XML=XML, XV=XV, W=W)
    LNH_LP = SHD["LNH_LP"]
    SHD["KAPPA"] = (LNH_LP[0] - _LN2[1] - (1 << LP), LNH_LP[1] - _LN2[0] - (1 << LP))
    # global parts of the remainder bounds (root part; all far roots have |z| <= 2/3, rho <= 1/3)
    sb = SHD["sum_b_roots"]
    N = NT
    third = Fraction(1, 3)
    SHD["R0R"] = sb * third ** (N + 1) * Fraction(3, 2) / (N + 1)
    SHD["R1R"] = sb * third ** N / H
    SHD["R2R"] = sb * Fraction(4, 9) / (H * H) * third ** (N - 1) * (Fraction(3, 2) * N + Fraction(3, 4))
    for m_ in (3, 5):
        q = Fraction(1, m_) ** N / (1 - Fraction(1, m_))
        SHD["B0_%d" % m_] = q / (2 * N * (N + 1))
        SHD["B1_%d" % m_] = q / N
        SHD["B2_%d" % m_] = q


def ivadd(a, b):
    return a[0] + b[0], a[1] + b[1]


@lru_cache(maxsize=1 << 20)
def GH(p, q):
    """G(H p/q) (q > 0), interval at scale P."""
    if p == 0:
        return 0, 0
    L = ln_q(abs(p), q)
    LNH = SHD["LNH"]
    return ivq((L[0] + LNH[0] - ONE, L[1] + LNH[1] - ONE), SHD["HN"] * p, SHD["HD"] * q)


def lnH_rel(p, q):
    """ln(H |p|/q) at scale P."""
    L = ln_q(abs(p), q)
    LNH = SHD["LNH"]
    return L[0] + LNH[0], L[1] + LNH[1]


# ----------------------------------------------------------------------------------------------------------------
# Direct (termwise) evaluation: U at a rational point or root enclosure, E4, f at a rational point
# ----------------------------------------------------------------------------------------------------------------

def U_direct(cert, lo, hi):
    """Interval of U(r) for r in [lo, hi] (Fractions; lo == hi allowed). The enclosure must not contain an edge
    or an edge +-1 unless lo == hi."""
    C, k_lo = cert.C, cert.k_lo
    tot_lo = tot_hi = 0
    for off, Ck in enumerate(C):
        if Ck == 0:
            continue
        ek = cert.edges[k_lo + off]
        if lo == hi:
            u = lo - ek
            g = G_rat(u.numerator, u.denominator)
        else:
            u1, u2 = lo - ek, hi - ek
            assert u1 * u2 > 0 and (u1 - 1) * (u2 - 1) > 0 and (u1 + 1) * (u2 + 1) > 0, "refine the enclosure"
            g1, g2 = G_rat(u1.numerator, u1.denominator), G_rat(u2.numerator, u2.denominator)
            g = (min(g1[0], g2[0]), max(g1[1], g2[1]))
        if Ck > 0:
            tot_lo += Ck * g[0]
            tot_hi += Ck * g[1]
        else:
            tot_lo += Ck * g[1]
            tot_hi += Ck * g[0]
    return tot_lo // cert.Dc, cdiv(tot_hi, cert.Dc)


def safe_enclosure(Q, lo, hi, cert):
    """Refine an enclosure until no edge e_k and no e_k +- 1 lies in it (r is irrational)."""
    while True:
        bad = False
        for k in range(len(cert.edges)):
            ek = cert.edges[k]
            for x in (ek, ek - 1, ek + 1):
                if lo <= x <= hi:
                    bad = True
        if not bad:
            return lo, hi
        mid = (lo + hi) / 2
        if _sgn(ev_frac(Q, mid)) == _sgn(ev_frac(Q, lo)):
            lo = mid
        else:
            hi = mid


def P_direct(cert, x):
    """sum_j b_j ln|Q_j(x)| at a rational x (no root of any Q_j), interval."""
    tot = (0, 0)
    for j, Q in enumerate(cert.A):
        if cert.b[j] == 0:
            continue
        v = ev_frac(Q, x)
        assert v != 0
        L = ln_q(abs(v.numerator), v.denominator)
        tot = ivadd(tot, ivq(L, cert.b[j].numerator, cert.b[j].denominator))
    return tot


def f_direct(cert, x):
    u = U_direct(cert, x, x)
    p = P_direct(cert, x)
    return u[0] - p[1], u[1] - p[0]


# ----------------------------------------------------------------------------------------------------------------
# E5: the per-cell model and the bisection
# ----------------------------------------------------------------------------------------------------------------

def horner_point(coefs, num, k):
    """Interval value of sum coefs[n] tau^n at tau = num / 2^k."""
    lo, hi = coefs[-1]
    for c in reversed(coefs[:-1]):
        if num >= 0:
            lo, hi = (lo * num) >> k, -((-(hi * num)) >> k)
        else:
            lo, hi = (hi * num) >> k, -((-(lo * num)) >> k)
        lo += c[0]
        hi += c[1]
    return lo, hi


def horner_range(coefs, n1, n2, k):
    """Interval enclosure of sum coefs[n] tau^n over tau in [n1, n2] / 2^k."""
    lo, hi = coefs[-1]
    for c in reversed(coefs[:-1]):
        a, b, cc, d = lo * n1, lo * n2, hi * n1, hi * n2
        lo, hi = min(a, b, cc, d) >> k, -((-max(a, b, cc, d)) >> k)
        lo += c[0]
        hi += c[1]
    return lo, hi


def dyadic(fr):
    """(num, k) with fr = num / 2^k; fr must have a power-of-two denominator."""
    d = fr.denominator
    k = d.bit_length() - 1
    assert d == 1 << k
    return fr.numerator, k


class CellModel:
    """f on support cell i: x = e_i + rho H, rho in [0, 1], tau = rho - 1/2."""

    def __init__(self, i):
        H, HN, HD, Dc = SHD["H"], SHD["HN"], SHD["HD"], SHD["Dc"]
        self.i = i
        self.ei = SHD["e0"] + i * H
        x0 = self.ei + H / 2
        k_lo, C, jmax = SHD["k_lo"], SHD["C"], SHD["jmax"]
        a = jmax - i + k_lo
        b = a + len(C)
        Ci, Ci1 = C[i - k_lo], C[i + 1 - k_lo]
        Csum = SHD["Csum"] - Ci - Ci1
        Cpos = SHD["Cpos"] - max(Ci, 0) - max(Ci1, 0)
        Cneg = SHD["Cneg"] - min(Ci, 0) - min(Ci1, 0)
        W = SHD["W"]
        CM = sum(map(mul, C, SHD["XM"][a:b]))
        CMa = sum(map(mul, SHD["Cabs"], SHD["XMa"][a:b]))
        DL = sum(map(mul, C, SHD["XL"][a:b]))
        DML = sum(map(mul, C, SHD["XML"][a:b]))
        K = SHD["KAPPA"]
        co = [None] * (NT + 1)
        kc = (K[0] * CM, K[1] * CM) if CM >= 0 else (K[1] * CM, K[0] * CM)
        co[0] = ivq((DML - W * CMa + kc[0], DML + W * CMa + kc[1]), HN, 2 * HD * Dc)
        K1 = (K[0] + (1 << LP), K[1] + (1 << LP))
        kc1 = (K1[0] * Csum, K1[1] * Csum) if Csum >= 0 else (K1[1] * Csum, K1[0] * Csum)
        co[1] = ivq((DL + W * Cneg + kc1[0], DL + W * Cpos + kc1[1]), HN, HD * Dc)
        for n in range(2, NT + 1):
            DV = sum(map(mul, C, SHD["XV"][n][a:b]))
            co[n] = ivq((DV + Cneg, DV + Cpos), (-1) ** n * 2 ** (n - 1) * HN, n * (n - 1) * HD * Dc)
        co[0] = ivadd(co[0], SHD["lcterm_LP"])
        # root part
        self.near = []
        lim = 3 * H / 2
        for (bj, lo, hi, j) in SHD["roots"]:
            dl, dh = x0 - hi, x0 - lo
            if dl < lim and dh > -lim:
                self.near.append((bj, (lo - self.ei) / H, (hi - self.ei) / H, j))
                continue
            sg = 1 if dl > 0 else -1
            amin, amax = (dl, dh) if sg > 0 else (-dh, -dl)
            l1 = ln_q_lp(amin.numerator, amin.denominator)
            l2 = ln_q_lp(amax.numerator, amax.denominator)
            co[0] = ivadd(co[0], ivq((l1[0], l2[1]), -bj.numerator, bj.denominator))
            zl = (HN * amax.denominator << LP) // (HD * amax.numerator)
            zh = cdiv(HN * amin.denominator << LP, HD * amin.numerator)
            pl, ph = zl, zh
            for n in range(1, NT + 1):
                if n > 1:
                    pl = (pl * zl) >> LP
                    ph = cdiv(ph * zh, 1 << LP)
                if (-sg) ** n > 0:
                    iv = (pl, ph)
                else:
                    iv = (-ph, -pl)
                co[n] = ivadd(co[n], ivq(iv, bj.numerator, bj.denominator * n))
        self.coef = [(lo >> SH, cdiv(hi, 1 << SH)) for lo, hi in co]
        self.d1 = [(n * self.coef[n][0], n * self.coef[n][1]) for n in range(1, NT + 1)]
        self.d2 = [(n * (n - 1) * self.coef[n][0], n * (n - 1) * self.coef[n][1]) for n in range(2, NT + 1)]
        # remainder bounds (value, x-derivative, second x-derivative), as integers at scale P
        cfr = SHD["cfr"]
        c3 = sum(abs(cfr[k]) for k in (i - 1, i + 2) if 0 <= k < len(cfr))
        call = Fraction(SHD["Cpos"] - SHD["Cneg"] - abs(Ci) - abs(Ci1), Dc)
        c5 = call - c3
        R0 = H * (c3 * SHD["B0_3"] + c5 * SHD["B0_5"]) + SHD["R0R"]
        R1 = c3 * SHD["B1_3"] + c5 * SHD["B1_5"] + SHD["R1R"]
        R2 = 2 / H * (c3 * SHD["B2_3"] + c5 * SHD["B2_5"]) + SHD["R2R"]
        self.R0, self.R1, self.R2 = (qiv(R)[1] for R in (R0, R1, R2))
        self.Rfr = (R0, R1, R2)
        # own edges
        self.cL, self.cR = cfr[i], cfr[i + 1]
        self.own = []
        if self.cL != 0:
            self.own.append(("L", self.cL, self.cL > 0))      # c G(u), u = rho H > 0: convex iff c > 0
        if self.cR != 0:
            self.own.append(("R", self.cR, self.cR < 0))      # u = (rho - 1) H < 0: convex iff c < 0
        # floats for the heuristics
        self.cf = [ivfl(c) for c in self.coef]
        self.d1f = [n * self.cf[n] for n in range(1, NT + 1)]
        self.d2f = [n * (n - 1) * self.cf[n] for n in range(2, NT + 1)]
        self.Hf = float(H)
        self.cLf, self.cRf = float(self.cL), float(self.cR)
        self.nearf = [(float(bb), float((lo_ + hi_) / 2)) for bb, lo_, hi_, _ in self.near]
        self.cache = {}

    # ---- floats (heuristics only) ----
    def f_float(self, rho, sig=None):
        if sig is None:
            sig = 1.0 - rho
        tau = rho - 0.5
        acc = 0.0
        for c in reversed(self.cf):
            acc = acc * tau + c
        Hf = self.Hf
        if self.cLf and rho > 0:
            u = rho * Hf
            acc += self.cLf * u * (math.log(u) - 1)
        if self.cRf and sig > 0:
            u = sig * Hf
            acc -= self.cRf * u * (math.log(u) - 1)
        for bb, r in self.nearf:
            d = abs(rho - r)
            if d == 0:
                return float("inf")
            acc -= bb * (math.log(d * Hf))
        return acc

    def g_float(self, rho, sig):
        """df/drho = H f'(x), floats."""
        tau = rho - 0.5
        acc = 0.0
        for c in reversed(self.d1f):
            acc = acc * tau + c
        Hf = self.Hf
        if self.cLf:
            acc += Hf * self.cLf * math.log(rho * Hf) if rho > 0 else -math.copysign(float("inf"), self.cLf)
        if self.cRf:
            acc += Hf * self.cRf * math.log(sig * Hf) if sig > 0 else -math.copysign(float("inf"), self.cRf)
        for bb, r in self.nearf:
            acc -= bb / (rho - r)
        return acc

    def choose_t(self, j, s):
        """A point of [j/2^s, (j+1)/2^s] near the minimiser of f (floats; any point gives a valid bound)."""
        ra, rb = Fraction(j, 1 << s), Fraction(j + 1, 1 << s)
        right = 2 * (j + 1) > (1 << s)               # work with sigma = 1 - rho in the right half
        if not right:
            va, vb = float(ra), float(rb)

            def g(v):
                return self.g_float(v, 1.0 - v)
            singular = j == 0 and self.cLf != 0
        else:
            va, vb = float(1 - rb), float(1 - ra)

            def g(v):
                return -self.g_float(1.0 - v, v)      # increasing in v as rho decreases: sign flipped
            singular = j + 1 == (1 << s) and self.cRf != 0
        # in terms of v (distance from the nearer cell edge), we look for the root of the increasing function
        # h(v) = g(rho(v)) * orientation, oriented so that h increases with rho in the left half
        # left half: h = g increasing in v; right half: h(v) = -g, and g increasing in rho = decreasing in v,
        # so -g is increasing in v. In both cases: find v* with h(v*) = 0; f decreases where h < 0 in rho order.
        lo, hi = va, vb
        try:
            hlo = g(lo) if not (singular and lo == 0.0) else None
            hhi = g(hi)
        except (ValueError, ZeroDivisionError):
            return (ra + rb) / 2
        # left half: h < 0 means f decreasing in rho; minimiser to the right
        # right half: h(v) = -g: h < 0 means g > 0, f increasing in rho, decreasing in v: minimiser at larger v
        if hlo is None:
            hlo = -float("inf") if (self.cLf > 0 if not right else self.cRf < 0) else float("inf")
        if hlo >= 0:
            v = lo
        elif hhi <= 0:
            v = hi
        else:
            for _ in range(200):
                if lo < hi * 1e-3:
                    mid = math.sqrt(max(lo, hi * 1e-300) * hi)
                else:
                    mid = 0.5 * (lo + hi)
                if not lo < mid < hi:
                    break
                try:
                    hm = g(mid)
                except (ValueError, ZeroDivisionError):
                    break
                if hm < 0:
                    lo = mid
                else:
                    hi = mid
                if hi - lo <= 1e-15 * hi:
                    break
            v = hi
        t = Fraction(v) if not right else 1 - Fraction(v)
        return min(max(t, ra), rb)

    # ---- rigorous pieces ----
    def S_pt(self, num, k):
        lo, hi = horner_point(self.coef, num, k)
        return lo - self.R0, hi + self.R0

    def S_at(self, j, s):
        while s > 0 and j % 2 == 0:
            j //= 2
            s -= 1
        key = (j, s)
        v = self.cache.get(key)
        if v is None:
            v = self.S_pt(2 * j - (1 << s), s + 1)
            self.cache[key] = v
        return v

    def dS_pt(self, num, k):
        iv = ivq(horner_point(self.d1, num, k), SHD["HD"], SHD["HN"])
        return iv[0] - self.R1, iv[1] + self.R1

    def d2S_rng(self, n1, n2, k):
        HN, HD = SHD["HN"], SHD["HD"]
        iv = ivq(horner_range(self.d2, n1, n2, k), HD * HD, HN * HN)
        return iv[0] - self.R2, iv[1] + self.R2

    def own_val(self, side, c, p, q):
        """c G(u) at rho = p/q; u = rho H (L) or (rho - 1) H (R)."""
        g = GH(p, q) if side == "L" else GH(p - q, q)
        return ivq(g, c.numerator, c.denominator)

    def near_val_lo(self, bb, rlo, rhi, x):
        """Lower bound of -b ln|x - r| (x a Fraction rho, r in [rlo, rhi], x outside)."""
        d = max(abs(x - rlo), abs(x - rhi))
        L = lnH_rel(d.numerator, d.denominator)
        return ivq(L, -bb.numerator, bb.denominator)[0]

    def tangent(self, t, convex, roots):
        """Lower end of the tangent-group value at t, and the interval of its x-derivative."""
        Phi = Dlo = Dhi = 0
        p, q = t.numerator, t.denominator
        for side, c in convex:
            vp = p if side == "L" else p - q
            if vp == 0:
                return None
            val = self.own_val(side, c, p, q)
            L = lnH_rel(vp, q)
            der = ivq(L, c.numerator, c.denominator)
            Phi += val[0]
            Dlo += der[0]
            Dhi += der[1]
        H = SHD["H"]
        for bb, rlo, rhi, _ in roots:
            d1, d2 = t - rlo, t - rhi
            Phi += self.near_val_lo(bb, rlo, rhi, t)
            v1, v2 = -bb / (H * d1), -bb / (H * d2)
            vl, vh = (v1, v2) if v1 <= v2 else (v2, v1)
            Dlo += qiv(vl)[0]
            Dhi += qiv(vh)[1]
        return Phi, Dlo, Dhi

    def piece(self, j, s, goal, hint=None):
        """Certified lower bound of f on rho in [j/2^s, (j+1)/2^s]; returns (lb, t)."""
        H, HN, HD = SHD["H"], SHD["HN"], SHD["HD"]
        ra, rb = Fraction(j, 1 << s), Fraction(j + 1, 1 << s)
        SA, SB = self.S_at(j, s), self.S_at(j + 1, s)
        conc_a = conc_b = 0
        convex = []
        for side, c, cvx in self.own:
            if cvx:
                convex.append((side, c))
            else:
                conc_a += self.own_val(side, c, j, 1 << s)[0]
                conc_b += self.own_val(side, c, j + 1, 1 << s)[0]
        termw = 0
        troots = []
        for item in self.near:
            bb, rlo, rhi, _ = item
            if rlo <= rb and rhi >= ra:
                d = max(rb - rlo, rhi - ra)
                L = lnH_rel(d.numerator, d.denominator)
                termw += ivq(L, -bb.numerator, bb.denominator)[0]
            else:
                troots.append(item)
        n1, n2 = 2 * j - (1 << s), 2 * j + 2 - (1 << s)
        d2 = self.d2S_rng(n1, n2, s + 1)
        M = max(d2[1], 0)
        lossA = cdiv(M * HN * HN, 8 * HD * HD * (1 << (2 * s)))
        LA, LB = SA[0] + conc_a, SB[0] + conc_b
        best, tbest = None, None
        tries = [hint, None] if hint is not None else [None]
        if not convex and not troots:
            best = min(LA, LB) - lossA + termw
            if best >= goal:
                return best, None
        for t0 in tries:
            if t0 is None:
                t = self.choose_t(j, s)
            else:
                t = min(max(t0, ra), rb)
            tg = self.tangent(t, convex, troots)
            if tg is None:
                continue
            Phi, Dlo, Dhi = tg
            if convex or troots:
                lbA = min(LA + Phi + frq((ra - t) * H, Dhi, 1), LB + Phi + frq((rb - t) * H, Dlo, 1)) - lossA + termw
                if best is None or lbA > best:
                    best, tbest = lbA, t
                if best >= goal:
                    return best, tbest
            num, k = dyadic(t - Fraction(1, 2))
            St = self.S_pt(num, k)
            dSt = self.dS_pt(num, k)
            m2 = min(d2[0], 0)
            PhB, DlB, DhB = Phi + St[0], Dlo + dSt[0], Dhi + dSt[1]
            ca = frq(((ra - t) * H) ** 2 / 2, m2, 1)
            cb = frq(((rb - t) * H) ** 2 / 2, m2, 1)
            lbB = min(conc_a + PhB + frq((ra - t) * H, DhB, 1) + ca,
                      conc_b + PhB + frq((rb - t) * H, DlB, 1) + cb) + termw
            if best is None or lbB > best:
                best, tbest = lbB, t
            if best >= goal:
                return best, tbest
        if best is None:                 # no usable tangent point: report nothing certified, so the piece is split
            return -(1 << 400), None
        return best, tbest

    def f_upper(self, rho):
        """Rigorous upper bound of f at the point rho (Fraction with power-of-two denominator)."""
        num, k = dyadic(rho - Fraction(1, 2))
        v = self.S_pt(num, k)[1]
        p, q = rho.numerator, rho.denominator
        for side, c, _ in self.own:
            v += self.own_val(side, c, p, q)[1]
        for bb, rlo, rhi, _ in self.near:
            if rlo <= rho <= rhi:
                return None
            d = min(abs(rho - rlo), abs(rho - rhi))
            L = lnH_rel(d.numerator, d.denominator)
            v += ivq(L, -bb.numerator, bb.denominator)[1]
        return v

    def f_interval(self, rho):
        num, k = dyadic(rho - Fraction(1, 2))
        lo, hi = self.S_pt(num, k)
        p, q = rho.numerator, rho.denominator
        for side, c, _ in self.own:
            o = self.own_val(side, c, p, q)
            lo += o[0]
            hi += o[1]
        for bb, rlo, rhi, _ in self.near:
            dmin = min(abs(rho - rlo), abs(rho - rhi))
            dmax = max(abs(rho - rlo), abs(rho - rhi))
            L1 = lnH_rel(dmin.numerator, dmin.denominator)
            L2 = lnH_rel(dmax.numerator, dmax.denominator)
            lo += ivq(L2, -bb.numerator, bb.denominator)[0]
            hi += ivq(L1, -bb.numerator, bb.denominator)[1]
        return lo, hi

    def float_scan(self):
        """Float estimate of min f on the cell, and where (for choosing the goal and the upper bound)."""
        pts = [k / 64.0 for k in range(65)]
        for e in range(7, 61, 2):
            pts += [2.0 ** -e]
        best = (float("inf"), 0.5, 0.5)
        vals = []
        for r in pts:
            v = self.f_float(r, 1.0 - r)
            vals.append((v, r, 1.0 - r))
        for e in range(7, 61, 2):
            s_ = 2.0 ** -e
            vals.append((self.f_float(1.0 - s_, s_), 1.0 - s_, s_))
        vals.sort()
        for v, r, sg in vals[:4]:
            if v < best[0]:
                best = (v, r, sg)
            # golden-section refinement around r, in the scale of its distance to the nearer edge
            w = max(min(r, sg) / 2, 1.0 / 128)
            a, b = max(0.0, r - w), min(1.0, r + w)
            gr = (math.sqrt(5) - 1) / 2
            x1, x2 = b - gr * (b - a), a + gr * (b - a)
            f1, f2 = self.f_float(x1), self.f_float(x2)
            for _ in range(60):
                if f1 < f2:
                    b, x2, f2 = x2, x1, f1
                    x1 = b - gr * (b - a)
                    f1 = self.f_float(x1)
                else:
                    a, x1, f1 = x1, x2, f2
                    x2 = a + gr * (b - a)
                    f2 = self.f_float(x2)
            for xx, ff in ((x1, f1), (x2, f2)):
                if ff < best[0]:
                    best = (ff, xx, 1.0 - xx)
        return best


def bisect_cell(i, goal, smax=56, max_pieces=400000):
    """Certify f >= goal on support cell i by bisection. Returns a summary dict."""
    t0 = time.time()
    cm = CellModel(i)
    fmin, rmin, smin = cm.float_scan()
    stack = [(0, 0, None)]
    npieces = 0
    minlb, minat = None, None
    fails = []
    maxs = 0
    while stack:
        j, s, hint = stack.pop()
        lb, t = cm.piece(j, s, goal, hint)
        npieces += 1
        if lb >= goal:
            if minlb is None or lb < minlb:
                minlb, minat = lb, (j, s)
            maxs = max(maxs, s)
            continue
        if s >= smax or npieces > max_pieces:
            fails.append((j, s, lb))
            continue
        h = t if t is not None else hint
        stack.append((2 * j + 1, s + 1, h))
        stack.append((2 * j, s + 1, h))
    # rigorous upper bound of f at the float minimiser
    rho = Fraction(rmin) if rmin <= 0.5 else 1 - Fraction(smin)
    rho = min(max(rho, Fraction(0)), Fraction(1))
    ub = cm.f_upper(rho)
    return {"cell": i, "pieces": npieces, "minlb": minlb, "minat": minat, "fails": fails, "maxdepth": maxs,
            "fmin_float": fmin, "rho_min": float(rho), "ub": ub, "time": time.time() - t0}


def _bisect_job(args):
    i, goal = args
    return bisect_cell(i, goal)


def run_bisection(cells, goal, workers):
    if workers <= 1:
        return [bisect_cell(i, goal) for i in cells]
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    order = sorted(cells, key=lambda i: (i * 7919) % len(cells))
    with ctx.Pool(workers) as pool:
        out = list(pool.imap_unordered(_bisect_job, [(i, goal) for i in order], chunksize=4))
    return sorted(out, key=lambda r: r["cell"])


def _scan_job(i):
    cm = CellModel(i)
    v, r, sg = cm.float_scan()
    return i, v, r, sg


def run_scan(cells, workers):
    if workers <= 1:
        return [_scan_job(i) for i in cells]
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        return sorted(pool.imap_unordered(_scan_job, cells, chunksize=8))


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


def dec(fr, digits, up=False):
    """fr rounded down (or up) to `digits` decimals, as a string."""
    q, r = divmod(fr.numerator * 10 ** digits, fr.denominator)
    if up and r:
        q += 1
    s = str(abs(q)).rjust(digits + 1, "0")
    return ("-" if q < 0 else "") + s[:-digits] + "." + s[-digits:]


def fsci(v, scale=P):
    return "%.6e" % fx(v, scale)


def enc(iv, digits=10):
    """An interval at scale P as [lo, hi] in scientific notation, outward."""
    lo = Fraction(iv[0], ONE)
    hi = Fraction(iv[1], ONE)
    return "[%s, %s]" % (sci_floor(lo, digits), sci_ceil(hi, digits))


def sci_floor(fr, digits):
    if fr == 0:
        return "0"
    e = math.floor(math.log10(abs(float(fr))))
    m = fr / Fraction(10) ** e
    q = math.floor(m * 10 ** (digits - 1))
    return "%se%d" % (_mant(q, digits), e)


def sci_ceil(fr, digits):
    if fr == 0:
        return "0"
    e = math.floor(math.log10(abs(float(fr))))
    m = fr / Fraction(10) ** e
    q = math.ceil(m * 10 ** (digits - 1))
    return "%se%d" % (_mant(q, digits), e)


def _mant(q, digits):
    s = str(abs(q))
    out = s[0] + "." + s[1:] if len(s) > 1 else s
    return ("-" if q < 0 else "") + out


# ----------------------------------------------------------------------------------------------------------------
# The conditions
# ----------------------------------------------------------------------------------------------------------------

def check_E1(rep, cert):
    e, d = cert.edges, cert.dens
    ok = len(e) == len(d) + 1 and len(cert.A) == len(cert.b)
    rep.check("E1 shapes: %d edges, %d densities, %d polynomials, %d weights" % (len(e), len(d), len(cert.A),
                                                                                len(cert.b)), ok)
    rep.check("E1 edges strictly increasing", all(a < b for a, b in zip(e, e[1:])))
    rep.check("E1 every density >= 0", all(x >= 0 for x in d), "%d nonzero cells" % len(cert.support))
    rep.check("E1 every weight b_j >= 0", all(x >= 0 for x in cert.b),
              "%d positive" % sum(1 for x in cert.b if x > 0))
    rep.check("E1 mass m > 0", cert.m > 0, "m = %s = %s" % (cert.m, dec(cert.m, 18)))
    lo, hi = e[cert.support[0]], e[cert.support[-1] + 1]
    rep.check("E1 min S > 0", lo > 0, "S lies in [%s, %s]" % (lo, hi))
    gaps = sum(1 for a, b in zip(cert.support, cert.support[1:]) if b != a + 1)
    rep.info("support", "%d nonzero cells of width %s in [%s, %s] (%d cells in the hull, %d gaps)" % (
        len(cert.support), e[1] - e[0], lo, hi, cert.support[-1] - cert.support[0] + 1, gaps))
    return lo, hi


def check_E2(rep, cert):
    """Returns the isolated real roots of each Q (or None when a condition fails)."""
    roots, allok = [], True
    seen = {}
    for j, Q in enumerate(cert.A):
        name = pstr(Q)
        deg = len(Q) - 1
        ok = deg >= 1 and Q[-1] != 0
        g = 0
        for c in Q:
            g = math.gcd(g, c)
        prim = g == 1
        nreal, sqf = real_root_count(Q) if ok else (0, False)
        sieve, sdet = irreducible_sieve(Q) if ok else (False, "")
        if ok and sqf and nreal == deg:
            try:
                rs = isolate_roots(Q)
                if deg >= 2 and abs(Q[-1]) == 1:
                    rirr, fac = irreducible_roots(Q, rs)
                else:
                    rirr, fac = (True if deg == 1 else None), None
            except ValueError as ex:
                rs, rirr, fac = None, False, str(ex)
        else:
            rs, rirr, fac = None, None, None
        key = tuple(Q)
        neg = tuple(-c for c in Q)
        assoc = key in seen or neg in seen
        seen[key] = j
        irr = (sieve or rirr is True) and rirr is not False
        good = ok and prim and irr and sqf and nreal == deg and not assoc
        allok &= good
        det = "deg %d, gcd %d, %s, roots-method %s, %d distinct real roots%s, b = %s" % (
            deg, g, sdet if sieve else "sieve: " + sdet, {True: "irreducible", False: "REDUCIBLE %s" % (fac,),
                                                          None: "n/a"}[rirr], nreal, ", ASSOCIATE" if assoc else "",
            cert.b[j])
        rep.check("E2 Q%d = %s: primitive, irreducible, all roots real, not associate to another" % (j, name), good,
                  det)
        roots.append(rs)
    return roots if allok else None


def check_E3(rep, cert):
    s = sum(b * (len(Q) - 1) for b, Q in zip(cert.b, cert.A))
    ok = s < cert.m
    rep.check("E3 sum_j b_j deg Q_j < m", ok, "sum = %s, sum/m = %s" % (dec(s, 12), dec(s / cert.m, 8)))
    return s / cert.m


def E4_values(cert, roots):
    """[(j, interval of int log|Q_j| dmu)], with the safe enclosures."""
    out = []
    for j, Q in enumerate(cert.A):
        tot = (0, 0)
        if abs(Q[-1]) != 1:
            L = ln_int(abs(Q[-1]))
            tot = ivq(L, cert.m.numerator, cert.m.denominator)
        for lo, hi in roots[j]:
            if lo != hi:
                lo, hi = safe_enclosure(Q, lo, hi, cert)
            tot = ivadd(tot, U_direct(cert, lo, hi))
        out.append((j, tot))
    return out


def check_E4(rep, cert, roots, stated_min=None):
    vals = E4_values(cert, roots)
    worst = min(vals, key=lambda v: v[1][0])
    for j, iv in vals:
        rep.info("E4 int log|Q%d| dmu, Q%d = %s" % (j, j, pstr(cert.A[j])), enc(iv, 12))
    ok = all(iv[0] > 0 for _, iv in vals)
    rep.check("E4 int log|Q_j| dmu > 0 for every j (rigorous lower ends)", ok,
              "smallest: Q%d = %s, in %s" % (worst[0], pstr(cert.A[worst[0]]), enc(worst[1], 12)))
    return vals


def check_E6(rep, cert, printed):
    b = cert.T / cert.m
    rep.info("T = int x dmu", "%s" % cert.T)
    rep.check("E6 T/m equals the certificate's bound exactly", b == cert.bound,
              "T/m = %s = %s..." % (b, dec(b, 22)))
    for where, s in printed:
        p = Fraction(s)
        rep.check("E6 printed upper bound %s (%s) >= T/m" % (s, where), p >= b,
                  "it is T/m rounded up at %d decimals" % (len(s.split(".")[1])) if p == Fraction(dec(b, len(s.split(".")[1]), up=True)) else
                  "difference %.3e" % float(p - b))
    return b


def energy(cert):
    """I(mu) = int int log|x - y| dmu dmu for equal cells: sum_j R(j) K(j), R the autocorrelation of the densities,
    K(j) = H^2 [((j+1)^2 ln|j+1| - 2 j^2 ln|j| + (j-1)^2 ln|j-1|)/2 + ln H - 3/2]."""
    ok, H = cert.uniform()
    s0, s1 = cert.support[0], cert.support[-1]
    D = [int(x * cert.Dc) for x in cert.dens[s0:s1 + 1]]
    n = len(D)
    R = [sum(map(mul, D[:n - j], D[j:])) for j in range(n)]
    LNH = ln_q(H.numerator, H.denominator)

    def xlx(k):     # k^2 ln k, k >= 0
        if k == 0:
            return 0, 0
        L = ln_int(k)
        return k * k * L[0], k * k * L[1]

    lo = hi = 0
    for j in range(n):
        a, b, c = xlx(j + 1), xlx(j), xlx(abs(j - 1))
        klo = (a[0] - 2 * b[1] + c[0]) // 2 + LNH[0] - 3 * ONE // 2
        khi = cdiv(a[1] - 2 * b[0] + c[1], 2) + LNH[1] - 3 * ONE // 2 + 1
        w = 1 if j == 0 else 2
        lo += w * R[j] * klo
        hi += w * R[j] * khi
    f = H * H / (cert.Dc * cert.Dc)
    return ivq((lo, hi), f.numerator, f.denominator)


# ----------------------------------------------------------------------------------------------------------------
# E5 driver
# ----------------------------------------------------------------------------------------------------------------

def check_E5(rep, cert, workers, claim, rel=Fraction(1, 10 ** 5)):
    """Certify f >= g on every support cell, with g just below the float estimate of min_S f (and g >= the claimed
    margin). Then tighten the cell of the minimum. Returns (LB, UB, x, results)."""
    t0 = time.time()
    scan = run_scan(cert.support, workers)
    gi, gv, gr, gs = min(scan, key=lambda r: r[1])
    rep.info("E5 float scan (heuristic only)", "smallest f found %.6e in cell %d at x = %.7f; %.1f s" % (
        gv, gi, float(cert.edges[gi]) + (gr if gr <= 0.5 else 1 - gs) * float(SHD["H"]), time.time() - t0))
    g = qiv(Fraction(gv) * (1 - rel))[0] if gv > 0 else 0
    gclaim = qiv(claim)[1]
    if g < gclaim:
        g = gclaim
    t1 = time.time()
    res = run_bisection(cert.support, g, workers)
    failed = [r for r in res if r["fails"]]
    if failed:      # retry those cells at the weaker goal (the claimed margin, then 0)
        for goal in (gclaim, 0):
            redo = run_bisection([r["cell"] for r in failed], goal, workers)
            byc = {r["cell"]: r for r in redo}
            res = [byc.get(r["cell"], r) for r in res]
            failed = [r for r in res if r["fails"]]
            if not failed:
                break
    npieces = sum(r["pieces"] for r in res)
    LB = min(r["minlb"] for r in res if r["minlb"] is not None) if any(r["minlb"] is not None for r in res) else None
    rep.check("E5 f >= %s on every support cell (%d cells, %d pieces, max depth %d)" % (
        sci_floor(Fraction(g, ONE), 6), len(res), npieces, max(r["maxdepth"] for r in res)), not failed,
        "%.1f s; %d cells failed" % (time.time() - t1, len(failed)))
    # tighten the cell of the float minimum, and take the rigorous upper bound there
    tight = bisect_cell(gi, qiv(Fraction(gv) * (1 - Fraction(1, 10 ** 5)))[0] if gv > 0 else 0)
    UB = min(r["ub"] for r in res if r["ub"] is not None)
    ubcell = min((r for r in res if r["ub"] is not None), key=lambda r: r["ub"])
    LBt = min([tight["minlb"]] + [r["minlb"] for r in res if r["cell"] != gi and r["minlb"] is not None])
    x = cert.edges[ubcell["cell"]] + Fraction(ubcell["rho_min"]) * SHD["H"]
    rep.info("E5 smallest certified piece bound (all cells, goal %s)" % sci_floor(Fraction(g, ONE), 6),
             sci_floor(Fraction(LB, ONE), 8) if LB is not None else "none")
    rep.check("E5 min over S of f enclosed", not tight["fails"] and LBt <= UB,
              "min_S f in [%s, %s]; upper end: rigorous f at x = %.9f (cell %d); lower end: cell %d re-run with goal "
              "(1 - 1e-5) x the float minimum, %d pieces" % (
                  sci_floor(Fraction(LBt, ONE), 7), sci_ceil(Fraction(UB, ONE), 7), float(x), ubcell["cell"], gi,
                  tight["pieces"]))
    return LBt, UB, x, res, failed


def gap_points(cert):
    """Midpoints of the zero-density cells inside the hull, and two exterior points."""
    pts = []
    s = cert.support
    for a, b in zip(s, s[1:]):
        if b != a + 1:
            pts.append((cert.edges[a + 1] + cert.edges[b]) / 2)
    lo, hi = cert.edges[s[0]], cert.edges[s[-1] + 1]
    pts += [lo / 2, hi + Fraction(1, 5), Fraction(10)]
    return pts


# ----------------------------------------------------------------------------------------------------------------
# Red controls
# ----------------------------------------------------------------------------------------------------------------

def modified(cert, dens=None, b=None, A=None, bound=None):
    raw = {"edges": cert.edges, "dens": dens if dens is not None else cert.dens, "A": A if A is not None else cert.A,
           "b": b if b is not None else cert.b, "bound": bound if bound is not None else cert.bound}
    return Cert(raw)


def refute_E5(cert2, candidates):
    """Rigorous f of the modified certificate at candidate points; the most negative upper bound found."""
    best = None
    for x in candidates:
        lo, hi = f_direct(cert2, x)
        if best is None or hi < best[1]:
            best = (x, hi)
    return best


def controls(rep, cert, roots, scan_pts, res_e5, LB, UB):
    """Each must come out wrong."""
    # 1. one density lowered by 1% (the largest), bound recomputed
    i = max(cert.support, key=lambda k: cert.dens[k])
    d2 = list(cert.dens)
    d2[i] = d2[i] * Fraction(99, 100)
    c2 = modified(cert, dens=d2, bound=None)
    vals = E4_values(c2, roots)
    neg4 = [j for j, iv in vals if iv[1] < 0]
    best = refute_E5(c2, scan_pts)
    rep.check("control: density of cell %d (the largest, %.4f) lowered by 1%%: E4 or E5 refuted" % (
        i, float(cert.dens[i])), bool(neg4) or best[1] < 0,
        "E4 negative for %s; f < %s at x = %.9f (rigorous upper bound)" % (
            ["Q%d" % j for j in neg4] or "none", sci_ceil(Fraction(best[1], ONE), 4), float(best[0])))
    # 2. the largest weight raised by 0.1%
    jmax = max(range(len(cert.b)), key=lambda j: cert.b[j])
    b2 = list(cert.b)
    b2[jmax] = b2[jmax] * Fraction(1001, 1000)
    c3 = modified(cert, b=b2)
    best = refute_E5(c3, scan_pts)
    rep.check("control: largest weight (Q%d = %s, b = %.6f) raised by 0.1%%: E5 refuted" % (
        jmax, pstr(cert.A[jmax]), float(cert.b[jmax])), best[1] < 0,
        "f < %s at x = %.9f (rigorous upper bound)" % (sci_ceil(Fraction(best[1], ONE), 4), float(best[0])))
    # 3. E2 refusals: reducible, not totally real, not primitive, associate
    for label, Q in (("x^2 - 3x + 2 = (x - 1)(x - 2)", [2, -3, 1]), ("x^2 + 1 (not totally real)", [1, 0, 1]),
                     ("2x^2 - 6x + 2 (not primitive)", [2, -6, 2]),
                     ("(x - 1)(x^2 - 3x + 1) = x^3 - 4x^2 + 4x - 1", [-1, 4, -4, 1]),
                     ("x^4 - 10x^2 + 1 (irreducible, reducible mod every prime)", [1, 0, -10, 0, 1])):
        A2 = list(cert.A) + [Q]
        c4 = modified(cert, A=A2, b=list(cert.b) + [Fraction(1, 10 ** 6)])
        quiet = Report.__new__(Report)
        quiet.fails = []
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            r = check_E2(quiet, c4)
        last = [ln for ln in buf.getvalue().splitlines() if ln.strip()][-1]
        verdict = last.split(": ", 1)[1]
        why = verdict[verdict.index("(") + 1:-1] if "(" in verdict else verdict
        refused = r is None
        expect = "irreducible, reducible mod every prime" not in label
        if expect:
            rep.check("control: E2 refuses %s" % label, refused, "E2 on it: " + why[:160])
        else:
            rep.check("control: E2 accepts %s only through the root method (the sieve cannot prove it)" % label,
                      not refused and "sieve:" in last, "E2 on it: " + why[:160])
    # 4. E3: weights scaled by 1.3
    s = sum(bb * Fraction(13, 10) * (len(Q) - 1) for bb, Q in zip(cert.b, cert.A))
    rep.check("control: weights scaled by 1.3 violate E3", not s < cert.m, "sum/m = %s" % dec(s / cert.m, 6))
    # 5. E6: the bound lowered by 1e-15
    b = cert.T / cert.m
    low = b - Fraction(1, 10 ** 15)
    rep.check("control: a bound 1e-15 below T/m is refused by E6", low != b and low < b)
    # 6. the bisection cannot certify a goal above a rigorous upper bound of f
    gi = min(res_e5, key=lambda r: r["ub"] if r["ub"] is not None else 1 << 400)["cell"]
    r = bisect_cell(gi, UB + qiv(Fraction(1, 10 ** 12))[1], smax=24, max_pieces=20000)
    rep.check("control: goal min_S f + 1e-12 is not certified in cell %d (bisection to depth 24)" % gi,
              bool(r["fails"]), "%d pieces, %d uncertified" % (r["pieces"], len(r["fails"])))


# ----------------------------------------------------------------------------------------------------------------
# Self-tests (run with every row)
# ----------------------------------------------------------------------------------------------------------------

def selftests(rep, cert, rng):
    getcontext().prec = 70
    bad = 0
    ns = [rng.randrange(1, 10 ** rng.randrange(1, 80)) for _ in range(600)] + list(range(1, 300))
    for n in ns:
        lo, hi = ln_int(n)
        ref = Fraction(Decimal(n).ln())
        tol = Fraction(1, 10 ** 60)
        if not Fraction(lo, ONE) - tol <= ref <= Fraction(hi, ONE) + tol:
            bad += 1
    rep.check("self-test: our ln encloses the correctly rounded decimal ln (70 digits) at %d integers" % len(ns),
              bad == 0)
    # G and U against numerical quadrature (floats, tolerance 1e-9): U at a point equals the sum over cells of
    # the integral of log|x - y|, here by 20-point Gauss-Legendre on 64 subcells (the cell holding x split at x)
    nodes, weights = _gauss_legendre(20)
    worst = 0.0
    for _ in range(3):
        i = rng.choice(cert.support)
        x = cert.edges[i] + Fraction(rng.randrange(1, 1000), 1000) * (cert.edges[i + 1] - cert.edges[i])
        xf = float(x)
        tot = 0.0
        for k in cert.support:
            a, b = float(cert.edges[k]), float(cert.edges[k + 1])
            parts = [(a, xf), (xf, b)] if a < xf < b else [(a, b)]
            s = sum(_log_integral(xf, p, q, nodes, weights) for p, q in parts)
            tot += float(cert.dens[k]) * s
        iv = U_direct(cert, x, x)
        worst = max(worst, abs(ivfl(iv) - tot))
    rep.check("self-test: closed-form U agrees with Gauss-Legendre quadrature at 3 random points", worst < 1e-9,
              "largest difference %.2e" % worst)


def _log_integral(x, p, q, nodes, weights):
    """int_p^q log|x - y| dy in floats, x outside (p, q): Gauss-Legendre in the distance from the end nearer to x,
    on a mesh graded geometrically towards that end (self-test only). |x - y| = d0 + (distance from the near end)."""
    near, far = (p, q) if abs(x - p) <= abs(x - q) else (q, p)
    L = abs(far - near)
    d0 = abs(x - near)
    cuts, w = [L], L
    while w > 1e-30 * L and w > 0.25 * d0:
        w /= 2
        cuts.append(w)
    pts = sorted(set(cuts + [0.0]))
    s = 0.0
    if d0 == 0.0:                       # [0, eps] analytically: int_0^eps log u du = eps (log eps - 1)
        eps = pts[1]
        s += eps * (math.log(eps) - 1)
        pts = pts[1:]
    for u, v in zip(pts, pts[1:]):
        h = (v - u) / 2
        for z, wt in zip(nodes, weights):
            s += wt * h * math.log(d0 + u + h * (z + 1))
    return s


def _gauss_legendre(n):
    """Nodes and weights on [-1, 1] by Newton on the Legendre recurrence (floats; used only in a self-test)."""
    xs, ws = [], []
    for k in range(1, n + 1):
        x = math.cos(math.pi * (k - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for m in range(2, n + 1):
                p0, p1 = p1, ((2 * m - 1) * x * p1 - (m - 1) * p0) / m
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-16:
                break
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def model_vs_direct(rep, cert, rng, count=24):
    worst = 0.0
    ok = True
    for _ in range(count):
        i = rng.choice(cert.support)
        cm = CellModel(i)
        rho = Fraction(rng.randrange(1, 1 << 30), 1 << 30)
        a = cm.f_interval(rho)
        b = f_direct(cert, cm.ei + rho * SHD["H"])
        ok &= a[0] <= b[1] and b[0] <= a[1]
        worst = max(worst, abs(ivfl(a) - ivfl(b)))
    rep.check("self-test: the cell Taylor model and direct termwise sums agree (enclosures overlap) at %d random "
              "points" % count, ok, "largest midpoint difference %.2e" % worst)


def piece_soundness(rep, cert, rng, count=240):
    """For random pieces (random support cells, depths 0-44, many touching a cell edge or holding a root), the lower
    bound of `piece` never exceeds rigorous upper bounds of f at about 120 points of the piece, crowded towards both
    ends (where the edge cusps are)."""
    bad, n, tight = 0, 0, []
    models = {}
    cells = list(cert.support)
    for trial in range(count):
        i = rng.choice(cells)
        cm = models.get(i) or models.setdefault(i, CellModel(i))
        s = rng.choice([0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 44])
        mode = trial % 4
        if mode == 0:
            j = 0
        elif mode == 1:
            j = (1 << s) - 1
        elif mode == 2 and cm.near:
            r = rng.choice(cm.near)[1]
            j = min(max(int(r * (1 << s)), 0), (1 << s) - 1)
        else:
            j = rng.randrange(1 << s)
        lb, t = cm.piece(j, s, 1 << 400)
        ra, rb = Fraction(j, 1 << s), Fraction(j + 1, 1 << s)
        w = rb - ra
        pts = [ra, rb] + [ra + w * Fraction(k, 32) for k in range(1, 32)]
        pts += [ra + w / 2 ** k for k in range(1, 46)] + [rb - w / 2 ** k for k in range(1, 46)]
        if t is not None:
            pts.append(t)
        best = None
        for p in pts:
            ub = cm.f_upper(p)
            if ub is None:
                continue
            n += 1
            if ub < lb:
                bad += 1
            best = ub if best is None else min(best, ub)
        if best is not None:
            tight.append(fx(best - lb))
    tight.sort()
    rep.check("self-test: piece lower bounds never exceed rigorous point values of f (%d random pieces, depths 0-44, "
              "%d points)" % (count, n), bad == 0,
              "%d violations; sampled minimum minus bound: median %.1e, max %.1e" % (
                  bad, tight[len(tight) // 2], tight[-1]))


# ================================================================================================================
# PR #298: a piecewise-linear density on non-uniform rational nodes, auxiliaries with nonreal roots
# ================================================================================================================
#
# The density rho interpolates the node heights h_k linearly between the nodes e_k and is zero outside them (the
# end heights are 0, so rho is continuous). Then rho'' = sum_k j_k delta(e_k), j_k the slope jump at e_k, and two
# integrations by parts give, for every z in C,
#     U(z) = int log|z - y| rho(y) dy = Re sum_k j_k F(z - e_k),     F(w) = w^2 log(w)/2 - 3 w^2/4  (F'' = log),
# with the principal logarithm (for Im z > 0 every z - e_k lies in the upper half-plane, so the branch is
# continuous in y); on the real line F(u) = G2(u) = u^2 log|u|/2 - 3u^2/4, G2(0) = 0. sum j_k = sum j_k e_k = 0 and
# sum j_k e_k^2/2 = m are checked exactly. Every logarithm and arctangent below is our own, outward-rounded.

ROOT_BITS = 300                 # root certification works at 2^-(Mb + ROOT_BITS), Mb bounding Horner's partials
KC = 160                        # nonreal root centres are rounded to the grid 2^-KC (the rounding joins the radius)


# ---- roots: our own Aberth iteration, then Gershgorin's theorem on the Weierstrass matrix (copied from c82) ----

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


def certified_discs(Q):
    """Every root of the integer polynomial Q (constant term first) in its own disc: pairwise disjoint discs
    (cx, cy, rad) as Fractions, each holding exactly one root. Returns (discs, working precision in bits)."""
    zs, prec, _ = approx_roots(Q, ROOT_BITS)
    if zs is None:
        raise ValueError("Aberth iteration did not converge")
    discs = gershgorin(Q, zs, prec)
    return [(Fraction(x, 1 << prec), Fraction(y, 1 << prec), Fraction(r, 1 << prec)) for x, y, r in discs], prec


def split_discs(discs, nreal):
    """Discs above, below and meeting the real axis. When there are nreal (the exact Sturm count) discs meeting the
    axis and the rest split evenly, the discs off the axis hold exactly the nonreal roots (each disc holds one
    root, and an off-axis disc cannot hold a real one)."""
    up = [d for d in discs if d[1] - d[2] > 0]
    down = [d for d in discs if d[1] + d[2] < 0]
    mid = [d for d in discs if not (d[1] - d[2] > 0 or d[1] + d[2] < 0)]
    return len(mid) == nreal and len(up) == len(down), up, down, mid


def round_disc(cx, cy, rad):
    """Centre rounded to the grid 2^-KC (integers CX, CY), radius enlarged by the rounding."""
    CX = (2 * cx.numerator * (1 << KC) + cx.denominator) // (2 * cx.denominator)
    CY = (2 * cy.numerator * (1 << KC) + cy.denominator) // (2 * cy.denominator)
    return CX, CY, rad + Fraction(1, 1 << KC)


# ---- arctangent and pi (rigorous) ----

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


_PI_FR = pi_machin(LP + 40)
PI_LP = ((_PI_FR[0].numerator << LP) // _PI_FR[0].denominator, cdiv(_PI_FR[1].numerator << LP, _PI_FR[1].denominator))


def _atan_euler(p, q):
    """atan(p/q) for integers 0 <= p <= q, q > 0, at scale 2^-LP, by Euler's series
    atan(x) = (x/(1 + x^2)) sum_n c_n y^n, y = x^2/(1 + x^2) <= 1/2, c_0 = 1, c_n = c_{n-1} 2n/(2n + 1):
    positive terms whose ratio is at most y <= 1/2, so the tail after a term is at most that term."""
    if p == 0:
        return 0, 0
    S = LP + 20
    B = p * p + q * q
    tl = ((p * q) << S) // B
    th = cdiv((p * q) << S, B)
    lo = hi = 0
    n = 0
    while True:
        lo += tl
        hi += th
        if th <= 1:
            hi += th                          # the tail
            break
        n += 1
        tl = (tl * 2 * n * p * p) // ((2 * n + 1) * B)
        th = cdiv(th * 2 * n * p * p, (2 * n + 1) * B)
    return lo >> 20, cdiv(hi, 1 << 20)


_ATAN_TAB = []


def atan_tab():
    if not _ATAN_TAB:
        _ATAN_TAB.extend(_atan_euler(j, 256) for j in range(257))      # atan(j/256)
    return _ATAN_TAB


def _atan_small(dn, dd):
    """atan(dn/dd) for |dn/dd| <= 1/512 (dd > 0), at scale 2^-LP: the alternating series with decreasing terms, so
    the error of a partial sum is at most the next term."""
    if dn == 0:
        return 0, 0
    neg = dn < 0
    dn = -dn if neg else dn
    S = LP + 16
    pl, ph = (dn << S) // dd, cdiv(dn << S, dd)
    d2l, d2h = ((dn * dn) << S) // (dd * dd), cdiv((dn * dn) << S, dd * dd)
    lo = hi = 0
    k = 0
    while True:
        tl, th = pl // (2 * k + 1), cdiv(ph, 2 * k + 1)
        if k % 2 == 0:
            lo += tl
            hi += th
        else:
            lo -= th
            hi -= tl
        k += 1
        pl, ph = (pl * d2l) >> S, cdiv(ph * d2h, 1 << S)
        if ph <= 1:                           # the next term, which bounds the error, is below one unit
            lo -= 2
            hi += 2
            break
    lo, hi = lo >> 16, cdiv(hi, 1 << 16)
    return (-hi, -lo) if neg else (lo, hi)


def atan_q(p, q):
    """atan(p/q) for integers p >= 0, q > 0, as an interval at scale 2^-LP."""
    if p == 0:
        return 0, 0
    if p > q:                                 # atan x = pi/2 - atan(1/x)
        lo, hi = atan_q(q, p)
        return (PI_LP[0] >> 1) - hi, cdiv(PI_LP[1], 2) - lo
    j = (512 * p + q) // (2 * q)              # nearest integer to 256 p/q, in [0, 256]
    dn, dd = 256 * p - j * q, 256 * q + j * p  # atan x = atan(j/256) + atan((x - j/256)/(1 + x j/256))
    T = atan_tab()[j]
    a = _atan_small(dn, dd)
    return T[0] + a[0], T[1] + a[1]


def atan2_lp(Y, n):
    """arg(n + iY) in (0, pi) for integers Y > 0, n, at scale 2^-LP."""
    if n > 0:
        return atan_q(Y, n)
    if n == 0:
        return PI_LP[0] >> 1, cdiv(PI_LP[1], 2)
    lo, hi = atan_q(Y, -n)
    return PI_LP[0] - hi, PI_LP[1] - lo


# ---- the certificate as data ----

class PLCert:
    """Nodes e_k, heights h_k (linear between nodes, zero outside), polynomials A (constant term first), weights b."""

    def __init__(self, e, h, A, b, bound):
        self.e, self.h, self.A, self.b, self.bound = list(e), list(h), [list(Q) for Q in A], list(b), bound
        N = len(self.e)
        e, h = self.e, self.h
        self.m = sum((h[k] + h[k + 1]) / 2 * (e[k + 1] - e[k]) for k in range(N - 1))
        self.T = sum((e[k + 1] - e[k]) * (h[k] * (2 * e[k] + e[k + 1]) + h[k + 1] * (e[k] + 2 * e[k + 1])) / 6
                     for k in range(N - 1))
        s = [(h[k + 1] - h[k]) / (e[k + 1] - e[k]) for k in range(N - 1)]
        self.jump = [s[0]] + [s[k] - s[k - 1] for k in range(1, N - 1)] + [-s[N - 2]]
        self.support = [k for k in range(N - 1) if h[k] > 0 or h[k + 1] > 0]


PL = {}            # module globals for E4 and E5 (forked workers inherit them)
_LNT_CACHE = {}


def pl_prepare(cert, rroots, croots):
    """rroots[j]: the real roots of Q_j as enclosures (lo, hi); croots[j]: the roots of Q_j in the upper half-plane
    as (CX, CY, rad), centre (CX + i CY)/2^KC, radius rad."""
    PL.clear()
    e, jump = cert.e, cert.jump
    Le = 1
    for x in e:
        Le = Le * x.denominator // math.gcd(Le, x.denominator)
    D = 2 * Le                                  # D e_k is even, so cell midpoints and half-widths are k/D
    DJ = 1
    for x in jump:
        DJ = DJ * x.denominator // math.gcd(DJ, x.denominator)
    EI = []
    for x in e:
        v = x * D
        assert v.denominator == 1
        EI.append(v.numerator)
    nz = [k for k in range(len(e)) if jump[k] != 0]
    NE = [EI[k] for k in nz]
    NJ = []
    for k in nz:
        v = jump[k] * DJ
        assert v.denominator == 1
        NJ.append(v.numerator)
    s0, s1 = cert.support[0], cert.support[-1]
    x0lo, x0hi = (EI[s0] + EI[s0 + 1]) // 2, (EI[s1] + EI[s1 + 1]) // 2
    nmax = max(max(NE) - x0lo, x0hi - min(NE))
    key = (D, nmax)
    if key not in _LNT_CACHE:
        _LNT_CACHE.clear()
        LNT, W = [0] * (nmax + 1), 0
        for n in range(1, nmax + 1):
            lo, hi = ln_int_lp(n)
            LNT[n] = lo
            W = max(W, hi - lo)
        _LNT_CACHE[key] = (LNT, W)
    LNT, W = _LNT_CACHE[key]
    rr, cr = [], []
    for j in range(len(cert.A)):
        if cert.b[j] == 0:
            continue
        for lo, hi in rroots[j]:
            rr.append((cert.b[j], (lo + hi) / 2, (hi - lo) / 2, j))
        for CX, CY, rad in croots[j]:
            cr.append((2 * cert.b[j], CX, CY, rad, j))
    lc = (0, 0)
    for j, Q in enumerate(cert.A):
        if abs(Q[-1]) != 1 and cert.b[j] != 0:
            t = ivq(ln_int_lp(abs(Q[-1])), -cert.b[j].numerator, cert.b[j].denominator)
            lc = (lc[0] + t[0], lc[1] + t[1])
    PL.update(D=D, DJ=DJ, EI=EI, NE=NE, NJ=NJ, LNT=LNT, LNTW=W, LND=ln_int_lp(D), rr=rr, cr=cr, lc=lc, m=cert.m,
              A=cert.A, b=cert.b, e=e, jabs=sum(abs(x) for x in jump), rroots=rroots, croots=croots, cert=cert)


def lnM_lp(k):
    """ln(D 2^k) at scale 2^-LP."""
    L = PL["LND"]
    return L[0] + k * _LN2[0], L[1] + k * _LN2[1]


# ---- the potential and f at single points (direct sums over every node) ----

def U_real_lp(p, q):
    """U(p/q) (q > 0) at scale 2^-LP: sum_k j_k G2(u_k), u_k = n_k/M with n_k = p D - E_k q, M = D q, so
    4 DJ M^2 U = 2 sum J_k n_k^2 ln|n_k| - (2 ln M + 3) sum J_k n_k^2 (J_k = DJ j_k, E_k = D e_k integers)."""
    D, DJ = PL["D"], PL["DJ"]
    alo = ahi = B = 0
    pD = p * D
    for Ek, Jk in zip(PL["NE"], PL["NJ"]):
        n = pD - Ek * q
        if n == 0:
            continue
        c = Jk * n * n
        L = ln_int_lp(n if n > 0 else -n)
        if c > 0:
            alo += c * L[0]
            ahi += c * L[1]
        else:
            alo += c * L[1]
            ahi += c * L[0]
        B += c
    M = D * q
    lm = ln_int_lp(M)
    xb = (lm[0] * B, lm[1] * B) if B >= 0 else (lm[1] * B, lm[0] * B)
    nlo = 2 * (alo - xb[1]) - 3 * B * (1 << LP)
    nhi = 2 * (ahi - xb[0]) - 3 * B * (1 << LP)
    den = 4 * DJ * M * M
    return nlo // den, cdiv(nhi, den)


def U_complex_lp(CX, CY, K):
    """U(c) at c = (CX + i CY)/2^K, CY > 0, at scale 2^-LP. With n_k = CX D - E_k 2^K, Y = CY D, M = D 2^K:
    4 DJ M^2 U = sum J_k (n_k^2 - Y^2) ln(n_k^2 + Y^2) - (2 ln M + 3) sum J_k (n_k^2 - Y^2) - 4 Y sum J_k n_k theta_k,
    theta_k = arg(n_k + iY) in (0, pi)."""
    D, DJ = PL["D"], PL["DJ"]
    Y = CY * D
    Y2 = Y * Y
    cxD = CX * D
    alo = ahi = B = clo = chi = 0
    for Ek, Jk in zip(PL["NE"], PL["NJ"]):
        n = cxD - (Ek << K)
        L = ln_int_lp(n * n + Y2)
        c = Jk * (n * n - Y2)
        if c > 0:
            alo += c * L[0]
            ahi += c * L[1]
        else:
            alo += c * L[1]
            ahi += c * L[0]
        B += c
        th = atan2_lp(Y, n)
        cn = Jk * n
        if cn > 0:
            clo += cn * th[0]
            chi += cn * th[1]
        else:
            clo += cn * th[1]
            chi += cn * th[0]
    M = D << K
    lm = ln_int_lp(M)
    xb = (lm[0] * B, lm[1] * B) if B >= 0 else (lm[1] * B, lm[0] * B)
    nlo = alo - 2 * xb[1] - 3 * B * (1 << LP) - 4 * Y * chi
    nhi = ahi - 2 * xb[0] - 3 * B * (1 << LP) - 4 * Y * clo
    den = 4 * DJ * M * M
    return nlo // den, cdiv(nhi, den)


def P_direct_lp(p, q):
    """sum_j b_j ln|Q_j(p/q)| at scale 2^-LP (exact polynomial values; p/q must not be a root)."""
    lo = hi = 0
    for Q, bb in zip(PL["A"], PL["b"]):
        if bb == 0:
            continue
        v = ev_int(Q, p, q)
        if v == 0:
            raise ValueError("a root of an auxiliary polynomial")
        t = ivq(ln_q_lp(abs(v), q ** (len(Q) - 1)), bb.numerator, bb.denominator)
        lo += t[0]
        hi += t[1]
    return lo, hi


def f_direct_pl(x):
    """f(x) = U(x) - sum_j b_j ln|Q_j(x)| at a rational x (Fraction), interval at scale 2^-P (direct sums)."""
    p, q = x.numerator, x.denominator
    u = U_real_lp(p, q)
    pp = P_direct_lp(p, q)
    return (u[0] - pp[1]) >> SH, cdiv(u[1] - pp[0], 1 << SH)


# ---- E4 ----

def pl_E4_one(j):
    """int log|Q_j| d mu = m log|lc| + sum over real roots U(r) + 2 sum over upper roots U(r), at scale 2^-P.
    A real root r in [c - rad, c + rad] costs |U(r) - U(c)| <= 5 rad sum_k |j_k| (|G2'(u)| <= 5 for |u| <= 6); a
    nonreal root in the disc D(c, rad) costs at most rad m/(Im c - rad) (|grad U| <= m/|Im z|)."""
    Q = PL["A"][j]
    m = PL["m"]
    lo = hi = 0
    if abs(Q[-1]) != 1:
        t = ivq(ln_int_lp(abs(Q[-1])), m.numerator, m.denominator)
        lo, hi = t
    lip = 5 * PL["jabs"]
    for rlo, rhi in PL["rroots"][j]:
        c, rad = (rlo + rhi) / 2, (rhi - rlo) / 2
        u = U_real_lp(c.numerator, c.denominator)
        w = cdiv((lip * rad).numerator << LP, (lip * rad).denominator) if rad else 0
        lo += u[0] - w
        hi += u[1] + w
    for CX, CY, rad in PL["croots"][j]:
        u = U_complex_lp(CX, CY, KC)
        cy = Fraction(CY, 1 << KC)
        bnd = rad * m / (cy - rad)
        w = cdiv(bnd.numerator << LP, bnd.denominator)
        lo += 2 * (u[0] - w)
        hi += 2 * (u[1] + w)
    return j, (lo >> SH, cdiv(hi, 1 << SH))


def _pl_E4_quad(j):
    """int log|Q_j| d mu in floats: brute-force quadrature of U at the centres of Q_j's root enclosures (self-test)."""
    nodes, weights = _gauss_legendre(20)
    cert = PL["cert"]
    Q = cert.A[j]
    v = float(cert.m) * math.log(abs(Q[-1]))
    for lo, hi in PL["rroots"][j]:
        v += pl_quad_U(cert, float((lo + hi) / 2), nodes, weights)
    for CX, CY, _ in PL["croots"][j]:
        v += 2 * pl_quad_U(cert, complex(CX / 2.0 ** KC, CY / 2.0 ** KC), nodes, weights)
    return j, v


def pl_E4_all(workers, fn=None):
    fn = fn or pl_E4_one
    jobs = list(range(len(PL["A"])))
    if workers <= 1:
        return [fn(j) for j in jobs]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(min(workers, len(jobs))) as pool:
        out = list(pool.imap_unordered(fn, sorted(jobs, key=lambda j: -len(PL["A"][j])), chunksize=1))
    return sorted(out)

# ---- E5: the per-cell model ----
#
# On support cell [e_i, e_{i+1}] write x = x0 + r sigma, sigma in [-1, 1]. Nodes and roots within LAM r of x0 are
# "near" and enter exactly on each piece; the others ("far") form S(sigma), a Taylor polynomial of degree NTP in
# sigma with explicit remainder bounds for S, S' and S'' (derivatives in sigma). For a far node, t = x0 - e_k,
# zeta = r/t: the coefficients of G2(t + r sigma) are G2(t), r t (ln|t| - 1), r^2 ln|t|/2 and
# (-1)^(n-1) r^2 zeta^(n-2)/(n (n-1) (n-2)) for n >= 3; Lagrange's remainder (|G2^(N+1)(v)| = (N-2)!/|v|^(N-1))
# is at most r^2 q^(N-1)/((N+1) N (N-1)), with q = r/(|t| - r), and r^2 q^(N-1)/(N (N-1)), r^2 q^(N-1)/(N-1) for
# the first two sigma-derivatives. For a far root c, w = x0 - c, z = r/w: -wt ln|w| - wt Re sum_n (-1)^(n-1) z^n
# sigma^n/n, remainders wt |z|^(N+1)/((N+1)(1-|z|)), wt |z|^(N+1)/(1-|z|) and
# wt |z|^(N+1)(N/(1-|z|) + |z|/(1-|z|)^2); the distance from the true root to the centre c adds
# wt rad/(d - rad), wt r rad/(d - rad)^2 and 3 wt r^2 rad/(d - rad)^3, d = |w| - r.
# Every near node has |x - e_k| < 1 on the cell (checked), where G2'' = ln|u| < 0: j_k G2 is concave for j_k > 0
# (it enters through its chord) and convex for j_k < 0 (through its tangent). A near real root outside a piece
# enters through its tangent (-ln|x - r| is convex), one that may lie in the piece through -wt ln(largest distance).

NTP = 24                        # degree of the per-cell polynomial (in sigma)
LAM = 4                         # near: |x0 - e_k| < LAM r (so |zeta| <= 1/4 and q <= 1/3 for far nodes)
KT = 24                         # extra binary digits of a tangent point


def iv_mul_lp(a, b):
    c = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(c) >> LP, -((-max(c)) >> LP)


def fix_lp(fr):
    return (fr.numerator << LP) // fr.denominator, cdiv(fr.numerator << LP, fr.denominator)


def up32(fr):
    """A rational upper bound of fr >= 0 with denominator 2^32."""
    return Fraction(cdiv(fr.numerator << 32, fr.denominator), 1 << 32)


def lo32(fr):
    """A rational lower bound of fr with denominator 2^32."""
    return Fraction((fr.numerator << 32) // fr.denominator, 1 << 32)


def sqrt_lo32(fr):
    """A rational lower bound of sqrt(fr), fr >= 0, denominator 2^32."""
    return Fraction(isqrt((fr.numerator << 64) // fr.denominator), 1 << 32)


def sqrt_up32(fr):
    return Fraction(isqrt(cdiv(fr.numerator << 64, fr.denominator)) + 1, 1 << 32)


class PLCell:
    def __init__(self, i):
        D, DJ, EI = PL["D"], PL["DJ"], PL["EI"]
        Ea, Eb = EI[i], EI[i + 1]
        X0, R = (Ea + Eb) // 2, (Eb - Ea) // 2
        assert 2 * X0 == Ea + Eb and 2 * R == Eb - Ea and R > 0
        self.i, self.X0, self.R = i, X0, R
        x0, r = Fraction(X0, D), Fraction(R, D)
        self.x0, self.r = x0, r
        lim = LAM * R
        self.nn, ns, Js = [], [], []
        for Ek, Jk in zip(PL["NE"], PL["NJ"]):
            n = X0 - Ek
            if -lim < n < lim:
                self.nn.append((Ek, Jk))
            else:
                ns.append(n)
                Js.append(Jk)
        self.supported = (LAM + 1) * R < D          # every near node within distance 1 of the cell: G2 concave
        N = NTP
        ONEL = 1 << LP
        # ---- far nodes: the logarithmic coefficients (exact integer sums of table logarithms) ----
        LNT, W = PL["LNT"], PL["LNTW"]
        Ls = [LNT[n if n > 0 else -n] for n in ns]
        lnD = PL["LND"]

        def lsum(cs):
            """Interval of sum_k cs_k ln|t_k| at scale 2^-LP, ln|t_k| = ln|n_k| - ln D, ln|n_k| in [L, L + W]."""
            S = sum(map(mul, cs, Ls))
            Bs = sum(cs)
            Ps = sum(c for c in cs if c > 0)
            xb = (lnD[0] * Bs, lnD[1] * Bs) if Bs >= 0 else (lnD[1] * Bs, lnD[0] * Bs)
            return (S + W * (Bs - Ps) - xb[1], S + W * Ps - xb[0]), Bs

        l0, B0 = lsum([J * n * n for J, n in zip(Js, ns)])
        l1, B1 = lsum([J * n for J, n in zip(Js, ns)])
        l2, _ = lsum(Js)
        DD = DJ * D * D
        co = [None] * (N + 1)
        co[0] = ivq((2 * l0[0] - 3 * B0 * ONEL, 2 * l0[1] - 3 * B0 * ONEL), 1, 4 * DD)
        co[1] = ivq((l1[0] - B1 * ONEL, l1[1] - B1 * ONEL), R, DD)
        co[2] = ivq(l2, R * R, 2 * DD)
        # ---- far nodes: power sums sum_k J_k zeta_k^p, zeta_k = R/n_k (outward fixed point) ----
        Zl = [(R << LP) // (n if n > 0 else -n) for n in ns]
        Zh = [z + 1 for z in Zl]
        Jep = [J if J > 0 else 0 for J in Js]
        Jen = [J if J < 0 else 0 for J in Js]
        Jo = [J if n > 0 else -J for J, n in zip(Js, ns)]
        Jop = [J if J > 0 else 0 for J in Jo]
        Jon = [J if J < 0 else 0 for J in Jo]
        pl, ph = Zl, Zh
        for p in range(1, N - 1):
            jp, jn = (Jop, Jon) if p % 2 else (Jep, Jen)
            P_ = (sum(map(mul, jp, pl)) + sum(map(mul, jn, ph)), sum(map(mul, jp, ph)) + sum(map(mul, jn, pl)))
            n = p + 2
            co[n] = ivq(P_, (1 if n % 2 else -1) * R * R, DD * n * (n - 1) * (n - 2))
            if p < N - 2:
                pl = [(a * z) >> LP for a, z in zip(pl, Zl)]
                ph = [-((-(a * z)) >> LP) for a, z in zip(ph, Zh)]
        # remainders of the far nodes, in three distance tiers
        s1 = s2 = s3 = 0
        for J, n in zip(Js, ns):
            an = n if n > 0 else -n
            aJ = J if J > 0 else -J
            if an >= 64 * R:
                s3 += aJ
            elif an >= 8 * R:
                s2 += aJ
            else:
                s1 += aJ
        tiers = s1 * Fraction(1, 3) ** (N - 1) + s2 * Fraction(1, 7) ** (N - 1) + s3 * Fraction(1, 63) ** (N - 1)
        base = Fraction(R * R, DD) * tiers
        Rm = [qiv(base / ((N + 1) * N * (N - 1)), LP)[1], qiv(base / (N * (N - 1)), LP)[1], qiv(base / (N - 1), LP)[1]]
        # ---- roots ----
        acc = [(0, 0)] * (N + 1)
        self.nr, self.nc = [], []
        for wt, c, rad, j in PL["rr"]:
            w = x0 - c
            aw = abs(w)
            if aw < LAM * r:
                self.nr.append((wt, c, rad))
                continue
            t = ivq(ln_q_lp(aw.numerator, aw.denominator), -wt.numerator, wt.denominator)
            acc[0] = (acc[0][0] + t[0], acc[0][1] + t[1])
            zi = fix_lp(r / w)
            pw = zi
            for n in range(1, N + 1):
                if n > 1:
                    pw = iv_mul_lp(pw, zi)
                t = ivq(pw, (-1 if n % 2 else 1) * wt.numerator, wt.denominator * n)
                acc[n] = (acc[n][0] + t[0], acc[n][1] + t[1])
            self._root_rem(Rm, wt, up32(r / aw), rad, lo32(aw - r))
        for wt, CX, CY, rad, j in PL["cr"]:
            cx, cy = Fraction(CX, 1 << KC), Fraction(CY, 1 << KC)
            a = x0 - cx
            g = a * a + cy * cy
            if g < (LAM * r) ** 2:
                self.nc.append((wt, cx, cy, rad))
                continue
            t = ivq(ln_q_lp(g.numerator, g.denominator), -wt.numerator, 2 * wt.denominator)
            acc[0] = (acc[0][0] + t[0], acc[0][1] + t[1])
            zr, zim = fix_lp(r * a / g), fix_lp(r * cy / g)
            pr, pim = zr, zim
            for n in range(1, N + 1):
                if n > 1:
                    a1, a2 = iv_mul_lp(pr, zr), iv_mul_lp(pim, zim)
                    b1, b2 = iv_mul_lp(pr, zim), iv_mul_lp(pim, zr)
                    pr, pim = (a1[0] - a2[1], a1[1] - a2[0]), (b1[0] + b2[0], b1[1] + b2[1])
                t = ivq(pr, (-1 if n % 2 else 1) * wt.numerator, wt.denominator * n)
                acc[n] = (acc[n][0] + t[0], acc[n][1] + t[1])
            self._root_rem(Rm, wt, sqrt_up32(r * r / g), rad, sqrt_lo32(g) - r)
        lc = PL["lc"]
        acc[0] = (acc[0][0] + lc[0], acc[0][1] + lc[1])
        coef = []
        for n in range(N + 1):
            lo, hi = co[n][0] + acc[n][0], co[n][1] + acc[n][1]
            coef.append((lo >> SH, cdiv(hi, 1 << SH)))
        self.coef = coef
        self.d1 = [(n * coef[n][0], n * coef[n][1]) for n in range(1, N + 1)]
        self.d2 = [(n * (n - 1) * coef[n][0], n * (n - 1) * coef[n][1]) for n in range(2, N + 1)]
        self.R0, self.R1, self.R2 = (cdiv(x, 1 << SH) + 1 for x in Rm)
        self.Rlp = Rm
        # floats for the heuristics (they choose points only)
        self.cf = [ivfl(c) for c in coef]
        self.rf = float(r)
        self.nnf = [(float(Fraction(Ek - X0, R)), Jk / float(DJ)) for Ek, Jk in self.nn]
        self.nrf = [(float(wt), float((c - x0) / r)) for wt, c, _ in self.nr]
        self.ncf = [(float(wt), float((cx - x0) / r), float(cy / r)) for wt, cx, cy, _ in self.nc]
        self.nn_cav = [(Ek, Jk) for Ek, Jk in self.nn if Jk > 0]      # j G2 concave: chords
        self.nn_cvx = [(Ek, Jk) for Ek, Jk in self.nn if Jk < 0]      # j G2 convex: tangents
        self.cache, self.ccache, self.lnM = {}, {}, {}

    def _root_rem(self, Rm, wt, qz, rad, d):
        """Add a far root's remainder bounds (value, sigma-derivative, second sigma-derivative) to Rm, integers at
        scale 2^-LP rounded up: qz >= |z|, rad >= the centre's distance to the root, 0 < d <= min over the cell of
        |x - centre|."""
        N = NTP
        if qz >= 1 or d <= rad:
            raise ValueError("a far root too close to the cell")
        a = wt * qz ** (N + 1)
        r = self.r
        for idx, v in enumerate((a / ((N + 1) * (1 - qz)) + wt * rad / (d - rad),
                                 a / (1 - qz) + wt * r * rad / (d - rad) ** 2,
                                 a * (N / (1 - qz) + qz / (1 - qz) ** 2) + 3 * wt * r * r * rad / (d - rad) ** 3)):
            Rm[idx] += cdiv(v.numerator << LP, v.denominator)

    # ---- exact pieces of f ----
    def xs(self, num, k):
        """x at sigma = num/2^k."""
        return Fraction((self.X0 << k) + self.R * num, PL["D"] << k)

    def _lnM(self, k):
        v = self.lnM.get(k)
        if v is None:
            v = self.lnM[k] = lnM_lp(k)
        return v

    def _lu(self, Ek, num, k):
        nu = ((self.X0 - Ek) << k) + self.R * num          # D 2^k (x - e_k)
        if nu == 0:
            return 0, None
        L = ln_int_lp(nu if nu > 0 else -nu)
        lm = self._lnM(k)
        return nu, (L[0] - lm[1], L[1] - lm[0])

    def node_val(self, Ek, Jk, num, k):
        """j_k G2(x - e_k) at sigma = num/2^k, scale 2^-P."""
        nu, lu = self._lu(Ek, num, k)
        if nu == 0:
            return 0, 0
        s2 = nu * nu
        M = PL["D"] << k
        v = ivq((s2 * (2 * lu[0] - 3 * (1 << LP)), s2 * (2 * lu[1] - 3 * (1 << LP))), Jk, 4 * M * M * PL["DJ"])
        return v[0] >> SH, cdiv(v[1], 1 << SH)

    def node_vd(self, Ek, Jk, num, k):
        """(value, sigma-derivative r j_k G2'(u)) of j_k G2(x - e_k) at sigma = num/2^k, scale 2^-P."""
        nu, lu = self._lu(Ek, num, k)
        if nu == 0:
            return (0, 0), (0, 0)
        s2 = nu * nu
        M = PL["D"] << k
        v = ivq((s2 * (2 * lu[0] - 3 * (1 << LP)), s2 * (2 * lu[1] - 3 * (1 << LP))), Jk, 4 * M * M * PL["DJ"])
        d = ivq((lu[0] - (1 << LP), lu[1] - (1 << LP)), Jk * self.R * nu, PL["DJ"] * PL["D"] * M)
        return (v[0] >> SH, cdiv(v[1], 1 << SH)), (d[0] >> SH, cdiv(d[1], 1 << SH))

    def conc_at(self, num, k):
        """Sum of the lower ends of the concave near terms (j_k > 0) at sigma = num/2^k (cached)."""
        while k > 0 and num % 2 == 0:
            num //= 2
            k -= 1
        key = (num, k)
        v = self.ccache.get(key)
        if v is None:
            v = self.ccache[key] = sum(self.node_val(Ek, Jk, num, k)[0] for Ek, Jk in self.nn_cav)
        return v

    def S_pt(self, num, k):
        lo, hi = horner_point(self.coef, num, k)
        return lo - self.R0, hi + self.R0

    def S_at(self, num, k):
        while k > 0 and num % 2 == 0:
            num //= 2
            k -= 1
        key = (num, k)
        v = self.cache.get(key)
        if v is None:
            v = self.cache[key] = self.S_pt(num, k)
        return v

    def dS_pt(self, num, k):
        lo, hi = horner_point(self.d1, num, k)
        return lo - self.R1, hi + self.R1

    def d2S_rng(self, n1, n2, k):
        lo, hi = horner_range(self.d2, n1, n2, k)
        return lo - self.R2, hi + self.R2

    @staticmethod
    def _neg_wt_ln_lo(wt, d):
        """Lower end of -wt ln d (d > 0 rational), scale 2^-P."""
        L = ln_q_lp(d.numerator, d.denominator)
        return ivq(L, -wt.numerator, wt.denominator)[0] >> SH

    @staticmethod
    def _neg_wt_ln_hi(wt, d):
        L = ln_q_lp(d.numerator, d.denominator)
        return cdiv(ivq(L, -wt.numerator, wt.denominator)[1], 1 << SH)

    def nc_lo(self, wt, cx, cy, rad, xa, xb):
        """Lower bound of -wt ln|x - alpha| on [xa, xb] for a nonreal near root: the largest distance."""
        dd = max((xa - cx) ** 2, (xb - cx) ** 2) + cy * cy
        return self._neg_wt_ln_lo(wt, sqrt_up32(dd) + rad)

    def tangent(self, tn, tk, convex, troots):
        """Lower end of the tangent group's value at t = tn/2^tk, and the interval of its sigma-derivative."""
        Phi = Dlo = Dhi = 0
        for Ek, Jk in convex:
            v, d = self.node_vd(Ek, Jk, tn, tk)
            Phi += v[0]
            Dlo += d[0]
            Dhi += d[1]
        if troots:
            x = self.xs(tn, tk)
            for wt, c, rad in troots:
                Phi += self._neg_wt_ln_lo(wt, abs(x - c) + rad)
                v1, v2 = -wt * self.r / (x - c + rad), -wt * self.r / (x - c - rad)
                Dlo += qiv(min(v1, v2))[0]
                Dhi += qiv(max(v1, v2))[1]
        return Phi, Dlo, Dhi

    def piece(self, j, s, goal, hint=None):
        """Certified lower bound of f on sigma in [na, nb]/2^s, na = 2j - 2^s, nb = na + 2. Returns (lb, t)."""
        na = 2 * j - (1 << s)
        nb = na + 2
        SA, SB = self.S_at(na, s), self.S_at(nb, s)
        conc_a, conc_b = self.conc_at(na, s), self.conc_at(nb, s)
        convex = self.nn_cvx
        termw = 0
        troots = []
        if self.nr or self.nc:
            xa, xb = self.xs(na, s), self.xs(nb, s)
            for wt, c, rad in self.nr:
                if c - rad <= xb and c + rad >= xa:
                    termw += self._neg_wt_ln_lo(wt, max(xb - c + rad, c + rad - xa))
                else:
                    troots.append((wt, c, rad))
            for wt, cx, cy, rad in self.nc:
                termw += self.nc_lo(wt, cx, cy, rad, xa, xb)
        d2 = self.d2S_rng(na, nb, s)
        M2 = max(d2[1], 0)
        lossA = cdiv(M2, 1 << (2 * s + 1))            # (width/2)^2 / 2 sup S'', width = 2^(1-s)
        LA, LB = SA[0] + conc_a, SB[0] + conc_b
        best, tbest = None, None
        if not convex and not troots:
            best = min(LA, LB) - lossA + termw
            if best >= goal:
                return best, None
        lo_s, hi_s = Fraction(na, 1 << s), Fraction(nb, 1 << s)
        tries = [hint, None] if hint is not None and lo_s <= hint <= hi_s else [None]
        m2 = min(d2[0], 0)
        for t in tries:
            if t is None:
                t = self.choose_t(na, nb, s)
            tn, tk = dyadic(t)
            Phi, Dlo, Dhi = self.tangent(tn, tk, convex, troots)
            K = max(tk, s)
            da = (na << (K - s)) - (tn << (K - tk))       # (sigma_a - t) 2^K <= 0
            db = (nb << (K - s)) - (tn << (K - tk))       # (sigma_b - t) 2^K >= 0
            if convex or troots:
                lbA = min(LA + Phi + ((Dhi * da) >> K), LB + Phi + ((Dlo * db) >> K)) - lossA + termw
                if best is None or lbA > best:
                    best, tbest = lbA, t
                if best >= goal:
                    return best, tbest
            St = self.S_pt(tn, tk)
            dSt = self.dS_pt(tn, tk)
            ca = (m2 * da * da) >> (2 * K + 1)
            cb = (m2 * db * db) >> (2 * K + 1)
            lbB = min(conc_a + St[0] + Phi + (((dSt[1] + Dhi) * da) >> K) + ca,
                      conc_b + St[0] + Phi + (((dSt[0] + Dlo) * db) >> K) + cb) + termw
            if best is None or lbB > best:
                best, tbest = lbB, t
            if best >= goal:
                return best, tbest
        return best, tbest

    def f_upper(self, num, k):
        """Rigorous upper bound of f at sigma = num/2^k (None at a possible root)."""
        v = self.S_pt(num, k)[1]
        for Ek, Jk in self.nn:
            v += self.node_val(Ek, Jk, num, k)[1]
        if self.nr or self.nc:
            x = self.xs(num, k)
            for wt, c, rad in self.nr:
                d = abs(x - c) - rad
                if d <= 0:
                    return None
                v += self._neg_wt_ln_hi(wt, d)
            for wt, cx, cy, rad in self.nc:
                d = sqrt_lo32((x - cx) ** 2 + cy * cy) - rad
                if d <= 0:
                    return None
                v += self._neg_wt_ln_hi(wt, d)
        return v

    def f_interval(self, num, k):
        lo, hi = self.S_pt(num, k)
        for Ek, Jk in self.nn:
            a = self.node_val(Ek, Jk, num, k)
            lo += a[0]
            hi += a[1]
        x = self.xs(num, k)
        for wt, c, rad in self.nr:
            lo += self._neg_wt_ln_lo(wt, abs(x - c) + rad)
            hi += self._neg_wt_ln_hi(wt, abs(x - c) - rad)
        for wt, cx, cy, rad in self.nc:
            dd = (x - cx) ** 2 + cy * cy
            lo += self._neg_wt_ln_lo(wt, sqrt_up32(dd) + rad)
            hi += self._neg_wt_ln_hi(wt, sqrt_lo32(dd) - rad)
        return lo, hi

    # ---- floats (heuristics only: they choose tangent points and goals) ----
    def f_float(self, sg):
        acc = 0.0
        for c in reversed(self.cf):
            acc = acc * sg + c
        rf = self.rf
        for sk, jf in self.nnf:
            u = rf * (sg - sk)
            if u != 0.0:
                acc += jf * u * u * (0.5 * math.log(abs(u)) - 0.75)
        for wf, sc in self.nrf:
            d = abs(sg - sc)
            if d == 0.0:
                return float("inf")
            acc -= wf * math.log(rf * d)
        for wf, sr, si in self.ncf:
            acc -= wf * 0.5 * math.log(rf * rf * ((sg - sr) ** 2 + si * si))
        return acc

    def _golden(self, a, b, iters=40):
        gr = (math.sqrt(5) - 1) / 2
        x1, x2 = b - gr * (b - a), a + gr * (b - a)
        f1, f2 = self.f_float(x1), self.f_float(x2)
        for _ in range(iters):
            if f1 < f2:
                b, x2, f2 = x2, x1, f1
                x1 = b - gr * (b - a)
                f1 = self.f_float(x1)
            else:
                a, x1, f1 = x1, x2, f2
                x2 = a + gr * (b - a)
                f2 = self.f_float(x2)
        return (x1, f1) if f1 < f2 else (x2, f2)

    def choose_t(self, na, nb, s):
        """A dyadic point of [na, nb]/2^s near the minimiser of f there (floats; any point gives a valid bound)."""
        a, b = na / float(1 << s), nb / float(1 << s)
        cands = [(self.f_float(a), a), (self.f_float(b), b)]
        x, fx_ = self._golden(a, b)
        cands.append((fx_, x))
        v = min(cands)[1]
        k = s + KT
        t = Fraction(int(round(v * (1 << k))), 1 << k)
        return min(max(t, Fraction(na, 1 << s)), Fraction(nb, 1 << s))

    def float_scan(self):
        """Float estimate of min f on the cell and where (for the goal and the upper bound)."""
        pts = [-1 + k / 128.0 for k in range(257)]
        for e_ in range(3, 53):
            pts += [-1 + 2.0 ** -e_, 1 - 2.0 ** -e_]
        vals = sorted((self.f_float(p), p) for p in pts)
        best = vals[0]
        for v, p in vals[:6]:
            w = max(min(1 + p, 1 - p) / 2, 1.0 / 128)
            x, fx_ = self._golden(max(-1.0, p - w), min(1.0, p + w), 60)
            if fx_ < best[0]:
                best = (fx_, x)
        return best

    def bisect(self, goal, smax=50, max_pieces=200000):
        stack = [(0, 0, None)]
        npieces = 0
        minlb, minat = None, None
        fails, maxs = [], 0
        while stack:
            j, s, hint = stack.pop()
            lb, t = self.piece(j, s, goal, hint)
            npieces += 1
            if lb >= goal:
                if minlb is None or lb < minlb:
                    minlb, minat = lb, (j, s)
                maxs = max(maxs, s)
                continue
            if s >= smax or npieces > max_pieces:
                fails.append((j, s, lb))
                continue
            h = t if t is not None else hint
            stack.append((2 * j + 1, s + 1, h))
            stack.append((2 * j, s + 1, h))
        return {"pieces": npieces, "minlb": minlb, "minat": minat, "fails": fails, "maxdepth": maxs}


def _pl_job(args):
    """Model, float scan and bisection of one support cell. goal None: certify f >= 0.9 x (float minimum) on the cell,
    falling back to 0.5 x and to 0."""
    i, goal = args
    t0 = time.time()
    cm = PLCell(i)
    if not cm.supported:
        return {"cell": i, "unsupported": True, "fails": [(0, 0, None)], "pieces": 0, "minlb": None, "maxdepth": 0,
                "v": None, "sigma": None, "ub": None, "goal": goal, "time": time.time() - t0, "build": 0.0}
    tb = time.time() - t0
    v, sg = cm.float_scan()
    if goal is not None:
        goals = [goal]
    elif v > 0:
        goals = [qiv(Fraction(v) * Fraction(9, 10))[0], qiv(Fraction(v) / 2)[0], 0]
    else:
        goals = [0]
    for g in goals:
        res = cm.bisect(g)
        if not res["fails"]:
            break
    sgf = Fraction(min(max(sg, -1.0), 1.0))
    num, k = dyadic(sgf)
    res.update(cell=i, v=v, sigma=sgf, x=cm.x0 + cm.r * sgf, ub=cm.f_upper(num, k), goal=g, build=tb,
               time=time.time() - t0, unsupported=False)
    return res


def pl_run_cells(jobs, workers):
    if workers <= 1:
        return [_pl_job(a) for a in jobs]
    import multiprocessing as mp
    order = sorted(jobs, key=lambda a: (a[0] * 7919) % 104729)
    with mp.get_context("fork").Pool(workers) as pool:
        out = list(pool.imap_unordered(_pl_job, order, chunksize=2))
    return sorted(out, key=lambda r: r["cell"])


# ---- E5 driver ----

def pl_check_E5(rep, cert, workers):
    """Every support cell: f >= 0.9 x (its float minimum) certified (else 0.5 x, else 0); then the cells whose
    certified bound lies below (1 - 1e-8) x the smallest float minimum are re-run with that goal. The enclosure of
    min_S f: lower end the smallest certified piece bound, upper end f at the float minimiser from the direct sums."""
    t0 = time.time()
    res = pl_run_cells([(i, None) for i in cert.support], workers)
    uns = [r["cell"] for r in res if r["unsupported"]]
    rep.check("E5 every support cell has all its near nodes within distance 1 (where G2 is concave)", not uns,
              "%d cells outside" % len(uns))
    if uns:
        return None
    rmin = min(res, key=lambda r: r["v"])
    vmin = rmin["v"]
    rep.info("E5 float scan (heuristic only)", "smallest f found %.7e in cell %d [%s, %s] at x = %.9f" % (
        vmin, rmin["cell"], cert.e[rmin["cell"]], cert.e[rmin["cell"] + 1], float(rmin["x"])))
    failed = [r for r in res if r["fails"]]
    fb = sum(1 for r in res if r["v"] > 0 and r["goal"] < qiv(Fraction(r["v"]) * Fraction(9, 10))[0])
    npieces = sum(r["pieces"] for r in res)
    build = sum(r["build"] for r in res)
    tot = sum(r["time"] for r in res)
    lb1 = min(r["minlb"] for r in res if r["minlb"] is not None)
    rep.check("E5 f >= 0.9 x (the cell's float minimum) certified on every support cell (%d cells, %d pieces, max "
              "depth %d)" % (len(res), npieces, max(r["maxdepth"] for r in res)), not failed,
              "%d cells failed, %d needed the weaker goal 0.5 x or 0; models %.0f s and bisection %.0f s summed over "
              "the cells' jobs, %.0f s wall" % (len(failed), fb, build, tot - build, time.time() - t0))
    rep.info("E5 smallest certified piece bound, first pass", sci_floor(Fraction(lb1, ONE), 8))
    # tighten the cells whose certified bound lies below (1 - 1e-8) x the smallest float minimum
    target = qiv(Fraction(vmin) * (1 - Fraction(1, 10 ** 8)))[0]
    redo = [r["cell"] for r in res if not r["fails"] and r["minlb"] < target]
    t1 = time.time()
    if redo:
        byc = {r["cell"]: r for r in pl_run_cells([(i, target) for i in redo], workers)}
        tight_fail = [i for i in redo if byc[i]["fails"]]
        res = [byc[r["cell"]] if r["cell"] in byc and not byc[r["cell"]]["fails"] else r for r in res]
    else:
        tight_fail = []
    LB = min(r["minlb"] for r in res if r["minlb"] is not None)
    rep.check("E5 re-run of the %d cells below (1 - 1e-8) x the smallest float minimum, at that goal" % len(redo),
              not tight_fail, "%.1f s; cells %s" % (time.time() - t1, ", ".join(str(i) for i in redo[:12])))
    xs = rmin["x"]                                # a dyadic point of the closed cell, so of S
    UBd = f_direct_pl(xs)
    UB = UBd[1]
    rep.check("E5 min over S of f enclosed", LB <= UB,
              "min_S f in [%s, %s]; upper end: f at x = %.12f (cell %d, [%s, %s]) from the direct sums over all "
              "nodes and exact polynomial values; lower end: the smallest certified piece bound over all cells" % (
                  sci_floor(Fraction(LB, ONE), 8), sci_ceil(Fraction(UB, ONE), 8), float(xs), rmin["cell"],
                  cert.e[rmin["cell"]], cert.e[rmin["cell"] + 1]))
    ok = not failed and LB > 0
    rep.check("E5: f > 0 on S (every closed support cell, cell ends included)", ok,
              "certified f >= %s" % sci_floor(Fraction(LB, ONE), 6))
    return {"LB": LB, "UB": UB, "x": xs, "rmin": rmin, "res": res, "ok": ok}


# ---- the exact conditions for this format ----

def p298_E2(rep, A, b, quiet=False):
    """Every Q_j with b_j > 0: primitive, irreducible (mod-p factor-degree sieve), positive leading coefficient, no
    two equal. Returns the list of failures."""
    bad = []
    seen = set()
    for j, Q in enumerate(A):
        if b[j] == 0:
            continue
        g = 0
        for c in Q:
            g = math.gcd(g, c)
        ok = len(Q) >= 2 and Q[-1] != 0
        sieve, sdet = irreducible_sieve(Q, maxp=2000) if ok else (False, "degree < 1")
        nreal, sqf = real_root_count(Q) if ok else (0, False)
        dup = tuple(Q) in seen or tuple(-c for c in Q) in seen
        seen.add(tuple(Q))
        ok = ok and Q[-1] > 0 and g == 1 and sieve and sqf and not dup
        if not ok:
            bad.append((j, "gcd %d, lc %d, %s%s" % (g, Q[-1], sdet, ", repeated" if dup else "")))
        if not quiet:
            rep.check("E2 Q%d (degree %d, lc %d): primitive, irreducible, positive lc, distinct" % (
                j, len(Q) - 1, Q[-1]), ok, "%s; %d of %d roots real" % (sdet, nreal, len(Q) - 1))
    return bad


def p298_roots(rep, cert):
    """Real roots by exact Sturm isolation (width 2^-112), all roots by Gershgorin discs; the two agree, and the
    discs off the real axis hold exactly the nonreal roots. Returns (rroots, croots, discs, nreal)."""
    rroots, croots, alld, nreals = [], [], [], []
    okall = True
    worst = Fraction(0)
    for j, Q in enumerate(cert.A):
        nreal, _ = real_root_count(Q)
        rr = isolate_roots(Q)
        discs, prec = certified_discs(Q)
        ok, up, down, mid = split_discs(discs, nreal)
        matched = len(rr) == nreal and all(sum(1 for d in mid if d[0] - d[2] <= hi and d[0] + d[2] >= lo) == 1
                                           for lo, hi in rr)
        conj = all(any(abs(u[0] - w[0]) <= u[2] + w[2] and abs(u[1] + w[1]) <= u[2] + w[2] for w in down) for u in up)
        okall &= ok and matched and conj
        worst = max([worst] + [d[2] for d in discs])
        rroots.append(rr)
        croots.append([round_disc(*u) for u in up])
        alld.append(discs)
        nreals.append(nreal)
    rep.check("roots: every Q_j has deg Q_j disjoint Gershgorin discs (one root each); exactly nreal of them meet the "
              "real axis, each holding one Sturm-isolated real root; the others pair up under conjugation", okall,
              "%d real roots (exact Sturm counts, enclosures of width <= 2^-112), %d nonreal in %d conjugate pairs; "
              "largest disc radius 2^%d" % (sum(nreals), sum(len(Q) - 1 for Q in cert.A) - sum(nreals),
                                             sum(len(c) for c in croots), worst.numerator.bit_length() -
                                             worst.denominator.bit_length()))
    return rroots, croots, alld, nreals


def vieta_pl(Q, rr, discs):
    """e1 = -a_{n-1}/a_n, e2 = a_{n-2}/a_n and prod |root|^2 = (a_0/a_n)^2 for the roots in the real enclosures and
    the discs off the axis, within the propagated radii (exact rationals)."""
    n = len(Q) - 1
    zs = [((lo + hi) / 2, Fraction(0), (hi - lo) / 2) for lo, hi in rr]
    zs += [d for d in discs if d[1] - d[2] > 0 or d[1] + d[2] < 0]
    if len(zs) != n:
        return False
    an = Q[n]
    R = sum(z[2] for z in zs)
    ok = abs(sum(z[0] for z in zs) + Fraction(Q[n - 1], an)) <= R and abs(sum(z[1] for z in zs)) <= R
    if n >= 2:
        A_ = sum(abs(z[0]) + abs(z[1]) for z in zs)
        ex, ey = sum(z[0] for z in zs), sum(z[1] for z in zs)
        px = sum(z[0] * z[0] - z[1] * z[1] for z in zs)
        py = sum(2 * z[0] * z[1] for z in zs)
        e2x, e2y = (ex * ex - ey * ey - px) / 2, (2 * ex * ey - py) / 2
        B = (A_ + R) ** 2 - A_ ** 2
        ok &= abs(e2x - Fraction(Q[n - 2], an)) <= B and abs(e2y) <= B
    plo, phi = Fraction(1), Fraction(1)
    for x, y, r in zs:
        m2 = x * x + y * y
        bb = abs(x) + abs(y)
        plo *= max(Fraction(0), m2 - 2 * r * bb)
        phi *= m2 + 2 * r * bb + r * r
    t = Fraction(Q[0], an) ** 2
    return ok and plo <= t <= phi


def _wlog(x, p, q, rp, rq, nodes, weights):
    """int_p^q log|x - y| rho(y) dy for rho linear from rp (at p) to rq (at q), x outside (p, q), in floats:
    Gauss-Legendre in the distance from the end nearer to x, on a mesh graded geometrically towards it."""
    if abs(x - p) <= abs(x - q):
        near, far, rn, rf_ = p, q, rp, rq
    else:
        near, far, rn, rf_ = q, p, rq, rp
    L = abs(far - near)
    d0 = abs(x - near)
    sl = (rf_ - rn) / L
    cuts, w = [L], L
    while w > 1e-30 * L and w > 0.25 * d0:
        w /= 2
        cuts.append(w)
    pts = sorted(set(cuts + [0.0]))
    s = 0.0
    if d0 == 0.0:
        eps = pts[1]
        le = math.log(eps)
        s += rn * eps * (le - 1) + sl * (eps * eps * le / 2 - eps * eps / 4)
        pts = pts[1:]
    for u, v in zip(pts, pts[1:]):
        h = (v - u) / 2
        for z, wt in zip(nodes, weights):
            t = u + h * (z + 1)
            s += wt * h * math.log(d0 + t) * (rn + sl * t)
    return s


def pl_quad_U(cert, z, nodes, weights):
    """U(z) by brute-force quadrature over every support cell, in floats (self-test only)."""
    tot = 0.0
    if isinstance(z, complex):
        zr, zi2 = z.real, z.imag ** 2
        for k in cert.support:
            a, b = float(cert.e[k]), float(cert.e[k + 1])
            ha, hb = float(cert.h[k]), float(cert.h[k + 1])
            h = (b - a) / 2
            for t, wt in zip(nodes, weights):
                y = a + h * (t + 1)
                tot += wt * h * 0.5 * math.log((y - zr) ** 2 + zi2) * (ha + (hb - ha) * (t + 1) / 2)
        return tot
    for k in cert.support:
        a, b = float(cert.e[k]), float(cert.e[k + 1])
        ha, hb = float(cert.h[k]), float(cert.h[k + 1])
        if a < z < b:
            hz = ha + (hb - ha) * (z - a) / (b - a)
            tot += _wlog(z, a, z, ha, hz, nodes, weights) + _wlog(z, z, b, hz, hb, nodes, weights)
        else:
            tot += _wlog(z, a, b, ha, hb, nodes, weights)
    return tot


def pl_selftests(rep, cert, rng, alld):
    """Our logarithms and arctangents, the jump identities, the closed-form U against quadrature, Vieta."""
    getcontext().prec = 70
    bad = 0
    ns = [rng.randrange(1, 10 ** rng.randrange(1, 80)) for _ in range(600)] + list(range(1, 300))
    for n in ns:
        lo, hi = ln_int(n)
        ref = Fraction(Decimal(n).ln())
        tol = Fraction(1, 10 ** 60)
        if not Fraction(lo, ONE) - tol <= ref <= Fraction(hi, ONE) + tol:
            bad += 1
    rep.check("self-test: our ln encloses the correctly rounded decimal ln (70 digits) at %d integers" % len(ns),
              bad == 0)
    # arctangent: 4 atan 1 and 4 (atan 1/2 + atan 1/3) enclose pi (pi from Machin's formula, exact series), the
    # addition formula atan a + atan b = atan((a + b)/(1 - ab)) at random rationals (both sides of atan 1), math.atan2
    one = 1 << LP
    t1 = atan_q(1, 1)
    a2, a3 = atan_q(1, 2), atan_q(1, 3)
    ok = 4 * t1[0] <= PI_LP[1] and PI_LP[0] <= 4 * t1[1]
    ok &= 4 * (a2[0] + a3[0]) <= PI_LP[1] and PI_LP[0] <= 4 * (a2[1] + a3[1])
    worst_w, worst_f = 0, 0.0
    for _ in range(400):
        q1, q2 = rng.randrange(2, 10 ** rng.randrange(2, 40)), rng.randrange(2, 10 ** rng.randrange(2, 40))
        p1, p2 = rng.randrange(0, q1), rng.randrange(0, q2)          # a = p1/q1, b = p2/q2 in [0, 1), ab < 1
        u, v = atan_q(p1, q1), atan_q(p2, q2)
        w = atan_q(p1 * q2 + p2 * q1, q1 * q2 - p1 * p2)
        ok &= u[0] + v[0] <= w[1] and w[0] <= u[1] + v[1]
        worst_w = max(worst_w, u[1] - u[0], w[1] - w[0])
        worst_f = max(worst_f, abs((u[0] + u[1]) / 2 / one - math.atan2(p1, q1)))
        n_ = rng.randrange(-10 ** 30, 10 ** 30)
        th = atan2_lp(q2, n_)
        worst_f = max(worst_f, abs((th[0] + th[1]) / 2 / one - math.atan2(q2, n_)))
        ok &= 0 < th[0] <= th[1] < PI_LP[1]
    rep.check("self-test: our arctangent: 4 atan 1 and 4 (atan 1/2 + atan 1/3) enclose pi (from Machin's formula), "
              "atan a + atan b encloses atan((a + b)/(1 - ab)) at 400 random rationals, agreement with math.atan2",
              ok and worst_f < 1e-15, "widest enclosure 2^%d, largest float difference %.1e" % (
                  worst_w.bit_length() - LP, worst_f))
    # the slope-jump identities
    jsum = sum(cert.jump)
    j1 = sum(j * x for j, x in zip(cert.jump, cert.e))
    j2 = sum(j * x * x for j, x in zip(cert.jump, cert.e)) / 2
    rep.check("self-test: the slope jumps satisfy sum j_k = 0, sum j_k e_k = 0, sum j_k e_k^2/2 = m exactly",
              jsum == 0 and j1 == 0 and j2 == cert.m, "%d nonzero jumps, sum |j_k| = %.6e" % (
                  sum(1 for x in cert.jump if x), float(sum(abs(x) for x in cert.jump))))
    # U against brute-force quadrature: random points of the support (one at a node), outside it, and complex points
    nodes, weights = _gauss_legendre(20)
    pts = []
    for _ in range(16):
        k = rng.choice(cert.support)
        pts.append(cert.e[k] + Fraction(rng.randrange(1, 1000), 1000) * (cert.e[k + 1] - cert.e[k]))
    pts += [cert.e[rng.choice(cert.support)] for _ in range(3)]
    pts += [Fraction(1, 40), Fraction(1), Fraction(6)]
    worst = 0.0
    for x in pts:
        iv = U_real_lp(x.numerator, x.denominator)
        worst = max(worst, abs((iv[0] + iv[1]) / 2 / (1 << LP) - pl_quad_U(cert, float(x), nodes, weights)))
    cpts = [(c[1], c[2]) for c in PL["cr"]]
    for _ in range(4):
        cpts.append((round(rng.uniform(0.1, 5.0) * (1 << KC)), round(rng.uniform(0.01, 0.5) * (1 << KC))))
    worstc = 0.0
    for CX, CY in cpts:
        iv = U_complex_lp(CX, CY, KC)
        z = complex(CX / 2.0 ** KC, CY / 2.0 ** KC)
        worstc = max(worstc, abs((iv[0] + iv[1]) / 2 / (1 << LP) - pl_quad_U(cert, z, nodes, weights)))
    rep.check("self-test: closed-form U agrees with brute-force Gauss-Legendre quadrature over every support cell at "
              "%d real points (three of them nodes, three outside S) and %d complex points (the %d nonreal roots in the "
              "upper half-plane and four random points)" % (len(pts), len(cpts), len(PL["cr"])),
              worst < 1e-10 and worstc < 1e-10,
              "largest differences %.1e (real), %.1e (complex)" % (worst, worstc))
    # Vieta
    vbad = [j for j, Q in enumerate(cert.A) if not vieta_pl(Q, PL["rroots"][j], alld[j])]
    rep.check("self-test: Vieta: for every Q_j the root enclosures give e1 = -a_{n-1}/a_n, e2 = a_{n-2}/a_n and "
              "prod |root|^2 = (a_0/a_n)^2 within their radii (exact rationals)", not vbad,
              "failed for %s" % vbad if vbad else "all %d polynomials" % len(cert.A))


def pl_model_tests(rep, cert, rng, count=48, pieces=160):
    """The cell model against the direct sums, and piece bounds against rigorous point values."""
    ok, worst = True, 0.0
    models = {}
    for _ in range(count):
        i = rng.choice(cert.support)
        cm = models.get(i) or models.setdefault(i, PLCell(i))
        num, k = rng.randrange(-(1 << 30), (1 << 30) + 1), 30
        a = cm.f_interval(num, k)
        b = f_direct_pl(cm.xs(num, k))
        ok &= a[0] <= b[1] and b[0] <= a[1]
        worst = max(worst, abs((a[0] + a[1]) / 2 - (b[0] + b[1]) / 2) / ONE)
    rep.check("self-test: the cell model (Taylor polynomial, remainders, near terms) and the direct sums over all "
              "nodes agree (enclosures overlap) at %d random points" % count, ok,
              "largest midpoint difference %.1e" % worst)
    bad, n, tight = 0, 0, []
    cells = list(cert.support)
    for trial in range(pieces):
        i = rng.choice(cells)
        cm = models.get(i) or models.setdefault(i, PLCell(i))
        s = rng.choice([0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40])
        mode = trial % 3
        j = 0 if mode == 0 else ((1 << s) - 1 if mode == 1 else rng.randrange(1 << s))
        lb, t = cm.piece(j, s, 1 << 400)
        na = 2 * j - (1 << s)
        pts = [(na, s), (na + 2, s)] + [((na << 5) + 2 * q, s + 5) for q in range(1, 32)]
        pts += [((na << (e_ + 1)) + 1, s + e_ + 1) for e_ in range(1, 40)]
        pts += [(((na + 2) << (e_ + 1)) - 1, s + e_ + 1) for e_ in range(1, 40)]
        best = None
        for num, k in pts:
            ub = cm.f_upper(num, k)
            if ub is None:
                continue
            n += 1
            bad += ub < lb
            best = ub if best is None else min(best, ub)
        if best is not None:
            tight.append((best - lb) / ONE)
    tight.sort()
    rep.check("self-test: piece lower bounds never exceed rigorous point values of f (%d random pieces, depths 0-40, "
              "%d points crowded towards both ends)" % (pieces, n), bad == 0,
              "%d violations; sampled minimum minus bound: median %.1e, max %.1e" % (
                  bad, tight[len(tight) // 2], tight[-1]))


# ---- red controls ----

def pl_controls(rep, cert, e5, rroots, croots, alld, workers):
    """Each must come out wrong."""
    xs, i_star = e5["x"], e5["rmin"]["cell"]
    UB = Fraction(e5["UB"], ONE)
    # 1. a weight raised just enough to break E5 at its minimiser: the one with the largest log|Q(x*)| > 0
    lq = []
    for j, Q in enumerate(cert.A):
        v = ev_frac(Q, xs)
        lq.append((math.log(abs(float(v))), j))
    lmax, jw = max(lq)
    # f_new(x*) = f(x*) - delta ln|Q(x*)|: delta = 2 UB / (a lower bound of ln|Q(x*)|) makes it <= -UB
    Lq = ln_q_lp(abs(ev_int(cert.A[jw], xs.numerator, xs.denominator)), xs.denominator ** (len(cert.A[jw]) - 1))
    lqlo = Fraction(Lq[0], 1 << LP)
    delta = Fraction(cdiv((2 * UB / lqlo).numerator * 10 ** 16, (2 * UB / lqlo).denominator), 10 ** 16)
    b2 = list(cert.b)
    b2[jw] += delta
    c2 = PLCert(cert.e, cert.h, cert.A, b2, cert.bound)
    pl_prepare(c2, rroots, croots)
    fv = f_direct_pl(xs)
    cm = PLCell(i_star)
    r = cm.bisect(0, smax=30, max_pieces=20000)
    s_ = sum(bb * (len(Q) - 1) for bb, Q in zip(b2, cert.A))
    rep.check("control: b of Q%d (the largest log|Q(x*)| = %.4f at the minimiser) raised by %s (relative %.2e, twice "
              "what breaks E5 there): E5 refuted at x*, and the bisection does not certify f >= 0 on cell %d" % (
                  jw, lmax, sci_ceil(delta, 3), float(delta / cert.b[jw]), i_star), fv[1] < 0 and bool(r["fails"]),
              "f(x*) <= %s (rigorous upper end); %d pieces, %d uncertified; E3 still holds (sum/m = %s)" % (
                  sci_ceil(Fraction(fv[1], ONE), 4), r["pieces"], len(r["fails"]), dec(s_ / cert.m, 8)))
    # 2. the highest node lowered by 1%: E4 must fail for some Q (and the stated bound is no longer T/m)
    k = max(range(len(cert.h)), key=lambda k: cert.h[k])
    h2 = list(cert.h)
    h2[k] = h2[k] * Fraction(99, 100)
    c3 = PLCert(cert.e, h2, cert.A, cert.b, cert.bound)
    pl_prepare(c3, rroots, croots)
    t0 = time.time()
    vals = pl_E4_all(workers)
    neg = [j for j, iv in vals if iv[1] < 0]
    worst = min(vals, key=lambda v: v[1][1])
    rep.check("control: the highest node (e = %s, height %.4f) lowered by 1%% (mass and trace recomputed): E4 refuted" % (
        cert.e[k], float(cert.h[k])), bool(neg) and c3.T / c3.m != cert.bound,
        "int log|Q| d mu < 0 rigorously for %d of %d polynomials, most negative Q%d in %s; E6 refuses the old bound; "
        "%.0f s" % (len(neg), len(vals), worst[0], enc(worst[1], 6), time.time() - t0))
    # 3. E3: weights scaled by 1.3; 4. E6: a bound 1e-15 low; 5. E1: a negative height, edges out of order
    s = sum(bb * Fraction(13, 10) * (len(Q) - 1) for bb, Q in zip(cert.b, cert.A))
    rep.check("control: weights scaled by 1.3 violate E3", not s < cert.m, "sum/m = %s" % dec(s / cert.m, 6))
    b = cert.T / cert.m
    low = b - Fraction(1, 10 ** 15)
    rep.check("control: a bound 1e-15 below T/m is refused by E6", low != b and low < b)
    kk = cert.support[len(cert.support) // 2]
    h3 = list(cert.h)
    h3[kk] = -Fraction(1, 10 ** 9)
    e3 = list(cert.e)
    e3[kk], e3[kk + 1] = e3[kk + 1], e3[kk]
    rep.check("control: a height -1e-9 and two swapped nodes are refused by E1",
              not all(x >= 0 for x in h3) and not all(x < y for x, y in zip(e3, e3[1:])))
    # 6. E2: reducible, non-primitive, repeated, negative leading coefficient
    quiet = Report.__new__(Report)
    quiet.fails = []
    for label, Q in (("x^2 - 3x + 2 = (x - 1)(x - 2)", [2, -3, 1]), ("2x^2 - 6x + 2 (not primitive)", [2, -6, 2]),
                     ("Q1 again, negated", [-c for c in cert.A[1]]), ("-x + 3 (negative leading coefficient)", [3, -1]),
                     ("(x^2 + 1)(x^2 - 3x + 1), no real factor", [1, -3, 2, -3, 1])):
        bad = p298_E2(quiet, list(cert.A) + [Q], list(cert.b) + [Fraction(1, 10 ** 6)], quiet=True)
        rep.check("control: E2 refuses %s" % label, any(j == len(cert.A) for j, _ in bad),
                  "; ".join(d for j, d in bad if j == len(cert.A))[:150])
    # 7. a nonreal root moved by 1e-12 (radius kept) is caught by the Vieta test
    pl_prepare(cert, rroots, croots)
    jc = max(range(len(cert.A)), key=lambda j: len(cert.A[j]) if croots[j] else 0)
    disc = alld[jc]
    k_ = next(i for i, d_ in enumerate(disc) if d_[1] - d_[2] > 0)
    moved = list(disc)
    moved[k_] = (disc[k_][0] + Fraction(1, 10 ** 12), disc[k_][1], disc[k_][2])
    rep.check("control: a nonreal root of Q%d moved by 1e-12 (radius kept) fails the Vieta test" % jc,
              vieta_pl(cert.A[jc], rroots[jc], disc) and not vieta_pl(cert.A[jc], rroots[jc], moved))
    # 8. the bisection cannot certify a goal above a rigorous upper bound of f
    cm = PLCell(i_star)
    r = cm.bisect(e5["UB"] + qiv(Fraction(1, 10 ** 12))[1], smax=24, max_pieces=20000)
    rep.check("control: goal min_S f + 1e-12 is not certified in cell %d (bisection to depth 24)" % i_star,
              bool(r["fails"]), "%d pieces, %d uncertified" % (r["pieces"], len(r["fails"])))


# ----------------------------------------------------------------------------------------------------------------
# Rows
# ----------------------------------------------------------------------------------------------------------------

PIN11 = [[0, 1], [1, -1], [1, -3, 1], [-1, 6, -5, 1], [1, -8, 14, -7, 1], [1, -7, 13, -7, 1], [-1, 15, -35, 28, -9, 1],
         [-1, 11, -29, 26, -9, 1], [-1, 13, -32, 27, -9, 1], [-1, 9, -6, 1], [-1, 12, -31, 27, -9, 1]]
PIN_PAGE_B = ["0.038759530583", "0.079868388122", "0.04631766723", "0.053369420893", "0.019753937664",
              "0.036545507753", "0.016049367094", "0.009731525491", "0.007106935573", "0.018100957401",
              "0.002893658474"]
KOS8 = [[-2, 1], [1, -4, 1], [2, -4, 1], [-3, 9, -6, 1], [1, -15, 53, -73, 43, -11, 1], [1, -14, 51, -72, 43, -11, 1],
        [-1, 18, -89, 172, -150, 64, -13, 1], [-1, 16, -78, 157, -143, 63, -13, 1]]
KOS_README_B = {  # the exact weights in the [Kos2026] README table
    (-1, 1): "79216914042837/1000000000000000", (0, 1): "494739586843/8000000000000",
    (-2, 1): "4041043491479/100000000000000", (1, -3, 1): "17462027687139/500000000000000",
    (1, -4, 1): "4580571915059/250000000000000", (2, -4, 1): "10103990454883/1000000000000000",
    (-1, 6, -5, 1): "17637921744613/500000000000000", (-1, 9, -6, 1): "110660999533/200000000000000",
    (-3, 9, -6, 1): "69448588053/500000000000000", (1, -7, 13, -7, 1): "17755144029567/1000000000000000",
    (1, -8, 14, -7, 1): "17418022936419/1000000000000000", (-1, 12, -31, 27, -9, 1): "12331565309643/1000000000000000",
    (-1, 15, -35, 28, -9, 1): "10368009092229/1000000000000000", (-1, 11, -29, 26, -9, 1): "4064736833123/500000000000000",
    (-1, 13, -32, 27, -9, 1): "724976585669/125000000000000",
    (1, -15, 53, -73, 43, -11, 1): "3614401969/1000000000000000",
    (1, -14, 51, -72, 43, -11, 1): "73408487/1000000000000000",
    (-1, 18, -89, 172, -150, 64, -13, 1): "472219232001/62500000000000",
    (-1, 16, -78, 157, -143, 63, -13, 1): "1521492933501/250000000000000"}

ROWS = {
    "kos": {
        "title": "ROW kos: [Kos2026] @90c57242, C86 <= 1.8084990*",
        "file": "kos/c86_A19d7_robust_1.8084990.json",
        "hashes": [("[Kos2026] SHA256SUMS", "1a9e108e781da18ae04b04303c96bd01285f81c004abc2aea966ea104697845c", ""),
                   ("PR #292 thread", "1a9e108e", "45c")],
        "cells": 633, "support": ("11/200", "519/100"),
        "m": "199999999999999703/200000000000000000",
        "bound": "144679917840593158601/79999999999999881200",
        "printed": [("86a page", "1.8084990")],
        "digits": [("86a page", "1.80849897300"), ("[Kos2026] README and FORMAT.md", "1.80849897300741")],
        "ratio": "0.83542", "margin": Fraction(1, 10 ** 7),
        "polys": PIN11 + KOS8, "readme_b": KOS_README_B, "search": 7,
        "e4_stated": [("86a page", Fraction(199, 10 ** 8))],
        "e4_rounded": ("PR #292, to 5 significant digits", "2.0000e-6"),
    },
    "pin": {
        "title": "ROW pin: [Pin26] @116a58f, C86 <= 1.8117251*",
        "file": "pin/cert_n1200_A11.json",
        "hashes": [("[Pin26] README", "722a287b74a9f84f8c6c2f505e7e6bb987ca634710b031287a5b01619bffe53d", ""),
                   ("PR #242 thread", "722a287b74a9", "")],
        "cells": 629, "support": ("1/20", "1037/200"),
        "m": "39999999999997/40000000000000",
        "bound": "144938007612277893/79999999999994000",
        "printed": [("86a page", "1.8117251")],
        "digits": [("86a page", "1.81172509515"), ("[Pin26] README", "1.81172509515360"),
                   ("PR #242", "1.81172509515360954")],
        "ratio": "0.82978", "margin": Fraction(0),
        "polys": PIN11, "page_b": PIN_PAGE_B,
        "e4_stated": [],
        "e4_about": ("86a page and [Pin26] README: each about 1e-4", Fraction(1, 10 ** 4)),
    },
    "vid": {
        "title": "ROW vid: [Pat2026] (PR #270) @fb80010 (= @7c6586f for this file), C86 <= 1.8077686*",
        "file": "vid/c86_deg8_h400_1.807768536863.json",
        "hashes": [("PR #270, 86a row", "7ca02d70d735", "")],
        "ncells": 2054, "width": "1/400",
        "bound": "578485931795944337659/319999999999999303200",
        "printed": [("PR #270, 86a row", "1.8077686"), ("file name and the row's VALID line", "1.807768536863")],
        "digits": [("PR #270, 86a row", "1.80776853686"), ("[Pat2026] README", "1.80776853686")],
        "margin": Fraction(0), "npos": 19, "deg8": 2, "search": 8,
        "e4_stated": [],
    },
}


def run_row(key, root, do_controls, workers):
    spec = ROWS[key]
    rep = Report(spec["title"])
    rng = random.Random(86)
    path = os.path.join(root, spec["file"])
    h = sha256(path)
    rep.info("sha256 " + spec["file"], h)
    rep.info("86a page", "states the commit only (no hash) for this row")
    for where, pre, suf in spec["hashes"]:
        rep.check("sha256 matches %s (%s%s)" % (where, pre, ("..." + suf) if suf else ("" if len(pre) == 64 else "...")),
                  h.startswith(pre) and h.endswith(suf))
    raw = load_cert(path)
    rep.info("keys read", "edges, densities, A, b, bound (all exact rationals or integers); other keys ignored: %s" %
             ", ".join(k for k in raw["raw"] if k not in ("edges", "densities", "A", "b", "bound")))
    cert = Cert(raw)
    # ---- E1
    lo, hi = check_E1(rep, cert)
    if "cells" in spec:
        rep.check("support and cell count as on the 86a page (%d cells in [%s, %s])" % (
            (spec["cells"],) + spec["support"]),
            len(cert.support) == spec["cells"] and (lo, hi) == tuple(Fraction(s) for s in spec["support"]))
    if "ncells" in spec:
        rep.check("%d cells of width %s, as the PR #270 row says" % (spec["ncells"], spec["width"]),
                  cert.n == spec["ncells"] and cert.uniform() == (True, Fraction(spec["width"])))
    if "m" in spec:
        rep.check("mass m = %s as on the 86a page" % spec["m"], cert.m == Fraction(spec["m"]))
    # ---- E2
    roots = check_E2(rep, cert)
    if roots is None:
        return rep.done()
    keyset = {tuple(Q) for Q in cert.A} | {tuple(-c for c in Q) for Q in cert.A}
    if "polys" in spec:
        rep.check("the %d polynomials are those listed on the 86a page (up to sign)" % len(spec["polys"]),
                  len(cert.A) == len(spec["polys"]) and all(tuple(Q) in keyset for Q in spec["polys"]))
    if "npos" in spec:
        npos = sum(1 for x in cert.b if x > 0)
        n8 = sum(1 for Q, x in zip(cert.A, cert.b) if x > 0 and len(Q) - 1 == 8)
        rep.check("%d polynomials with positive weight, %d of degree 8, as the row says" % (spec["npos"], spec["deg8"]),
                  npos == spec["npos"] and n8 == spec["deg8"], "%d and %d" % (npos, n8))
    if "page_b" in spec:
        okb = all(Fraction(bs) == cert.b[cert.A.index(Q)] for Q, bs in zip(spec["polys"], spec["page_b"]))
        rep.check("the weights b_Q printed on the 86a page are the certificate's exact weights", okb)
    if "readme_b" in spec:
        okb = all(Fraction(spec["readme_b"][tuple(Q)]) == cert.b[j] for j, Q in enumerate(cert.A))
        rep.check("the weights equal the exact table of the [Kos2026] README", okb)
    if "search" in spec:
        info = []
        for Q, rs in zip(cert.A, roots):
            tr = -Q[-2] if len(Q) > 1 else 0
            info.append(all(0 <= a and b <= Fraction(53, 10) for a, b in rs) and Q[-1] == 1 and
                        Fraction(tr, len(Q) - 1) <= 2 and len(Q) - 1 <= spec["search"])
        rep.check("each Q is monic of degree <= %d, all roots in [0, 5.3], trace <= 2 deg (the stated search space)"
                  % spec["search"], all(info))
    # ---- E3
    ratio = check_E3(rep, cert)
    if "ratio" in spec:
        rounded = dec(ratio + Fraction(5, 10 ** 6), 5)
        rep.check("sum b deg / m rounds to %s, as on the 86a page" % spec["ratio"], rounded == spec["ratio"])
    # ---- E6
    b = check_E6(rep, cert, spec["printed"])
    rep.check("T/m = %s as stated for this row" % spec["bound"], b == Fraction(spec["bound"]))
    for where, s in spec["digits"]:
        n = len(s.split(".")[1])
        rep.check("T/m begins %s... (%s)" % (s, where), dec(b, n) == s)
    # ---- E4
    prepare(cert, roots)
    selftests(rep, cert, rng)
    e4 = check_E4(rep, cert, roots)
    smallest = min(iv[0] for _, iv in e4)
    for where, val in spec["e4_stated"]:
        rep.check("smallest E4 >= %s (%s)" % (sci_floor(val, 3), where), Fraction(smallest, ONE) >= val)
    if "e4_rounded" in spec:
        where, s = spec["e4_rounded"]
        lo_, hi_ = min((iv for _, iv in e4), key=lambda v: v[0])
        val = Fraction(s.split("e")[0]) * Fraction(10) ** int(s.split("e")[1])
        half = Fraction(5, 10 ** (len(s.split("e")[0].split(".")[1]) + 1)) * Fraction(10) ** int(s.split("e")[1])
        rep.check("smallest E4 rounds to %s (%s)" % (s, where),
                  val - half <= Fraction(lo_, ONE) and Fraction(hi_, ONE) < val + half,
                  "exact enclosure %s" % enc((lo_, hi_), 12))
    if "e4_about" in spec:
        where, val = spec["e4_about"]
        dev = max(abs(Fraction(iv[0], ONE) / val - 1) for _, iv in e4)
        rep.check("every E4 within 1%% of %s (%s)" % (sci_floor(val, 1), where), dev < Fraction(1, 100),
                  "largest relative deviation %.2e" % float(dev))
    # ---- E5
    model_vs_direct(rep, cert, rng)
    piece_soundness(rep, cert, rng)
    LB, UB, xmin, res, failed = check_E5(rep, cert, workers, spec["margin"])
    if spec["margin"] > 0:
        rep.check("E5 margin claimed on the 86a page: f >= %s on S" % sci_floor(spec["margin"], 2),
                  not failed and Fraction(LB, ONE) >= spec["margin"],
                  "certified f >= %s" % sci_floor(Fraction(LB, ONE), 6))
    rep.check("E5: f >= 0 on S", not failed and LB >= 0)
    # ---- consequences, checked separately
    I = energy(cert)
    sb = (0, 0)
    for j, iv in e4:
        sb = ivadd(sb, ivq(iv, cert.b[j].numerator, cert.b[j].denominator))
    intf = (I[0] - sb[1], I[1] - sb[0])
    rep.info("energy I(mu) = int int log|x - y| dmu dmu", enc(I, 12))
    rep.check("consistency: int f dmu = I(mu) - sum_j b_j E4_j >= m min_S f (E4, E5 and the energy are computed "
              "by different code)", intf[0] >= frq(cert.m, LB, 1),
              "int f dmu in %s, m x (certified min f) = %s" % (enc(intf, 8), sci_floor(cert.m * Fraction(LB, ONE), 8)))
    rep.check("energy of nu = mu/m is positive, as Smith's criterion needs (a consequence of E4, E5)", I[0] > 0,
              "I(nu) = I(mu)/m^2")
    worst = None
    pts = gap_points(cert)
    for x in pts:
        lo_, hi_ = f_direct(cert, x)
        if worst is None or lo_ < worst[1]:
            worst = (x, lo_)
    rep.check("minimum principle, sampled: f > 0 at the %d gap midpoints and 3 exterior points (not needed; a "
              "consistency check)" % (len(pts) - 3), worst[1] > 0,
              "smallest f >= %s at x = %s" % (sci_floor(Fraction(worst[1], ONE), 6), float(worst[0])))
    if do_controls:
        scan_pts = [cert.edges[r["cell"]] + Fraction(r["rho_min"]) * SHD["H"] for r in res]
        controls(rep, cert, roots, scan_pts, res, LB, UB)
    rep.info("what rests on the cited theorems",
             "f >= 0 on S gives f >= 0 on C (minimum principle); E2-E5 make mu/m arithmetic ([Smi24] Prop. 5.15); "
             "Smith's Thm 1.5 on [min S, max S] (capacity %s/4 > 1) then gives C86 <= T/m" % (hi - lo))
    return rep.done()


def run_p298(root, do_controls, workers):
    """PR #298 (rpinho/math-certificates @0d82162, c86/adaptive-cone): a piecewise-linear density on 5,375
    rational nodes, 28 auxiliaries, ten with nonreal roots. Every condition, E1 to E6."""
    rep = Report("ROW p298: PR #298 @0d82162 (c86/adaptive-cone), C86 <= 1.80560940*")
    rng = random.Random(298)
    rel = "p298/certificate.json"
    path = os.path.join(root, rel)
    h = sha256(path)
    rep.info("sha256 " + rel, h)
    full = "10d68ec249f47c0504182d419fd8ddd1ca1c917e42a4bba59d944b9d1de0dad5"
    rep.check("sha256 matches the README and SHA256SUMS of c86/adaptive-cone (%s)" % full, h == full)
    d = json.load(open(path))
    e = [rat(x) for x in d["edges"]]
    hgt = [rat(x) for x in d["node_densities"]]
    A = [list(Q) for Q in d["A"]]
    assert all(isinstance(c, int) and not isinstance(c, bool) for Q in A for c in Q)
    b = [rat(x) for x in d["b"]]
    bound = rat(d["bound"])
    rep.info("keys read", "edges, node_densities, A, b, bound; ignored: %s" % ", ".join(
        k for k in d if k not in ("edges", "node_densities", "A", "b", "bound")))
    # ---- E1
    rep.check("E1 %d nodes, strictly increasing, %d heights, %d polynomials and weights" % (len(e), len(hgt), len(A)),
              len(e) == len(hgt) and len(A) == len(b) and all(x < y for x, y in zip(e, e[1:])))
    rep.check("E1 every node height >= 0 and every weight >= 0", all(x >= 0 for x in hgt) and all(x >= 0 for x in b),
              "%d positive heights, %d positive weights" % (sum(1 for x in hgt if x > 0), sum(1 for x in b if x > 0)))
    cert = PLCert(e, hgt, A, b, bound)
    m, T = cert.m, cert.T
    lo, hi = e[cert.support[0]], e[cert.support[-1] + 1]
    rep.check("E1 m > 0 and min S > 0 (density zero outside the nodes)", m > 0 and lo > 0 and hgt[0] == 0 and
              hgt[-1] == 0, "m = %s, S in [%s, %s]" % (dec(m, 18), lo, hi))
    gaps = sum(1 for a_, c_ in zip(cert.support, cert.support[1:]) if c_ != a_ + 1)
    widths = sorted(set(e[k + 1] - e[k] for k in cert.support))
    rep.info("support", "%d cells with a positive end height, widths %s to %s, in [%s, %s] (%d gaps)" % (
        len(cert.support), widths[0], widths[-1], lo, hi, gaps))
    # ---- E2
    bad = p298_E2(rep, A, b)
    if bad:
        return rep.done()
    # ---- E3
    s = sum(bb * (len(Q) - 1) for bb, Q in zip(b, A))
    rep.check("E3 sum_j b_j deg Q_j < m", s < m, "sum/m = %s" % dec(s / m, 10))
    # ---- E6
    rep.check("E6 T/m equals the certificate's bound exactly", T / m == bound, "T/m = %s = %s..." % (T / m,
                                                                                                    dec(T / m, 20)))
    rep.check("E6 the README's 288897505229428473/160000001279369600", bound == Fraction(288897505229428473,
                                                                                         160000001279369600))
    rep.check("E6 printed upper bound 1.80560940 (PR #298 title) >= T/m", Fraction("1.80560940") >= T / m,
              "difference %.3e" % float(Fraction("1.80560940") - T / m))
    # ---- roots, and the data for E4 and E5
    t0 = time.time()
    rroots, croots, alld, nreals = p298_roots(rep, cert)
    pl_prepare(cert, rroots, croots)
    rep.info("prepared", "%d nodes with a nonzero slope jump, D = %d (D e_k integers), lcm of the jump denominators "
             "2^%d; roots %.1f s" % (len(PL["NE"]), PL["D"], PL["DJ"].bit_length(), time.time() - t0))
    pl_selftests(rep, cert, rng, alld)
    # ---- E4
    t0 = time.time()
    span = max(abs(c - x) + rad for _, c, rad, _ in PL["rr"] for x in (cert.e[min(cert.support)],
                                                                       cert.e[max(cert.support) + 1]))
    rep.check("E4 every real root within distance 6 of every node with a nonzero jump (the Lipschitz bound used)",
              span <= 6, "largest distance %.4f" % float(span))
    e4 = pl_E4_all(workers)
    for j, iv in e4:
        rep.info("E4 int log|Q%d| dmu, Q%d = %s" % (j, j, pstr(A[j])), enc(iv, 12))
    worst = min(e4, key=lambda v: v[1][0])
    rep.check("E4 int log|Q_j| dmu > 0 for every j (rigorous lower ends)", all(iv[0] > 0 for _, iv in e4),
              "smallest: Q%d = %s, in %s; %.1f s" % (worst[0], pstr(A[worst[0]]), enc(worst[1], 12), time.time() - t0))
    t0 = time.time()
    quad = dict(pl_E4_all(workers, _pl_E4_quad))
    dq = max(abs((iv[0] + iv[1]) / 2 / ONE - quad[j]) for j, iv in e4)
    rep.check("self-test: every E4 value agrees with brute-force quadrature of U at its %d roots (floats)" % sum(
        len(Q) - 1 for Q in A), dq < 1e-11, "largest difference %.1e; %.1f s" % (dq, time.time() - t0))
    # ---- E5
    inS = sum(1 for _, c, rad, _ in PL["rr"] for k in cert.support if e[k] <= c - rad and c + rad <= e[k + 1])
    rep.info("E5 roots on the support", "%d of the %d real roots lie in S, where f tends to +infinity; a piece that "
             "may hold one takes -b log(its largest distance to the root), the others the tangent" % (inS, len(PL["rr"])))
    pl_model_tests(rep, cert, rng)
    e5 = pl_check_E5(rep, cert, workers)
    if e5 is None:
        return rep.done()
    # ---- consequences, checked separately
    pts = []
    sup = cert.support
    for a_, c_ in zip(sup, sup[1:]):
        if c_ != a_ + 1:
            pts.append((e[a_ + 1] + e[c_]) / 2)
    pts += [lo / 2, hi + Fraction(1, 5), Fraction(10)]
    worstg = None
    for x in pts:
        lo_, hi_ = f_direct_pl(x)
        if worstg is None or lo_ < worstg[1]:
            worstg = (x, lo_)
    rep.check("minimum principle, sampled: f > 0 at the %d gap midpoints and 3 exterior points (not needed; a "
              "consistency check)" % (len(pts) - 3), worstg[1] > 0,
              "smallest f >= %s at x = %s" % (sci_floor(Fraction(worstg[1], ONE), 6), float(worstg[0])))
    if do_controls:
        pl_controls(rep, cert, e5, rroots, croots, alld, workers)
    rep.info("what rests on the cited theorems",
             "f >= 0 on S gives f >= 0 on C (minimum principle: f is harmonic off S and the roots, tends to +inf at "
             "the roots and at infinity by E3); E2-E5 make mu/m arithmetic ([Smi24] Prop. 5.15); Smith's Thm 1.5 on "
             "an interval containing S of capacity > 1 ([1/1000, 13], as PROOF.md says) then gives C86 <= T/m")
    tms = os.times()
    rep.info("CPU time", "%.0f s (this process and its workers)" % (tms[0] + tms[1] + tms[2] + tms[3]))
    ok = rep.done()
    print("  ROW STATUS: FULL (E1-E6)")
    return ok


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in list(ROWS) + ["p298"]:
        print(__doc__)
        sys.exit(2)
    workers = os.cpu_count() or 1
    for a in sys.argv:
        if a.startswith("--workers="):
            workers = int(a.split("=")[1])
    if sys.argv[1] == "p298":
        ok = run_p298(sys.argv[2], "--controls" in sys.argv, workers)
    else:
        ok = run_row(sys.argv[1], sys.argv[2], "--controls" in sys.argv, workers)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
