#!/usr/bin/env python3
"""Exact replay of the 76a lower-bound row [BS2007] (teorth/optimizationproblems, main 87b4ebb):
l_0(d) >= 3d^2, "achieved by the classical family phi(x,y) - psi(z,t) = 0", hence C_76 = limsup l_0(d)/d^2 >= 3.

We take the member phi = x^d + y^d, psi = -(z^d + t^d): the Fermat surface F_d: x^d + y^d + z^d + t^d = 0 in P^3, and
check for d = 3..12, exactly in the cyclotomic field Q(w), w a primitive 2d-th root of unity (elements are rational
polynomials modulo the cyclotomic polynomial Phi_2d, which we compute):
  1. the 3d^2 lines {x = a y, z = b t}, {x = a z, y = b t}, {x = a t, y = b z}, a, b running over the d roots of
     u^d = -1 (the odd powers of w): each line, spanned by two explicit points P, Q, lies on F_d, checked by
     F_d(P + kQ) = 0 for k = 0..d and F_d(Q) = 0 (d + 2 points; a form of degree d vanishing at d + 1 points of a line
     vanishes on it);
  2. they are pairwise distinct: Plucker coordinates of (P, Q), normalised by the first nonzero coordinate (inverse in
     the field by the extended Euclidean algorithm), give 3d^2 different vectors;
  3. F_d is smooth: its partial derivatives d x^(d-1), ..., d t^(d-1) vanish together only at 0 (characteristic 0).
So l_0(d) >= 3d^2 for every d >= 3 (the lines are explicit for every d; we check d <= 12), and C_76 >= 3. That the count
in the family is exactly 3d^2 except in degrees 4, 6, 8, 12, 20 is read, not checked (not needed for the bound).
Red controls (--controls): a line with a^d = +1 is not on F_d; one line listed twice is caught by the distinctness test
(3d^2 - 1 distinct among 3d^2); a claim of 3d^2 + 1 lines in this family is refused.
usage: python3 -I fermat_lines_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


# ---------------------------------------------------------------- polynomials over Q (lists, low degree first)
def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(a, b):
    if not a or not b:
        return []
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return trim(r)


def pdivmod(a, b):
    a, b = trim(a), trim(b)
    q = [F(0)] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]
        k = len(a) - len(b)
        q[k] = c
        for i, y in enumerate(b):
            a[i + k] -= c * y
        a = trim(a)
    return trim(q), a


def cyclotomic(m):
    p = [F(-1)] + [F(0)] * (m - 1) + [F(1)]          # x^m - 1
    for k in range(1, m):
        if m % k == 0:
            p, r = pdivmod(p, cyclotomic(k))
            assert not r
    return p


class Field:
    def __init__(self, m):
        self.phi = cyclotomic(m)
        self.deg = len(self.phi) - 1

    def red(self, a):
        return tuple(trim(pdivmod(a, self.phi)[1]))

    def mul(self, a, b):
        return self.red(pmul(list(a), list(b)))

    def add(self, a, b):
        n = max(len(a), len(b))
        return tuple(trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]))

    def neg(self, a):
        return tuple(-x for x in a)

    def power_of_w(self, k):
        return self.red([F(0)] * k + [F(1)])

    def pw(self, a, e):
        r = (F(1),)
        for _ in range(e):
            r = self.mul(r, a)
        return r

    def inv(self, a):
        """Inverse by the extended Euclidean algorithm (invariant s_i a = r_i mod Phi), checked by multiplication."""
        r0, r1 = list(self.phi), list(a)
        s0, s1 = [], [F(1)]
        while trim(r1):
            q, r = pdivmod(r0, r1)
            r0, r1 = r1, r
            qs = pmul(q, s1)
            n = max(len(s0), len(qs))
            s0, s1 = s1, trim([(s0[i] if i < len(s0) else 0) - (qs[i] if i < len(qs) else 0) for i in range(n)])
        g = trim(r0)
        assert len(g) == 1, "not invertible"
        res = self.red([x / g[0] for x in s0])
        assert self.mul(res, a) == (F(1),), "inverse check failed"
        return res


def lines(K, d, a_exp_parity=1, duplicate=False):
    """All 3 d^2 lines as pairs of points (P, Q); a, b = w^(2i+parity)."""
    one, zero = (F(1),), ()
    roots = [K.power_of_w(2 * i + a_exp_parity) for i in range(d)]
    L = []
    for a in roots:
        for b in roots:
            L.append(((a, one, zero, zero), (zero, zero, b, one)))     # x = a y, z = b t
            L.append(((a, zero, one, zero), (zero, b, zero, one)))     # x = a z, y = b t
            L.append(((a, zero, zero, one), (zero, b, one, zero)))     # x = a t, y = b z
    if duplicate:
        L[-1] = L[0]
    return L


def fermat(K, d, v):
    s = ()
    for c in v:
        s = K.add(s, K.pw(c, d))
    return s


def on_surface(K, d, P, Q):
    if fermat(K, d, Q):
        return False
    for k in range(d + 1):
        pt = tuple(K.add(p, K.mul((F(k),), q)) for p, q in zip(P, Q))
        if fermat(K, d, pt):
            return False
    return True


def plucker(K, P, Q):
    co = []
    for i in range(4):
        for j in range(i + 1, 4):
            co.append(K.add(K.mul(P[i], Q[j]), K.neg(K.mul(P[j], Q[i]))))
    k = next(i for i, c in enumerate(co) if c)
    iv = K.inv(co[k])
    return tuple(K.mul(c, iv) for c in co)


def check_degree(d, parity=1, duplicate=False):
    K = Field(2 * d)
    w = K.power_of_w(1)
    prim = K.pw(w, d) == (F(-1),)                     # w^d = -1
    L = lines(K, d, parity, duplicate)
    on = all(on_surface(K, d, P, Q) for P, Q in L)
    distinct = len({plucker(K, P, Q) for P, Q in L})
    return prim, len(L), on, distinct


def main():
    ok = True
    print("76a row [BS2007]: l_0(d) >= 3d^2 from the Fermat surfaces x^d + y^d + z^d + t^d = 0 (smooth), so C_76 >= 3")
    for d in range(3, 13):
        prim, n, on, distinct = check_degree(d)
        good = prim and n == 3 * d * d and on and distinct == 3 * d * d
        print("  d = %2d: Q(w), w^%d = -1, [Q(w):Q] = %d; lines listed %d = 3d^2; all on F_d: %s; distinct: %d; %s"
              % (d, d, len(cyclotomic(2 * d)) - 1, n, on, distinct, "ok" if good else "FAIL"))
        ok &= good
    print("  smoothness: the partials d x^(d-1), d y^(d-1), d z^(d-1), d t^(d-1) vanish together only at 0")
    print("RESULT [BS2007]: for d = 3..12 the Fermat surface carries 3d^2 distinct lines, so l_0(d) >= 3d^2 and "
          "C_76 >= 3: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        K = Field(10)
        one, zero = (F(1),), ()
        a = K.power_of_w(0)                       # a^5 = +1
        b = K.power_of_w(1)
        bad = on_surface(K, 5, (a, one, zero, zero), (zero, zero, b, one))
        print("  d = 5, the line x = y, z = w t (1^5 = +1): %s" % ("on F_5 (UNEXPECTED)" if bad else "not on F_5, as expected"))
        c_ok &= not bad
        _, n, _, distinct = check_degree(4, duplicate=True)
        print("  d = 4 with one line listed twice: %d distinct among %d: %s" % (distinct, n,
              "caught, as expected" if distinct == n - 1 else "NOT caught (UNEXPECTED)"))
        c_ok &= distinct == n - 1
        _, n, _, distinct = check_degree(6)
        claim = distinct >= 3 * 36 + 1
        print("  a claim of 3d^2 + 1 = 109 lines in this family at d = 6: %s" % ("holds (UNEXPECTED)" if claim else
                                                                             "refused, as expected (%d)" % distinct))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
