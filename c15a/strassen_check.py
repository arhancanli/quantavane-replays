#!/usr/bin/env python3
"""Exact replay of the 15a upper-bound row [S1969] (teorth/optimizationproblems, main 87b4ebb):
omega <= log_2 7 ~ 2.8074, "recursively uses an algorithm for 2 x 2 matrices with 7 multiplication operations".

Strassen's paper (Numer. Math. 13 (1969), 354-356) was not accessible to us, so the seven products below are taken in
the form in which Strassen's algorithm is commonly stated (M1..M7 and the four output combinations), not from the
paper. What is checked, exactly:
  1. the bilinear identity: for each of the four entries of C = AB, the combination of M1..M7 has the same coefficient
     tensor (indexed by an entry of A on the left and an entry of B on the right) as the true entry; the order of the
     factors is kept, so the identity holds over every ring, commutative or not, which is what the recursion needs;
  2. the same identity evaluated on 2 x 2 block matrices whose blocks are random 3 x 3 integer matrices (a
     non-commutative ring), exactly in integers (fixed seed);
  3. the recursion: Strassen's algorithm applied recursively to 2^k x 2^k integer matrices, k = 1..6, equals the
     ordinary product, and uses exactly 7^k scalar multiplications (counted), so omega <= log_2 7;
  4. log_2 7 = 2.80735492205760410744... (enclosed exactly: 2^2.8073 < 7 < 2^2.8074 by integer powers); the printed
     2.8074 stands next to the exact expression log_2 7.
Red controls (--controls): one sign changed in M7, one product dropped from C11, and C12 = M3 - M5 are each refused by
the identity; a claimed omega <= 2.8073 is not given by 7 products (2^2.8073 < 7).
usage: python3 -I strassen_check.py [--controls]      (Python 3.9+, standard library only)
"""
import random
import sys
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING

CONTROLS = "--controls" in sys.argv
A_ENT = ["a11", "a12", "a21", "a22"]
B_ENT = ["b11", "b12", "b21", "b22"]


def lin(**kw):
    """A linear form as a dict entry -> integer coefficient."""
    return {k: v for k, v in kw.items() if v}


def prod(la, lb):
    """Coefficient tensor of (linear form in A) * (linear form in B), order kept: dict (a, b) -> coefficient."""
    t = {}
    for a, x in la.items():
        for b, y in lb.items():
            t[(a, b)] = t.get((a, b), 0) + x * y
    return {k: v for k, v in t.items() if v}


def comb(terms):
    """Sum of c * tensor over (c, tensor) pairs."""
    t = {}
    for c, m in terms:
        for k, v in m.items():
            t[k] = t.get(k, 0) + c * v
    return {k: v for k, v in t.items() if v}


def strassen_products(m7_sign=1):
    M = {}
    M[1] = prod(lin(a11=1, a22=1), lin(b11=1, b22=1))
    M[2] = prod(lin(a21=1, a22=1), lin(b11=1))
    M[3] = prod(lin(a11=1), lin(b12=1, b22=-1))
    M[4] = prod(lin(a22=1), lin(b21=1, b11=-1))
    M[5] = prod(lin(a11=1, a12=1), lin(b22=1))
    M[6] = prod(lin(a21=1, a11=-1), lin(b11=1, b12=1))
    M[7] = prod(lin(a12=1, a22=-m7_sign), lin(b21=1, b22=1))
    return M


def outputs(M, variant=None):
    C11 = [(1, M[1]), (1, M[4]), (-1, M[5]), (1, M[7])]
    if variant == "drop":
        C11 = C11[:-1]
    C12 = [(1, M[3]), (1, M[5])] if variant != "c12minus" else [(1, M[3]), (-1, M[5])]
    C21 = [(1, M[2]), (1, M[4])]
    C22 = [(1, M[1]), (-1, M[2]), (1, M[3]), (1, M[6])]
    return {"c11": comb(C11), "c12": comb(C12), "c21": comb(C21), "c22": comb(C22)}


TRUE = {"c11": {("a11", "b11"): 1, ("a12", "b21"): 1}, "c12": {("a11", "b12"): 1, ("a12", "b22"): 1},
        "c21": {("a21", "b11"): 1, ("a22", "b21"): 1}, "c22": {("a21", "b12"): 1, ("a22", "b22"): 1}}


def identity_holds(m7_sign=1, variant=None):
    out = outputs(strassen_products(m7_sign), variant)
    return all(out[k] == TRUE[k] for k in TRUE), out


# ---------------------------------------------------------------- integer matrices
def madd(X, Y, s=1):
    return [[X[i][j] + s * Y[i][j] for j in range(len(X))] for i in range(len(X))]


def mmul_naive(X, Y, counter=None):
    n = len(X)
    if counter is not None:
        counter[0] += n * n * n
    return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def split(X):
    h = len(X) // 2
    return ([r[:h] for r in X[:h]], [r[h:] for r in X[:h]], [r[:h] for r in X[h:]], [r[h:] for r in X[h:]])


def join(C11, C12, C21, C22):
    return [a + b for a, b in zip(C11, C12)] + [a + b for a, b in zip(C21, C22)]


def strassen(X, Y, counter, base=1):
    n = len(X)
    if n == base:
        return mmul_naive(X, Y, counter)
    a11, a12, a21, a22 = split(X)
    b11, b12, b21, b22 = split(Y)
    r = lambda P, Q: strassen(P, Q, counter, base)
    M1 = r(madd(a11, a22), madd(b11, b22))
    M2 = r(madd(a21, a22), b11)
    M3 = r(a11, madd(b12, b22, -1))
    M4 = r(a22, madd(b21, b11, -1))
    M5 = r(madd(a11, a12), b22)
    M6 = r(madd(a21, a11, -1), madd(b11, b12))
    M7 = r(madd(a12, a22, -1), madd(b21, b22))
    C11 = madd(madd(madd(M1, M4), M5, -1), M7)
    C12 = madd(M3, M5)
    C21 = madd(M2, M4)
    C22 = madd(madd(madd(M1, M2, -1), M3), M6)
    return join(C11, C12, C21, C22)


def block_check(rng, d=3):
    """2 x 2 block matrices with d x d integer blocks: Strassen's seven block products against the block product."""
    blk = lambda: [[rng.randint(-9, 9) for _ in range(d)] for _ in range(d)]
    A = [blk() for _ in range(4)]
    B = [blk() for _ in range(4)]
    big_A = [A[0][i] + A[1][i] for i in range(d)] + [A[2][i] + A[3][i] for i in range(d)]
    big_B = [B[0][i] + B[1][i] for i in range(d)] + [B[2][i] + B[3][i] for i in range(d)]
    cnt = [0]
    S = strassen(big_A, big_B, cnt, base=d)       # one Strassen level, blocks multiplied naively
    return S == mmul_naive(big_A, big_B), cnt[0]


def main():
    ok = True
    print("15a row [S1969]: omega <= log_2 7 ~ 2.8074 (Strassen's 2 x 2 algorithm with 7 multiplications, recursively)")
    good, out = identity_holds()
    print("  1. the four entries of AB as combinations of M1..M7 (coefficient tensors, factor order kept): %s" % good)
    for k in ("c11", "c12", "c21", "c22"):
        print("     %s = %s" % (k, " + ".join("%s*%s" % ab for ab in sorted(out[k]))))
    ok &= good
    rng = random.Random(20261112)
    bres = [block_check(rng) for _ in range(20)]
    bok = all(b for b, _ in bres)
    print("  2. 20 random 2 x 2 block matrices with 3 x 3 integer blocks (a non-commutative ring): equal to the block "
          "product in all 20: %s (each uses 7 block products = %d scalar multiplications)" % (bok, bres[0][1]))
    ok &= bok
    for k in range(1, 7):
        n = 2 ** k
        X = [[rng.randint(-5, 5) for _ in range(n)] for _ in range(n)]
        Y = [[rng.randint(-5, 5) for _ in range(n)] for _ in range(n)]
        cnt = [0]
        S = strassen(X, Y, cnt)
        eq = S == mmul_naive(X, Y)
        print("  3. k = %d (n = %d): recursive Strassen equals the ordinary product: %s; scalar multiplications %d = 7^%d: %s"
              % (k, n, eq, cnt[0], k, cnt[0] == 7 ** k))
        ok &= eq and cnt[0] == 7 ** k
    # log_2 7 enclosed: 2^(28073/10000) < 7 < 2^(28074/10000)  <=>  2^28073 < 7^10000 < 2^28074
    lo_ok = 2 ** 28073 < 7 ** 10000
    hi_ok = 7 ** 10000 < 2 ** 28074
    print("  4. 2^2.8073 < 7 < 2^2.8074 (exactly: 2^28073 < 7^10000 < 2^28074): %s, %s; so log_2 7 = 2.8073549..., "
          "printed 2.8074 (next to the exact log_2 7)" % (lo_ok, hi_ok))
    ok &= lo_ok and hi_ok
    # 30 digits of log_2 7: ln 7 and ln 2 correctly rounded at 60 digits, widened by 1e-50, then divided outward
    getcontext().prec = 60
    eps = Decimal(10) ** -50
    l7, l2 = Decimal(7).ln(), Decimal(2).ln()
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        lo = (l7 - eps) / (l2 + eps)
        c.rounding = ROUND_CEILING
        hi = (l7 + eps) / (l2 - eps)
    print("     log_2 7 in [%s, %s]" % (str(lo)[:32], str(hi)[:32]))
    ok &= Decimal("2.8073") < lo and hi < Decimal("2.8074")
    print("RESULT [S1969]: Strassen's seven products give a bilinear algorithm for 2 x 2 matrices over any ring; "
          "recursively 7^k multiplications for 2^k x 2^k, so omega <= log_2 7 = 2.80735...: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        for name, kw in (("M7 = (a12 + a22)(b21 + b22) (one sign changed)", dict(m7_sign=-1)),
                         ("C11 without M7", dict(variant="drop")),
                         ("C12 = M3 - M5", dict(variant="c12minus"))):
            g, _ = identity_holds(**kw)
            print("  %s: identity %s" % (name, "holds (UNEXPECTED)" if g else "fails, as expected"))
            c_ok &= not g
        claim = 7 ** 10000 < 2 ** 28073
        print("  a claimed omega <= 2.8073 from 7 products (needs 7 <= 2^2.8073): %s" %
              ("holds (UNEXPECTED)" if claim else "fails, as expected"))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
