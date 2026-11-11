#!/usr/bin/env python3
"""Exact replay of the 26b lower-bound row [DMPSS2014] (teorth/optimizationproblems, main 87b4ebb): B_{R,m} >=
2^(1 - 1/m) for every m >= 2, so C_26b >= 2, from the m-linear forms of Diniz, Munoz-Fernandez, Pellegrino and
Seoane-Sepulveda (arXiv 1111.3253v2, section 4): T_2(x, y) = x^1 y^1 + x^1 y^2 + x^2 y^1 - x^2 y^2 and
    T_m(x_1, ..., x_m) = (x_m^1 + x_m^2) T_{m-1}(x_1, ..., x_{m-1})
                         + (x_m^1 - x_m^2) T_{m-1}(B^{2^{m-2}} x_1, B^{2^{m-2}} x_2, B^{2^{m-3}} x_3, ..., B^2 x_{m-1}),
B the backward shift ((B^s x)^n = x^{n+s}).

For m = 2, ..., MMAX the coefficients T_m(e_{i_1}, ..., e_{i_m}) are built exactly by this recursion: every nonzero one
must be +-1, and there must be 4^(m-1) of them.  The norm ||T_m|| = sup over the cube (a multi-affine function, so a
maximum over sign vectors) is computed exhaustively for m <= MNORM (the last argument by sum_j |T(x_1, ..., e_j)|):
it must be 2^(m-1).  For larger m we check T_m(1, ..., 1) = 2^(m-1) (so ||T_m|| >= 2^(m-1)) and read the paper's
four-line induction for ||T_m|| <= 2^(m-1).  The ratio (sum |T(e)|^(2m/(m+1)))^((m+1)/(2m)) / ||T_m|| is then
(4^(m-1))^((m+1)/(2m)) / 2^(m-1) = 2^((m-1)/m) exactly (the exponent as a rational).  Integers only; no download.
usage: python3 -I dmpss_check.py [--controls]
"""
import itertools
import sys
from fractions import Fraction

CONTROLS = "--controls" in sys.argv
MMAX = 7
MNORM = 4


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def build(mmax, flip_sign=False):
    """{m: {index tuple (1-based): coefficient}}; flip_sign: the control, + instead of - in the second term."""
    T = {2: {(1, 1): 1, (1, 2): 1, (2, 1): 1, (2, 2): -1}}
    for m in range(3, mmax + 1):
        shifts = [2 ** (m - 2)] + [2 ** (m - k) for k in range(2, m)]      # for x_1, x_2, ..., x_{m-1}
        out = {}
        for key, c in T[m - 1].items():
            shifted = tuple(i + s for i, s in zip(key, shifts))
            for k2, v in ((key + (1,), c), (key + (2,), c), (shifted + (1,), c), (shifted + (2,), c if flip_sign else -c)):
                out[k2] = out.get(k2, 0) + v
        T[m] = {k: v for k, v in out.items() if v}
    return T


def norm(Tm, m):
    """max over sign vectors of |T_m(x_1, ..., x_m)|, the last argument by sum_j |T(x_1, ..., x_{m-1}, e_j)|."""
    dims = [max(k[a] for k in Tm) for a in range(m)]
    terms = list(Tm.items())
    best = 0
    for xs in itertools.product(*[itertools.product((1, -1), repeat=dims[a]) for a in range(m - 1)]):
        acc = {}
        for key, c in terms:
            v = c
            for a in range(m - 1):
                v *= xs[a][key[a] - 1]
            acc[key[-1]] = acc.get(key[-1], 0) + v
        best = max(best, sum(abs(v) for v in acc.values()))
    return best


def main():
    ok = True
    T = build(MMAX)
    print("m  coordinates per argument      nonzero coefficients  all +-1  ||T_m||                ratio = 2^((m-1)/m)")
    for m in range(2, MMAX + 1):
        Tm = T[m]
        dims = [max(k[a] for k in Tm) for a in range(m)]
        unit = all(abs(v) == 1 for v in Tm.values())
        at_ones = sum(Tm.values())
        if m <= MNORM:
            nm = norm(Tm, m)
            ntext = f"{nm} (exhaustive)"
            good_n = nm == 2 ** (m - 1)
        else:
            ntext = f">= {at_ones} at (1,..,1)"
            good_n = at_ones == 2 ** (m - 1)
        p = Fraction(2 * m, m + 1)
        expo = Fraction(2 * (m - 1)) / p - (m - 1)      # log2 of (4^(m-1))^(1/p) / 2^(m-1)
        good = unit and len(Tm) == 4 ** (m - 1) and good_n and expo == Fraction(m - 1, m) and max(dims) == 2 ** (m - 1)
        ok &= good
        print(f"{m}  {str(dims):<28} {len(Tm):<21} {str(unit):<8} {ntext:<22} 2^({expo}) = {2 ** float(expo):.12f}"
              f"{'' if good else '   FAILS'}")
    print(f"RESULT [DMPSS2014]: B_R,m >= 2^(1-1/m) for m = 2..{MMAX} ({'exhaustive norms for m <= ' + str(MNORM)}; "
          f"above, the norm's upper bound is the paper's induction), increasing to 2, so C_26b >= 2: "
          f"{'HOLDS' if ok else 'FAILS'}")
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        Tb = build(3, flip_sign=True)[3]
        nb = norm(Tb, 3)
        print(f"  T_3 with + for - in the second term: {len(Tb)} nonzero coefficients, norm {nb}, ratio "
              f"{(sum(abs(v) ** 1.5 for v in Tb.values())) ** (2 / 3) / nb:.6f} (not 2^(2/3) = {2 ** (2 / 3):.6f})")
        res.append(("T_3 with + for -", nb != 4))
        T2 = dict(T[2])
        T2[(2, 2)] = 1                                   # x^1 y^1 + x^1 y^2 + x^2 y^1 + x^2 y^2 = (x^1 + x^2)(y^1 + y^2)
        res.append(("T_2 with +x^2y^2", norm(T2, 2) != 2))
        print(f"  T_2 with + x^2 y^2: norm {norm(T2, 2)}, ratio {(4 ** 0.75) / norm(T2, 2):.6f} (< sqrt 2)")
        res.append(("a claimed exponent m/(m+1) at m = 3", Fraction(4) / Fraction(3, 2) - 2 != Fraction(3, 4)))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
