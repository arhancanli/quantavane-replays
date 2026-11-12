#!/usr/bin/env python3
"""Replay of two more 3a rows (teorth/optimizationproblems, main 87b4ebb), one on each side:

  [GHR2007] lower, 1.1078...: the lemma the page states, C_3a >= 1 + log(|U-U|/|U+U|)/log(2 max U + 1) for a finite
            set U of non-negative integers containing 0 with |U-U| < 2 max U + 1, applied to U = {0,1,3,6,13,17,21}
            as printed in the row (|U+U| = 26, |U-U| = 39, q = 43). Counted exactly by two routes (Python sets; sorted
            lists of all pair values), the value enclosed in 60-digit decimal interval arithmetic. The lemma is read.
  [K2026c]  upper, (13524e - 7359)/(9451e - 3286) = 1.312373302115...: the page states the two-set exponent
            lambda = (9451e - 3286)/(5378e + 787) = 1.454277448906... and the transfer C_3a <= 2 - 1/lambda. We check
            the transfer identity 2 - 1/lambda = (13524e - 7359)/(9451e - 3286) exactly (as linear forms in e), both
            decimals in interval arithmetic (e = exp(1) at 60 digits, widened by 1e-50), that lambda < 3/2 and the bound
            < 4/3 ([GHR2007]'s upper bound), and the README cell: 1.31237330 is below the value, so as an upper bound
            it is rounded the wrong way; PR #364's 1.31237331 is above it. The theorem |X - Y| <= |X + Y|^lambda and
            the multiscale transfer are read, not checked (the row records a Lean formalization).
Red controls (--controls): U with 21 replaced by 22 gives a different count and a value below 1.1078; a claim 1.1079
is refused; the identity with 13523 for 13524 is refused; 1.31237330 as an upper bound is refused.
usage: python3 -I ghr_k2026c_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING

getcontext().prec = 60
EPS = Decimal(10) ** -50
CONTROLS = "--controls" in sys.argv


def counts(U):
    U = sorted(set(U))
    s = {x + y for x in U for y in U}
    d = {x - y for x in U for y in U}
    pairs_s = sorted(x + y for i, x in enumerate(U) for y in U[i:])
    pairs_d = sorted(x - y for x in U for y in U)
    nd = lambda xs: sum(1 for i, x in enumerate(xs) if i == 0 or x != xs[i - 1])
    assert len(s) == nd(pairs_s) and len(d) == nd(pairs_d), "the two counting routes disagree"
    return len(s), len(d), 2 * max(U) + 1


def lemma(U):
    s, d, q = counts(U)
    ok = 0 in U and min(U) == 0 and d < q
    num = (Decimal(d) / Decimal(s)).ln()
    den = Decimal(q).ln()
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        lo = 1 + (num - EPS) / (den + EPS)
        c.rounding = ROUND_CEILING
        hi = 1 + (num + EPS) / (den - EPS)
    return s, d, q, ok, lo, hi


def lin_e(a, b):
    """a e + b as a pair of integers."""
    return (a, b)


def ratio_interval(p, q, e_lo, e_hi):
    """(p0 e + p1)/(q0 e + q1) enclosed for e in [e_lo, e_hi] (both linear forms positive there; monotone in e)."""
    vals = []
    for e in (e_lo, e_hi):
        vals.append((Decimal(p[0]) * e + p[1]) / (Decimal(q[0]) * e + q[1]))
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        lo = min(vals) - EPS
        c.rounding = ROUND_CEILING
        hi = max(vals) + EPS
    return lo, hi


def k2026c(c13524=13524):
    lam_num, lam_den = lin_e(9451, -3286), lin_e(5378, 787)
    # 2 - 1/lambda = 2 - lam_den/lam_num = (2 lam_num - lam_den)/lam_num
    lhs = (2 * lam_num[0] - lam_den[0], 2 * lam_num[1] - lam_den[1])
    row = lin_e(c13524, -7359)
    ident = lhs == row
    e = Decimal(1).exp()
    e_lo, e_hi = e - EPS, e + EPS
    lam = ratio_interval(lam_num, lam_den, e_lo, e_hi)
    val = ratio_interval(row, lam_num, e_lo, e_hi)
    return ident, lhs, lam, val


def main():
    ok = True
    U = [0, 1, 3, 6, 13, 17, 21]
    s, d, q, hyp, lo, hi = lemma(U)
    print("3a row [GHR2007] lower: U = %s" % U)
    print("  |U+U| = %d, |U-U| = %d, q = 2 max U + 1 = %d (the row: 26, 39, 43): %s; 0 in U and |U-U| < q: %s"
          % (s, d, q, (s, d, q) == (26, 39, 43), hyp))
    print("  1 + log(39/26)/log(43) in [%s, %s]" % (str(lo)[:24], str(hi)[:24]))
    ok1 = (s, d, q) == (26, 39, 43) and hyp and Decimal("1.1078") <= lo and hi < Decimal("1.1079")
    print("RESULT [GHR2007] 1.1078...: C_3a >= 1.10780... (the printed 1.1078... is its truncation): %s"
          % ("HOLDS" if ok1 else "FAILS"))
    ok &= ok1
    ident, lhs, lam, val = k2026c()
    print("3a row [K2026c] upper: lambda = (9451e - 3286)/(5378e + 787), bound 2 - 1/lambda")
    print("  2 - 1/lambda = (%de - %d)/(9451e - 3286) and the row's (13524e - 7359)/(9451e - 3286): identical: %s"
          % (lhs[0], -lhs[1], ident))
    print("  lambda in [%s, %s] (the page: 1.454277448906...); lambda < 3/2: %s"
          % (str(lam[0])[:26], str(lam[1])[:26], lam[1] < Decimal("1.5")))
    print("  bound in [%s, %s] (the row: 1.312373302115...); below 4/3: %s"
          % (str(val[0])[:26], str(val[1])[:26], val[1] < Decimal(4) / 3))
    digits_ok = str(val[0])[:14] == "1.312373302115" and str(val[1])[:14] == "1.312373302115" and \
        str(lam[0])[:14] == "1.454277448906" and str(lam[1])[:14] == "1.454277448906"
    readme_old = Decimal("1.31237330") >= val[1]
    readme_new = Decimal("1.31237331") >= val[1]
    print("  printed digits 1.312373302115 and 1.454277448906 are truncations: %s" % digits_ok)
    print("  README cell 1.31237330 >= the bound (valid as an upper bound): %s; PR #364's 1.31237331: %s"
          % (readme_old, readme_new))
    ok2 = ident and digits_ok and lam[1] < Decimal("1.5") and val[1] < Decimal(4) / 3 and (not readme_old) and readme_new
    print("RESULT [K2026c] 1.312373302115...: the transfer identity and the decimals hold; the README's 1.31237330 is "
          "below the bound by %.4e (rounded the wrong way), 1.31237331 is valid: %s"
          % (val[0] - Decimal("1.31237330"), "HOLDS" if ok2 else "FAILS"))
    ok &= ok2
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        s2, d2, q2, hyp2, lo2, hi2 = lemma([0, 1, 3, 6, 13, 17, 22])
        print("  U with 22 for 21: |U+U| = %d, |U-U| = %d, q = %d, value %s < 1.1078: %s"
              % (s2, d2, q2, str(hi2)[:10], hi2 < Decimal("1.1078")))
        c_ok &= hi2 < Decimal("1.1078")
        claim = Decimal("1.1079") <= lo
        print("  a claimed 1.1079 from U: %s" % ("holds (UNEXPECTED)" if claim else "refused, as expected"))
        c_ok &= not claim
        ident2, _, _, _ = k2026c(13523)
        print("  the transfer identity with 13523e - 7359: %s" % ("holds (UNEXPECTED)" if ident2 else "refused, as expected"))
        c_ok &= not ident2
        print("  1.31237330 as an upper bound: %s" % ("holds (UNEXPECTED)" if readme_old else "refused, as expected"))
        c_ok &= not readme_old
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
