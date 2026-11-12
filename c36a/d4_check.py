#!/usr/bin/env python3
"""Exact replay of the 36a lower-bound row [dLOV2014] (teorth/optimizationproblems, main 87b4ebb):
Delta_4 >= pi^2/16 ~ 0.616850, the D_4 lattice packing (center density 1/8, Table 1 of de Laat, Oliveira and
Vallentin, as the row quotes it; the lattice is the classical D_4).

What is checked, exactly:
  1. D_4 = {x in Z^4 : x_1 + x_2 + x_3 + x_4 even}; the basis (1,1,0,0), (1,-1,0,0), (0,1,-1,0), (0,0,1,-1) lies in
     D_4 and has determinant +-2 = [Z^4 : D_4] (index 2: D_4 is the kernel of the sum mod 2), so it is a basis and the
     covolume is 2 (Gram determinant 4);
  2. the minimal norm is 2: a nonzero integer vector of norm 1 is +-e_i, whose sum is odd; the 24 vectors of norm 2
     (+-e_i +- e_j) are listed by enumeration over [-2, 2]^4;
  3. balls of radius r = sqrt(2)/2 (half the minimal distance) centred at D_4 do not overlap, so the density is
     vol(B^4) r^4 / covolume = (pi^2/2)(1/4)/2 = pi^2/16, center density r^4/covolume = 1/8 = 0.12500;
  4. pi^2/16 enclosed through pi from Machin's formula with alternating-series bounds (exact rationals):
     0.61685027506808491367...; the printed 0.616850 stands next to the exact pi^2/16.
Also (not a replay of those rows; for PR #364): the three upper-bound rows convert center densities read from the
cited tables, Delta_4 <= (pi^2/2) delta_4 with delta_4 = 0.13127, 0.13126 [CE2003] and 0.130587 [dLOV2014]; the
products are enclosed and compared with the printed decimals and with the README cell 0.644421 (PR #364 proposes
0.644422).
Red controls (--controls): Z^4 (minimal norm 1, covolume 1) gives pi^2/32, below the row; a claim 0.6169 (above
pi^2/16) is refused; the vector (1,0,0,0) added to the basis set is refused as not in D_4; 0.644421 as an upper bound
for (pi^2/2) 0.130587 is refused.
usage: python3 -I d4_check.py [--controls]      (Python 3.9+, standard library only)
"""
import itertools
import sys
from fractions import Fraction as F

CONTROLS = "--controls" in sys.argv


def det(M):
    M = [[F(x) for x in row] for row in M]
    n, d = len(M), F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d


def arctan_bounds(x, terms):
    """lo <= arctan(x) <= hi for 0 < x < 1 (alternating series, exact rationals)."""
    s, lo, hi = F(0), None, None
    for k in range(terms):
        s += F((-1) ** k) * x ** (2 * k + 1) / (2 * k + 1)
        if k % 2 == 0:
            hi = s
        else:
            lo = s
    return lo, hi


def pi_bounds():
    l5, h5 = arctan_bounds(F(1, 5), 40)
    l239, h239 = arctan_bounds(F(1, 239), 14)
    return 16 * l5 - 4 * h239, 16 * h5 - 4 * l239


def dec(x, n=20):
    """floor of x to n decimals, as a string."""
    q = (x.numerator * 10 ** n) // x.denominator
    s = str(q).rjust(n + 1, "0")
    return s[:-n] + "." + s[-n:]


def in_lattice(v):
    return sum(v) % 2 == 0


def lattice_data(member, basis):
    pts = [v for v in itertools.product(range(-2, 3), repeat=4) if any(v) and member(v)]
    norms = sorted({sum(x * x for x in v) for v in pts})
    mn = norms[0]
    minimal = [v for v in pts if sum(x * x for x in v) == mn]
    return mn, len(minimal), abs(det(basis))


def main():
    ok = True
    print("36a row [dLOV2014] (lower): Delta_4 >= pi^2/16 ~ 0.616850 from the D_4 lattice")
    basis = [(1, 1, 0, 0), (1, -1, 0, 0), (0, 1, -1, 0), (0, 0, 1, -1)]
    in_ok = all(in_lattice(b) for b in basis)
    mn, nmin, cov = lattice_data(in_lattice, basis)
    print("  1. basis vectors in D_4: %s; |det| = %s = index [Z^4 : D_4] = 2: %s; Gram determinant %s"
          % (in_ok, cov, cov == 2, cov * cov))
    ok &= in_ok and cov == 2
    one_norm_odd = all(sum(v) % 2 == 1 for v in itertools.product(range(-1, 2), repeat=4) if sum(x * x for x in v) == 1)
    print("  2. every integer vector of norm 1 has odd sum (not in D_4): %s; least norm found %d, with %d minimal vectors "
          "(+-e_i +- e_j: 24): %s" % (one_norm_odd, mn, nmin, mn == 2 and nmin == 24))
    ok &= one_norm_odd and mn == 2 and nmin == 24
    r4 = F(mn, 4) ** 2                       # r^4 with r^2 = mn/4
    center = r4 / cov
    plo, phi = pi_bounds()
    dens_lo, dens_hi = plo * plo / 2 * r4 / cov, phi * phi / 2 * r4 / cov
    print("  3. r^2 = %s, r^4 = %s, center density r^4/covolume = %s (the table's 0.12500): %s"
          % (F(mn, 4), r4, center, center == F(1, 8)))
    ok &= center == F(1, 8)
    print("  4. pi in [%s, %s]; density (pi^2/2) r^4 / covolume = pi^2/16 in [%s, %s]"
          % (dec(plo, 22), dec(phi, 22), dec(dens_lo), dec(dens_hi)))
    ok &= dec(dens_lo, 6) == "0.616850" and dens_lo > F(616850, 10 ** 6)
    print("     printed 0.616850 is pi^2/16 truncated to 6 decimals: %s" % (dec(dens_lo, 6) == "0.616850"))
    print("RESULT [dLOV2014] lower: the D_4 packing has density exactly pi^2/16 = %s..., so Delta_4 >= pi^2/16: %s"
          % (dec(dens_lo, 15), "HOLDS" if ok else "FAILS"))
    # conversions of the upper-bound rows (PR #364)
    print("Upper-bound rows, conversion Delta_4 <= (pi^2/2) delta_4 only (the delta_4 values are read from the tables):")
    for key, delta, printed in (("[CE2003] Rogers", "0.13127", "0.647791"), ("[CE2003] new", "0.13126", "0.647742"),
                                ("[dLOV2014]", "0.130587", "0.644421")):
        d = F(delta)
        lo, hi = plo * plo / 2 * d, phi * phi / 2 * d
        tr = dec(lo, 6) == printed and dec(hi, 6) == printed
        print("  %-16s (pi^2/2) %s in [%s, %s]; printed %s... is its truncation: %s"
              % (key, delta, dec(lo, 15), dec(hi, 15), printed, tr))
        ok_conv = tr
        if key == "[dLOV2014]":
            readme_ok = F(644421, 10 ** 6) >= hi
            pr_ok = F(644422, 10 ** 6) >= hi
            print("  README cell 0.644421 >= the value (valid as an upper bound): %s; PR #364's 0.644422 >= the value: %s"
                  % (readme_ok, pr_ok))
            ok_conv &= (not readme_ok) and pr_ok
        ok &= ok_conv
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        mnz, nz, covz = lattice_data(lambda v: True, [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)])
        dz = phi * phi / 2 * F(mnz, 4) ** 2 / covz
        print("  Z^4: least norm %d, covolume %s, density <= %s < pi^2/16: %s" % (mnz, covz, dec(dz, 6), dz < dens_lo))
        c_ok &= dz < dens_lo
        claim = F(6169, 10000) <= dens_lo
        print("  a claimed Delta_4 >= 0.6169 from D_4: %s" % ("holds (UNEXPECTED)" if claim else "refused, as expected"))
        c_ok &= not claim
        bad = in_lattice((1, 0, 0, 0))
        print("  (1,0,0,0) in D_4: %s" % ("yes (UNEXPECTED)" if bad else "no, refused as expected"))
        c_ok &= not bad
        r = F(644421, 10 ** 6) >= phi * phi / 2 * F("0.130587")
        print("  0.644421 as an upper bound for (pi^2/2) 0.130587: %s" % ("holds (UNEXPECTED)" if r else "refused, as expected"))
        c_ok &= not r
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
