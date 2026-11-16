#!/usr/bin/env python3
"""71a, lower bound [Hod2017] beta(1/2) > 6.4547837, program B (written separately from program A): the Fibonacci
closed forms of the paper. With b_m the Fibonacci numbers (b_m = (Phi^-m - (-Phi)^m)/sqrt 5, Phi = (sqrt 5 - 1)/2),
the paper gives q_m = (b_(m-1) z + b_m)/(b_m z + b_(m+1)) and pi_m = 1/(z b_m + b_(m+1)); at z = 1/2 these are
q_m = F_(m+2)/F_(m+3) and pi_m = 2/F_(m+3), and the two extra terms pi_(-2) = 1/(1 - z) = 2 = 2/F_1 and
pi_(-1) = 1/z = 2 = 2/F_2 fit the same pattern, so (n = k + 3)
    beta(1/2) = 4 (1/2)(1/2) sum_(n >= 1) 2 h(F_(n+1)/F_(n+2)) / F_n.
Checked: the closed forms against the recursion q_(m+1) = 1/(1 + q_m) for m <= 60 (exact). The sum in 70-digit
decimals (ln from the decimal module, correctly rounded, results moved outward by 10^-60), 240 terms, the tail bounded
by sum_(n > 240) 2/F_n <= 2/F_241 * 1/(1 - 1/phi) (h <= 1, F_n/F_(n+1) <= 2/3 for n >= 2). Red controls (--control):
the claim 6.4547838; the off-by-one F_(n+2)/F_(n+3) in place of F_(n+1)/F_(n+2) (an easy slip; the sum drops to about
6.38); the factor 4 z (1 - z) dropped (taken as 1/2). Standard library only.
"""
import sys
from decimal import Decimal as D, getcontext
from fractions import Fraction as Q

getcontext().prec = 70
PAD = D(10) ** -60
CONTROL = sys.argv[sys.argv.index("--control") + 1] if "--control" in sys.argv else None


def main():
    Fb = [0, 1]
    for _ in range(300):
        Fb.append(Fb[-1] + Fb[-2])
    # closed forms against the recursion at z = 1/2
    z = Q(1, 2); q = z; ok_cf = True; prod = Q(1)
    for m in range(1, 61):
        q = 1 / (1 + q)
        prod *= q
        ok_cf &= q == Q(Fb[m + 2], Fb[m + 3]) and prod == Q(2, Fb[m + 3])
        ok_cf &= q == (Fb[m - 1] * z + Fb[m]) / (Fb[m] * z + Fb[m + 1]) and prod == 1 / (z * Fb[m] + Fb[m + 1])
    print("closed forms q_m = F_(m+2)/F_(m+3), pi_m = 2/F_(m+3) agree with the recursion for m <= 60: %s" % ok_cf)
    ln2 = D(2).ln()
    def h(p):
        return -(p * p.ln() + (1 - p) * (1 - p).ln()) / ln2
    N = 240; s = D(0)
    sh = 1 if CONTROL == "shift" else 0
    for n in range(1, N + 1):
        p = D(Fb[n + 1 + sh]) / D(Fb[n + 2 + sh])
        s += 2 * h(p) / Fb[n]
    phi = (1 + D(5).sqrt()) / 2
    tail = D(2) / Fb[N + 1] / (1 - 1 / phi)
    c = D(1) / 2 if CONTROL == "factor" else D(1)          # 4 z (1 - z) = 1 at z = 1/2
    lo, hi = c * s - PAD, c * (s + tail) + PAD
    print("beta(1/2) in [%s, %s] (tail <= %.1e)" % (format(lo, ".30f"), format(hi, ".30f"), float(tail)))
    claim = D("6.4547838") if CONTROL == "claim" else D("6.4547837")
    ok = ok_cf and lo > claim
    print("beta > %s (the row): %s (margin %.3e)" % (claim, lo > claim, float(lo - claim)))
    print("key %s" % format(lo, ".40f")[:27])
    print("VERDICT " + (("REFUSED" if not ok else "NOT REFUSED") if CONTROL else ("PASS" if ok else "FAILED")))


if __name__ == "__main__":
    main()
