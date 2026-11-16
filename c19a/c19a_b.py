#!/usr/bin/env python3
"""19a [E1956]: Esseen's lower bound c_E ~ 0.4097321837, program B (written separately from program A).

The same definition (Esseen's lattice limit for normalized two-point laws, divided by beta_3), maximised without
calculus by a branch-and-bound over p in (0, 1/2] (the law with 1 - p is the mirror image): on [0, 1/2] the ratio
times sqrt(2 pi) is h(p) = (2 - p) / (3 (2p^2 - 2p + 1)). Each box gets the mean-value bound
h(p) <= h(m) + max|h'| |p - m|, with h' = (2p^2 - 8p + 3) / (3 (2p^2 - 2p + 1)^2) enclosed by monotonicity; boxes
whose bound is below the best value found are dropped. Exact Fractions throughout; pi by the
Bailey-Borwein-Plouffe series with its tail bound, square roots by the decimal module (correctly rounded, widened by
one unit). Standard library only.
"""
import sys
from fractions import Fraction as F
from decimal import Decimal, getcontext, Context, ROUND_FLOOR, ROUND_CEILING

getcontext().prec = 70
DN = Context(prec=70, rounding=ROUND_FLOOR)
UPc = Context(prec=70, rounding=ROUND_CEILING)


def ratio(p, mode):
    """sqrt(2 pi) times (Esseen's limit / beta_3) for the two-point law with P(X=1) = p <= 1/2"""
    u = 1 - 2 * p
    lat = F(1, 2) if mode != "lattice" else F(1)
    num = u / 6 + lat
    den = p * p + (1 - p) * (1 - p) if mode != "beta" else F(1)
    return num / den


def deriv_bounds(a, b, mode):
    """enclosure of the derivative of ratio(., mode) on [a, b] (by sampling-free monotone pieces)"""
    if mode == "main":
        # h = (2 - p)/(3 D), D = 2p^2 - 2p + 1; h' = (2p^2 - 8p + 3)/(3 D^2); on [0, 1/2]: numerator and D decrease
        n_lo, n_hi = 2 * b * b - 8 * b + 3, 2 * a * a - 8 * a + 3
        d_lo, d_hi = 2 * b * b - 2 * b + 1, 2 * a * a - 2 * a + 1
        c = [n_lo / (3 * d_lo ** 2), n_lo / (3 * d_hi ** 2), n_hi / (3 * d_lo ** 2), n_hi / (3 * d_hi ** 2)]
        return min(c), max(c)
    if mode == "beta":
        return F(-1, 3), F(-1, 3)
    # lattice: h = (u/6 + 1)/D = (7 - 2p)/(6 D)... derivative ( -2 D - (7 - 2p)(4p - 2) )/(6 D^2) = (4p^2 - 28p + 12)/(6 D^2)
    n_lo, n_hi = 4 * b * b - 28 * b + 12, 4 * a * a - 28 * a + 12
    d_lo, d_hi = 2 * b * b - 2 * b + 1, 2 * a * a - 2 * a + 1
    c = [n_lo / (6 * d_lo ** 2), n_lo / (6 * d_hi ** 2), n_hi / (6 * d_lo ** 2), n_hi / (6 * d_hi ** 2)]
    return min(c), max(c)


def maximise(mode, tol):
    boxes = [(F(0), F(1, 2))]
    best = max(ratio(F(0), mode), ratio(F(1, 2), mode))
    best_p = None
    upper = None
    it = 0
    while True:
        it += 1
        nxt = []
        upper = best
        for a, b in boxes:
            m = (a + b) / 2
            hm = ratio(m, mode)
            if hm > best:
                best, best_p = hm, m
            dl, dh = deriv_bounds(a, b, mode)
            ub = hm + max(abs(dl), abs(dh)) * (b - a) / 2
            if ub > best:
                nxt.append((a, m)); nxt.append((m, b))
                upper = max(upper, ub)
        boxes = [(a, b) for a, b in nxt if True]
        if upper - best < tol or not boxes:
            return best, upper, best_p, it, len(boxes)


def pi_iv(n=60):
    s = F(0)
    for k in range(n):
        s += F(1, 16 ** k) * (F(4, 8 * k + 1) - F(2, 8 * k + 4) - F(1, 8 * k + 5) - F(1, 8 * k + 6))
    tail = F(4, 16 ** n) * F(16, 15)          # the terms are positive and below 4/(8k+1)/16^k
    return s, s + tail


def sqrt_iv(q_lo, q_hi):
    a = Decimal(q_lo.numerator) / Decimal(q_lo.denominator)
    b = Decimal(q_hi.numerator) / Decimal(q_hi.denominator)
    lo = DN.next_minus(DN.next_minus(a).sqrt())
    hi = UPc.next_plus(UPc.next_plus(b).sqrt())
    return F(lo), F(hi)


def main():
    args = sys.argv[1:]
    ctrl = args[args.index("--control") + 1] if "--control" in args else None
    mode = ctrl if ctrl in ("lattice", "beta") else "main"
    printed = "0.4097321838" if ctrl == "print" else "0.4097321837"
    print("program B (branch-and-bound over p with mean-value bounds, exact Fractions)%s" % ((" CONTROL " + ctrl) if ctrl else ""))
    # self-test: the derivative enclosure contains a finite difference on a few boxes
    for a, b in ((F(0), F(1, 4)), (F(1, 10), F(1, 5)), (F(3, 10), F(1, 2))):
        dl, dh = deriv_bounds(a, b, "main")
        fd = (ratio(b, "main") - ratio(a, "main")) / (b - a)
        assert dl <= fd <= dh, "derivative enclosure"
    print("self-test: derivative enclosures contain the difference quotients: ok")
    best, upper, p, it, nb = maximise(mode, F(1, 10 ** 30))
    print("max over p of sqrt(2 pi) * ratio in [%s, %s] (%d levels, %d boxes left), at p ~ %s" %
          (dec(best, 32), dec(upper, 32), it, nb, dec(p, 20) if p is not None else "an end point"))
    p0, p1 = pi_iv()
    r_lo, r_hi = sqrt_iv(2 * p0, 2 * p1)
    c_lo, c_hi = best / r_hi, upper / r_lo
    print("c_E in [%s, %s]" % (dec(c_lo, 30), dec(c_hi, 30)))
    # the closed form (sqrt 10 + 3)/6 against the maximum
    s10_lo, s10_hi = sqrt_iv(F(10), F(10))
    cf_lo, cf_hi = (s10_lo + 3) / 6, (s10_hi + 3) / 6
    inside = best - F(1, 10 ** 30) <= cf_lo and cf_hi <= upper + F(1, 10 ** 30)
    print("closed form (sqrt 10 + 3)/6 = %s... %s the maximum's enclosure" % (dec(cf_lo, 30), "inside" if inside else "OUTSIDE"))
    d = F(printed); ulp = F(1, 10 ** 10)
    rnd = abs(c_lo - d) <= ulp / 2 and abs(c_hi - d) <= ulp / 2
    tr = d <= c_lo and c_hi < d + ulp
    print("row 19a [E1956]: c_E = %s...; page %s: %s" % (dec(c_lo, 30), printed,
          ", ".join([x for x, y in (("truncation", tr), ("correct rounding", rnd)) if y]) or "NONE"))
    good = rnd and (inside or ctrl)
    print("VERDICT " + (("REFUSED" if not good else "NOT REFUSED") if ctrl else ("PASS" if good else "FAILED")))


def dec(q, k):
    n = q.numerator * 10 ** k // q.denominator
    s = str(n).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


if __name__ == "__main__":
    main()
