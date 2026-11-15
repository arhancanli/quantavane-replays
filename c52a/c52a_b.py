#!/usr/bin/env python3
"""52a, program B (written separately from program A): [FP1983] 5.191 from the page's model.

Route B: for each assignment a of n variables, count the clauses (three distinct variables, any signs) that a
satisfies; with m clauses drawn independently and uniformly, E[#satisfying assignments] = sum_a (s_a / N)^m,
N = 8 C(n,3). The program checks s_a = 7 C(n,3) for every a and n = 3..9, so the sum is 2^n (7/8)^m for every m.
Then r1 = log 2 / log(8/7) in decimal intervals (40 digits; ln correctly rounded in the decimal module, widened by
one unit in the last place) and r1 <= claim. Standard library only.
Usage: c52a_b.py [--claim R] [--variant two-literal]
"""
import sys
from decimal import Decimal as D, localcontext, Context, ROUND_HALF_EVEN, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as Q
from math import comb

def satisfied_count(n, a, variant):
    """number of clauses satisfied by the assignment a (bit i of a = value of variable i)"""
    cnt = 0
    if variant == "two-literal":              # control: clauses on two distinct variables
        for i in range(n):
            for j in range(i + 1, n):
                for sg in range(4):
                    cnt += ((a >> i) & 1) == (sg & 1) or ((a >> j) & 1) == ((sg >> 1) & 1)
        return cnt, 4 * comb(n, 2)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for sg in range(8):
                    lits = ((i, sg & 1), (j, (sg >> 1) & 1), (k, (sg >> 2) & 1))
                    cnt += any(((a >> v) & 1) == s for v, s in lits)
    return cnt, 8 * comb(n, 3)

def main():
    a = sys.argv[1:]
    opt = dict(zip(a[0::2], a[1::2]))
    claim = opt.get("--claim", "5.191")
    variant = opt.get("--variant")
    if variant:
        print("CONTROL: clause model", variant)
    print("52a program B: [FP1983] r_3 <= %s" % claim)
    ok = True
    for n in range(3, 10):
        vals = set()
        for asg in range(2 ** n):
            s, N = satisfied_count(n, asg, variant)
            vals.add(Q(s, N))
        good = vals == {Q(7, 8)}
        ok &= good
        print("n = %d: every assignment satisfies a fraction %s of the clauses: %s"
              % (n, ", ".join(str(v) for v in sorted(vals)), good))
        if not good:
            break
    if not ok:
        print("VERDICT REFUSED: the probability that a random clause is satisfied is not 7/8"); return
    with localcontext(Context(prec=40, rounding=ROUND_HALF_EVEN)) as c:
        l2 = D(2).ln(c); l87 = (D(8) / D(7)).ln(c)
        l2 = (l2.next_minus(c), l2.next_plus(c)); l87 = (l87.next_minus(c), l87.next_plus(c))
    with localcontext(Context(prec=40, rounding=ROUND_FLOOR)):
        r_lo = l2[0] / l87[1]
    with localcontext(Context(prec=40, rounding=ROUND_CEILING)):
        r_hi = l2[1] / l87[0]
    print("r1 = log 2 / log(8/7) in [%s, %s]" % (r_lo, r_hi))
    ok = r_hi < D(claim)
    print("row [FP1983] r_3 <= r1 <= %s: %s" % (claim, ok))
    print("VERDICT PASS" if ok else "VERDICT REFUSED: r1 is not below %s" % claim)

if __name__ == "__main__":
    main()
