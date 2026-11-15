#!/usr/bin/env python3
"""32a, program B (written separately from program A): the [Li25] closed form ~ 0.76 (see program A).
Route B: decimal arithmetic at 60 digits with directed rounding; log by the decimal module's ln (correctly rounded,
widened by one unit in the last place), divided by ln 2 to get log2; terms k <= 300; tail for k > 300 bounded by
sum_{k>300} k^2 2^(-k-1) <= 301^2 2^(-301) (ratio of consecutive terms below 0.51). Standard library only.
Usage: c32a_b.py [--print D] [--variant no-half]
"""
import sys
from decimal import Decimal as D, localcontext, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

P = 60
def ln_bounds(x):
    with localcontext(Context(prec=P, rounding=ROUND_HALF_EVEN)) as c:
        v = D(x).ln(c)
        return v.next_minus(c), v.next_plus(c)

def main():
    arg = sys.argv[1:]
    opt = dict(zip(arg[0::2], arg[1::2]))
    printed = opt.get("--print", "0.76"); variant = opt.get("--variant")
    if variant:
        print("CONTROL: 2^(-k) in place of 2^(-k-1)")
    print("32a program B")
    l2 = ln_bounds(2)
    e = 0 if variant == "no-half" else 1
    lo = D(0); hi = D(0)
    for k in range(2, 301):
        lk = ln_bounds(k)
        with localcontext(Context(prec=P, rounding=ROUND_FLOOR)):
            lo += D(k) * lk[0] / l2[1] / D(2) ** (k + e)
        with localcontext(Context(prec=P, rounding=ROUND_CEILING)):
            hi += D(k) * lk[1] / l2[0] / D(2) ** (k + e)
    with localcontext(Context(prec=P, rounding=ROUND_CEILING)):
        hi += D(301) ** 2 / D(2) ** 301
    # log2 log2 e = -ln(ln 2)/ln 2; ln 2 in (0.69, 0.70): ln(ln 2) negative, decreasing to the left
    ll_lo = ln_bounds(l2[0])[0]; ll_hi = ln_bounds(l2[1])[1]
    with localcontext(Context(prec=P, rounding=ROUND_FLOOR)):
        a_lo = (-ll_hi) / l2[1]
    with localcontext(Context(prec=P, rounding=ROUND_CEILING)):
        a_hi = (-ll_lo) / l2[0]
    with localcontext(Context(prec=P, rounding=ROUND_FLOOR)):
        v_lo = lo - a_hi
    with localcontext(Context(prec=P, rounding=ROUND_CEILING)):
        v_hi = hi - a_lo
    print("value in [%s, %s]" % (str(v_lo)[:30], str(v_hi)[:30]))
    d = len(printed.split(".")[1]); p = D(printed); h = D(5) / D(10) ** (d + 1)
    ok = p - h <= v_lo and v_hi <= p + h
    print("row [Li25] ~ %s: %s" % (printed, ok))
    print("VERDICT PASS" if ok else "VERDICT REFUSED: the printed value does not match")

if __name__ == "__main__":
    main()
