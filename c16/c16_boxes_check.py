"""The 16a upper bound C16 <= 2.1562256466* (maximusrufus, added in PR #246), from its certificate's delta-boxes.

Our own code (Quantavane replays): the Cauchy-Schwarz bound and the lattice sums of c16_check.py, which shares no code
with the certificate's checkers (verify_c16_upper.py and certify_c16_upper.py were neither read nor run). For every
mu > 0 and delta >= 1, the bound U(delta) = G(mu)(1 + mu delta)/pi - ln delta - ln(1 + ln delta) is convex in delta, so
on a box [d1, d2] it is at most its larger endpoint value (evaluated at d1 rounded down and d2 rounded up, a larger box).

The certificate's boxes [d1, d2, mu, bound] (d1, d2 exact rationals, mu a decimal string, read exactly) must tile
[1, D0] exactly.  Every box below 2^18 is checked with its own mu; for delta >= 2^18 (the boxes above, the certificate's
own tail beyond D0, and everything else) our own lattice-count tail bound of c16_check.py is used instead, so the claim
is decided without the certificate's tail argument.

usage: python3 -I c16_boxes_check.py upper_certificate.json [--controls]
"""
import json
import os
import sys
import time
from decimal import Decimal as D, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)            # for c16_check.py next to this file (python3 -I leaves the script's folder out)
from c16_check import G_bracket, U_hi, n_rows, tail_bound  # noqa: E402

TWO18 = Fraction(2 ** 18)


def dec(fr, rounding):
    with localcontext() as ctx:
        ctx.prec = 50
        ctx.rounding = rounding
        return D(fr.numerator) / D(fr.denominator)


def read(cert):
    boxes = []
    for d1, d2, mu, _bound in cert['boxes']:
        boxes.append((Fraction(d1), Fraction(d2), D(mu)))
    return boxes, Fraction(cert['D0'])


def check(boxes, D0, T, max_rows=1 << 17, width=D('1e-19')):
    s = sorted(boxes)
    if s[0][0] != 1 or s[-1][1] != D0:
        return False, 'boxes do not run from 1 to D0', None
    for (a1, a2, _), (b1, b2, _) in zip(s, s[1:]):
        if a2 != b1:
            return False, 'gap or overlap at %s' % float(a2), None
    worst, used = None, 0
    for i, (d1, d2, mu) in enumerate(s):
        if d1 >= TWO18:
            continue
        if not (d1 < d2 and mu > 0):
            return False, 'bad box %d' % i, None
        right = min(d2, TWO18)
        lo_d, hi_d = dec(d1, ROUND_FLOOR), dec(right, ROUND_CEILING)
        n = n_rows(1 / mu, 2048)
        while True:
            lo, hi = G_bracket(mu, n)
            u = max(U_hi(hi, mu, lo_d), U_hi(hi, mu, hi_d))
            # boxes clearly below T stop early; the others are refined until G is enclosed to `width`, so the
            # largest endpoint value is reproduced to about that, not only shown to be below T
            if u <= T - D('1e-6') or hi - lo < width or n >= max_rows:
                break
            n *= 2
        used += 1
        if worst is None or u > worst[0]:
            worst = (u, i, float(d1), float(right), mu, n)
        if u > T:
            return False, 'box %d [%s, %s] mu %s: bound %s > %s' % (i, float(d1), float(right), mu, u, T), worst
    tb, A = tail_bound()
    if tb > T:
        return False, 'tail bound %s > T' % tb, worst
    return True, '%d boxes that start below 2^18 checked with their own mu (up to 2^18); delta >= 2^18 by our tail bound, U <= %.12f' % (used, tb), worst


def main():
    t0 = time.time()
    cert = json.load(open(sys.argv[1]))
    boxes, D0 = read(cert)
    T = D('2.1562256466')
    print('certificate: %d boxes on [1, %s]; the stated U = %s, rational %s/%s'
          % (len(boxes), cert['D0'], cert['U'], cert['U_rational'][0][:12] + '...', cert['U_rational'][1][:12] + '...'))
    ok, msg, worst = check(boxes, D0, T)
    u, i, d1, d2, mu, n = worst
    if not ok:
        print('  ' + msg)
        print('RESULT 16a upper (maximusrufus): C16 <= 2.1562256466 FAILED')
        sys.exit(1)
    print('  ' + msg)
    print('  largest endpoint bound %.20f, box [%.9f, %.9f] (mu = %s, %d rows); margin below 2.1562256466: %.3e'
          % (u, d1, d2, mu, n, T - u))
    print('RESULT 16a upper (maximusrufus): C16 <= 2.1562256466 %s' % ('CERTIFIED' if ok else 'FAILED'))
    if '--controls' in sys.argv:
        r1 = check(boxes, D0, u - D('1e-15'))[0]
        print('control: the bound 1e-15 below our largest endpoint value refused:', not r1)
        s = sorted(boxes)
        bad = [(a, b, (m * D('1.01')) if j == i else m) for j, (a, b, m) in enumerate(s)]
        r2 = check(bad, D0, T)[0]
        print('control: mu of box %d raised by 1%% refused:' % i, not r2)
        r3 = check(s[:i] + s[i + 1:], D0, T)[0]
        print('control: box %d removed (a gap) refused:' % i, not r3)
        ok = ok and not r1 and not r2 and not r3
    print('run time %.0f s' % (time.time() - t0))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
