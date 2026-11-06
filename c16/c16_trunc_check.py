"""Independent check of the C16 lower bound 2.1562247706 (16a page, maximusrufus/c16-bgw-certificates, PR #326's row),
from the trial function the page states: u = sum_{0 < |k| <= R} e^{ik.x} / (|k|^2 (1 + mu |k|^2)), R = 4000, mu = 0.1221.
The certificate's own checker was neither read nor run; its R and mu are the only inputs, and its printed sums and bound
are compared at the end.

  python3 -I c16_trunc_check.py [CERT] [--controls]

As in c16_check.py (the definition on the 16a page): all c_k > 0, so |u|_inf = u(0) = S0 = sum c_k, and with
S1 = sum |k|^2 c_k^2, S2 = sum |k|^4 c_k^2, delta = S2/S1, L >= S0^2/(pi S1) - ln(delta) - ln(1 + ln(delta)).
With mu = p/q, c_k = q/(n d) for n = |k|^2 and d = q + p n, so the terms are q/(n d), q^2/(n d^2) and q^2/d^2. The three
sums run over the disc 0 < x^2 + y^2 <= R^2 by its 8-fold symmetry (0 <= y <= x), each term floored and ceiled in
integers at the scale 2^256, so each sum is enclosed exactly; the bound is then evaluated as in c16_check.py (decimal at
50 digits, the enclosure widened by 1e-30). Standard library only."""
import json
import os
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from c16_check import PI, SLACK  # noqa: E402  (pi by Machin's formula; the same slack as our other C16 bounds)

getcontext().prec = 50
K = 1 << 256


def sums(R, p, q):
    """[(lo, hi)] of K * S0, K * S1, K * S2 over 0 < x^2 + y^2 <= R^2, every term rounded down and up."""
    T0, T1 = K * q, K * q * q
    lo = [0, 0, 0]
    hi = [0, 0, 0]
    R2 = R * R
    for y in range(0, R + 1):
        for x in range(max(y, 1), R + 1):
            n = x * x + y * y
            if n > R2:
                break
            m = 4 if (y == 0 or y == x) else 8
            d = q + p * n
            for i, (num, den) in enumerate(((T0, n * d), (T1, n * d * d), (T1, d * d))):
                f, r = divmod(num, den)
                lo[i] += m * f
                hi[i] += m * (f + (1 if r else 0))
    return list(zip(lo, hi))


def bound(enc):
    (a0, b0), (a1, b1), (a2, b2) = [(D(lo) / D(K), D(hi) / D(K)) for lo, hi in enc]
    dlo, dhi = a2 / b1, b2 / a1
    B_lo = a0 ** 2 / (PI * b1) - dhi.ln() - (1 + dhi.ln()).ln() - SLACK
    B_hi = b0 ** 2 / (PI * a1) - dlo.ln() - (1 + dlo.ln()).ln() + SLACK
    return B_lo, B_hi, ((a0, b0), (a1, b1), (a2, b2))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cert = json.load(open(args[0])) if args else None
    R, (p, q) = (cert["R"], cert["mu"]) if cert else (4000, (1221, 10000))
    t0 = time.time()
    enc = sums(R, p, q)
    B_lo, B_hi, S = bound(enc)
    print(f"trial function: 0 < |k| <= R = {R}, mu = {p}/{q}; {time.time() - t0:.0f} s")
    for name, (a, b) in zip(("S0 = sum c_k", "S1 = sum |k|^2 c_k^2", "S2 = sum |k|^4 c_k^2"), S):
        print(f"  {name} in [{a:.40f}, {b:.40f}]")
    print(f"  L >= {B_lo:.20f}   (enclosure [{B_lo:.20f}, {B_hi:.20f}])")
    ok = B_lo > D("2.1562247706")
    print(f"  check L >= 2.1562247706, the printed bound: {'PASS' if ok else 'FAIL'}")
    if cert:
        scale = 1 << cert["PREC"]
        theirs = [(int(cert[f"S{i}lo_num"]), int(cert[f"S{i}hi_num"])) for i in (1, 2, 3)]
        agree = all(Fraction(lo, K) <= Fraction(th, scale) and Fraction(tl, scale) <= Fraction(hi, K)
                    for (lo, hi), (tl, th) in zip(enc, theirs))
        Lc = Fraction(*cert["L_lower_bound_rational"])
        print(f"  the certificate's three sum enclosures (scale 2^{cert['PREC']}) overlap ours: {'yes' if agree else 'NO'}")
        print(f"  its stated lower bound {float(Lc):.15f} lies in our enclosure: "
              f"{'yes' if D(Lc.numerator) / D(Lc.denominator) <= B_hi and D(Lc.numerator) / D(Lc.denominator) >= B_lo - D('1e-12') else 'NO'}")
        ok &= agree
    if "--controls" in sys.argv:
        print("Red controls (the same code, other trial functions):")
        for what, (Rc, pc, qc) in (("R = 400 (truncated ten times shorter)", (400, p, q)),
                                   ("mu = 0.14 instead of 0.1221", (R, 14, 100))):
            b = bound(sums(Rc, pc, qc))[0]
            refused = b < D("2.1562247706")
            print(f"  {what}: L >= {b:.12f}, {'below the printed bound, as expected' if refused else 'ABOVE (control failed)'}")
            ok &= refused
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
