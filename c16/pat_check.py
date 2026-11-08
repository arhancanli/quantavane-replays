"""Independent check of the 16a row [Pat2026], C_16 >= 2.156224 (#226): the trial function the page states,
u = sum_{0 < |k|^2 <= 4000^2} c_k e^{ik.x} with c_k = 1/(|k|^2 (|k|^2 + mu)), mu = 41/5.

Since 1/(n (n + mu)) = (1/mu) / (n (1 + n/mu)), u is a multiple of the trial function of c16_trunc_check.py with its mu
replaced by 1/mu = 5/41; the bound L >= S0^2/(pi S1) - ln(delta) - ln(1 + ln(delta)), delta = S2/S1, does not change
under scaling of the coefficients. The three sums are enclosed exactly by the same code (every term floored and ceiled
in integers at the scale 2^256, over the disc by its 8-fold symmetry). Standard library only.

  python3 -I pat_check.py [--controls] [--cutoffs]
"""
import os
import sys
import time
from decimal import Decimal as D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from c16_trunc_check import sums, bound  # noqa: E402

CLAIM = D("2.156224")


def main():
    t0 = time.time()
    R, p, q = 4000, 5, 41
    enc = sums(R, p, q)
    B_lo, B_hi, S = bound(enc)
    print(f"16a row [Pat2026]: trial function 0 < |k| <= R = {R}, c_k = 1/(|k|^2 (|k|^2 + 41/5)), i.e. mu' = {p}/{q} "
          f"in c16_trunc_check's normalization; {time.time() - t0:.0f} s")
    for name, (a, b) in zip(("S0 = sum c_k", "S1 = sum |k|^2 c_k^2", "S2 = sum |k|^4 c_k^2"), S):
        print(f"  {name} (scaled by mu) in [{a:.30f}, {b:.30f}]")
    (a1, b1), (a2, b2) = S[1], S[2]
    dlo, dhi = a2 / b1, b2 / a1
    print(f"  delta = S2/S1 in [{dlo:.12f}, {dhi:.12f}] (the page: approximately 3.93181): "
          f"{D('3.931805') <= dlo and dhi < D('3.931815')}")
    print(f"  L >= {B_lo:.20f}   (enclosure [{B_lo:.20f}, {B_hi:.20f}])")
    ok = B_lo >= CLAIM and D("3.931805") <= dlo and dhi < D("3.931815") and str(B_lo).startswith("2.1562247")
    print(f"  check L >= 2.156224, the row's bound, and the page's 2.1562247...: {'PASS' if ok else 'FAIL'}")
    if "--controls" in sys.argv:
        print("Red controls (the same code, other trial functions):")
        for what, (Rc, pc, qc) in (("R = 400 (truncated ten times shorter)", (400, p, q)),
                                   ("mu' = 0.14 in place of 5/41", (R, 14, 100))):
            b = bound(sums(Rc, pc, qc))[0]
            refused = b < CLAIM
            print(f"  {what}: L >= {b:.12f}, {'below the row bound, as expected' if refused else 'ABOVE (control failed)'}")
            ok &= refused
    if "--cutoffs" in sys.argv:
        # the row's last sentence: "Increasing the cutoff approaches the numerical value 2.15627"
        print("Other cutoffs, the same mu = 41/5 (each value rigorous):")
        vals = {R: B_lo}
        for Rc in (2000, 8000):
            vals[Rc] = bound(sums(Rc, p, q))[0]
        for Rc in (2000, 4000, 8000):
            print(f"  R = {Rc}: L >= {vals[Rc]:.12f}")
        g1, g2 = vals[4000] - vals[2000], vals[8000] - vals[4000]
        print(f"  gaps {g1:.4e} and {g2:.4e}, ratio {g1 / g2:.3f} (the tails of S0 and S2 are of order R^-2)")
        print(f"  extrapolated limit (an estimate, not a bound): {vals[8000] + g2 / 3:.10f}, "
              f"far below 2.15627 and below the certified upper bound 2.1562256466 on the same page")
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
