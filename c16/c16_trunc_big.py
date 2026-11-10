"""The 16a row [Pat2026b] of PR #270, C16 >= 2.156225504* (vidigoat), from the trial function its row states: the same
truncated conditional extremal as the [Pat2026] row (c_k = 1/(|k|^2 (|k|^2 + 41/5)), which our replay on PR #226 checked
at the cutoff |k| <= 4000), with the cutoff 0 < |k| <= 40000 (5,026,547,528 lattice points). The certificate's C and
Python checkers were neither read nor run; only R and mu from the row are used.

  python3 -I c16_trunc_big.py BINARY [--controls]      (BINARY: c16_trunc_fast, built from c16_trunc_fast.c)

The sums are those of c16_trunc_check.py (same definitions; up to a constant factor, which the bound does not see,
c_k = q/(n (q + p n)) with p/q = 5/41), computed by c16_trunc_fast.c at the scale 2^100, every term rounded down and up,
so each sum is enclosed exactly; the bound S0^2/(pi S1) - ln(delta) - ln(1 + ln(delta)), delta = S2/S1, is evaluated as
in c16_check.py (decimal at 50 digits, pi by Machin's formula, the enclosure widened by 1e-30).
Self-test first: the C sums equal the Python function sums() of c16_trunc_check.py, integer for integer, at the scale
2^100 for R = 300 and R = 1200. Standard library only (and a C compiler)."""
import os
import subprocess
import sys
import time
from decimal import Decimal as D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c16_trunc_check as T  # noqa: E402
from c16_check import PI, SLACK  # noqa: E402

E = 100
P, Q = 5, 41          # mu = 5/41 in c16_trunc_check's form; c_k = 41/(n (41 + 5 n)) = (41/5) / (n (n + 41/5))
PRINTED = D("2.156225504")


def c_sums(binary, R, p, q):
    out = subprocess.run([binary, str(R), str(p), str(q), str(E)], capture_output=True, text=True, check=True)
    rows = [tuple(int(x) for x in line.split()) for line in out.stdout.split("\n") if line.strip()]
    assert len(rows) == 3
    return rows, out.stderr.strip()


def bound(enc, scale):
    (a0, b0), (a1, b1), (a2, b2) = [(D(lo) / D(scale), D(hi) / D(scale)) for lo, hi in enc]
    dlo, dhi = a2 / b1, b2 / a1
    B_lo = a0 ** 2 / (PI * b1) - dhi.ln() - (1 + dhi.ln()).ln() - SLACK
    B_hi = b0 ** 2 / (PI * a1) - dlo.ln() - (1 + dlo.ln()).ln() + SLACK
    return B_lo, B_hi, ((a0, b0), (a1, b1), (a2, b2)), (dlo, dhi)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    binary = sys.argv[1]
    controls = "--controls" in sys.argv[2:]
    ok = True
    t0 = time.time()
    T.K = 1 << E           # the Python reference at the same scale
    for R in (300, 1200):
        py = T.sums(R, P, Q)
        cc, _ = c_sums(binary, R, P, Q)
        same = [tuple(x) for x in py] == cc
        print(f"self-test R = {R}: C sums equal the Python sums at the scale 2^{E}, integer for integer: {same}")
        ok &= same
    R = 40000
    t1 = time.time()
    enc, info = c_sums(binary, R, P, Q)
    B_lo, B_hi, S, (dlo, dhi) = bound(enc, 1 << E)
    print(f"[Pat2026b]: c_k = 1/(|k|^2 (|k|^2 + 41/5)) on 0 < |k| <= {R}; {info}; {time.time() - t1:.0f} s")
    for name, (a, b) in zip(("S0", "S1", "S2"), S):
        print(f"  {name} (times 41/5 scaling) in [{a:.30f}, {b:.30f}]")
    print(f"  delta = S2/S1 in [{dlo:.20f}, {dhi:.20f}]")
    print(f"  L >= {B_lo:.20f}   (enclosure [{B_lo:.20f}, {B_hi:.20f}])")
    good = B_lo > PRINTED
    print(f"  L >= {PRINTED}, the printed bound: {'PASS' if good else 'FAIL'} (margin {B_lo - PRINTED:.3e}); the row's "
          f"2.15622550442903... lies in our enclosure: {D('2.15622550442903') <= B_hi and D('2.15622550442904') >= B_lo}")
    ok &= good
    # the K = 4000 row [Pat2026] of PR #226, by the same C path: our earlier replay printed 2.15622475379391336305
    enc4, _ = c_sums(binary, 4000, P, Q)
    b4 = bound(enc4, 1 << E)[0]
    same4 = abs(b4 - D("2.15622475379391336305")) < D("1e-19")
    print(f"  cross-check, R = 4000 ([Pat2026], our #226 replay 2.15622475379391336305): L >= {b4:.20f}: {same4}")
    ok &= same4
    if controls:
        print("Red controls (each must be refused):")
        up = B_hi + D("1e-12")
        r1 = not (B_lo > up)
        print(f"  the printed bound raised to {up:.15f} (just above our enclosure): refused: {r1}")
        ok &= r1
        encc, _ = c_sums(binary, 20000, P, Q)
        bc = bound(encc, 1 << E)[0]
        r2 = not (bc > PRINTED)
        print(f"  the cutoff halved (R = 20000): L >= {bc:.15f}, below {PRINTED}: {r2}")
        ok &= r2
        encm, _ = c_sums(binary, R, 6, 41)
        bm = bound(encm, 1 << E)[0]
        r3 = not (bm > PRINTED)
        print(f"  mu = 6/41 instead of 5/41 (R = 40000): L >= {bm:.15f}, below {PRINTED}: {r3}")
        ok &= r3
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
