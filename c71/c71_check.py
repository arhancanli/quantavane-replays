#!/usr/bin/env python3
"""Independent replay of the three explicit lower bounds for the Fourier Entropy-Influence constant C_71 on the 71a page
([MI2026] n = 14, [MI2026b] n = 17, [Num2026] n = 18), from the functions alone. No program of any certificate was
read or run; truth tables and spectra are read as data.

For f: {-1,1}^n -> {-1,1} (index x in [0, 2^n), bit i of x set for x_i = -1), F_S = sum_x f(x) (-1)^|S & x| (fast
Walsh-Hadamard transform, integers), fhat(S)^2 = F_S^2 / 4^n, the spectral entropy H = sum_S fhat(S)^2 log2(1/fhat(S)^2)
(all S, as on the page) and the total influence I = sum_S |S| fhat(S)^2. The page's rule (O'Donnell-Tan, arXiv:1304.1347,
Lemma 5.1; for balanced f with H > 0): C_71 >= H/(I - 1). We check exactly: every value is +1 or -1 (for [MI2026b],
whose certificate gives the spectrum and not the table, the table is rebuilt by the inverse transform and must come
out in {+1, -1}), Parseval, balance (F_empty = 0), I as an exact rational, and monotonicity (f or -f non-decreasing
in every coordinate, over every comparable pair of neighbours); H is enclosed in 60-digit decimal (logarithms
correctly rounded, the sum widened by 1e-50).

  python3 -I c71_check.py ROW [CERTS] [--controls]       rows: mi2026, mi2026b, num2026
Python 3.9+, standard library only."""
import json
import os
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction as Fr

getcontext().prec = 60
WIDEN = D("1e-50")
HERE = os.path.dirname(os.path.abspath(__file__))
LN2 = D(2).ln()

ROWS = {
    "mi2026": dict(n=14, claim="6.4901128435233943", I=Fr(521, 256), file=("mi2026", "c71_input_table_n14.json"),
                   where="71a row [MI2026] (PR #94): g14 on 14 variables, the gist's c71_input_table_n14.json"),
    "mi2026b": dict(n=17, claim="6.514326913930565372", I=Fr(261, 128),
                    file=("mi2026b", "ridge_n17_g4_best2_iter14.json"),
                    where="71a row [MI2026b] (PR #124): 17 variables, the spectrum in ridge_n17_g4_best2_iter14.json"),
    "num2026": dict(n=18, claim="6.521845710923046575", I=Fr(261, 128), file=("num2026", "fei_c71_n18_artifact.json"),
                    where="71a row [Num2026] (PR #130): 18 variables, true_hex of fei_c71_n18_artifact.json"),
}


def wht(v):
    """in-place fast Walsh-Hadamard transform of a list of integers"""
    n = len(v)
    h = 1
    while h < n:
        for i in range(0, n, 2 * h):
            for j in range(i, i + h):
                a, b = v[j], v[j + h]
                v[j], v[j + h] = a + b, a - b
        h *= 2
    return v


def load(key, certs):
    row = ROWS[key]
    path = os.path.join(certs, *row["file"])
    d = json.load(open(path, encoding="utf-8"))
    n = row["n"]
    if key == "mi2026":
        f = list(d["truth_table"])
        assert len(f) == 1 << n
        return f, None
    if key == "num2026":
        assert d["n"] == n
        h = d["true_hex"]
        bits = bin(int(h, 16))[2:].zfill(4 * len(h))
        assert len(bits) == 1 << n
        f = [-1 if b == "1" else 1 for b in bits]        # TRUE -> -1; the spectrum's magnitudes do not depend on it
        return f, None
    # mi2026b: the integer spectrum F_S, S as a bit mask; rebuild the table
    assert d["n"] == n
    F = [0] * (1 << n)
    for k, v in d["spectrum_nonzero"].items():
        F[int(k)] = int(v)
    g = wht(list(F))                                   # inverse transform up to the factor 2^n
    N = 1 << n
    assert all(x % N == 0 for x in g)
    f = [x // N for x in g]
    return f, d


def analyse(f, n):
    N = 1 << n
    vals = set(f)
    boolean = vals <= {1, -1}
    F = wht(list(f))
    parseval = sum(x * x for x in F) == N * N
    balanced = F[0] == 0
    I = Fr(sum(bin(S).count("1") * x * x for S, x in enumerate(F) if x), N * N)
    s = D(0)
    for x in F:
        if x:
            w = Fr(x * x, N * N)
            s += (D(w.numerator) / D(w.denominator)) * (D(w.denominator).ln() - D(w.numerator).ln())
    H = ((s - WIDEN) / LN2, (s + WIDEN) / LN2)
    nz = sum(1 for x in F if x)
    return dict(boolean=boolean, parseval=parseval, balanced=balanced, I=I, H=H, nonzero=nz, F=F)


def monotone(f, n):
    """f or -f non-decreasing in every coordinate (bit set = x_i = -1); returns (+1, -1 or 0)"""
    up = down = True
    for x in range(1 << n):
        fx = f[x]
        for i in range(n):
            if not x >> i & 1:
                fy = f[x | 1 << i]
                if fy < fx:
                    up = False
                elif fy > fx:
                    down = False
        if not up and not down:
            return 0
    return 1 if up else (-1 if down else 0)


def ratio(A):
    I1 = A["I"] - 1
    lo = A["H"][0] / (D(I1.numerator) / D(I1.denominator))
    hi = A["H"][1] / (D(I1.numerator) / D(I1.denominator))
    return lo, hi


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args or args[0] not in ROWS:
        print(__doc__)
        return 2
    key = args[0]
    certs = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    controls = "--controls" in sys.argv
    row = ROWS[key]
    n = row["n"]
    t0 = time.time()
    print(row["where"] + f": C_71 >= {row['claim']}")
    f, extra = load(key, certs)
    ok = True
    if extra is not None:
        print(f"  the inverse transform of the certificate's {len(extra['spectrum_nonzero'])} nonzero coefficients is a "
              f"table of 2^{n} values in {{+1, -1}}: {set(f) <= {1, -1}}")
        ok &= set(f) <= {1, -1}
    A = analyse(f, n)
    good = A["boolean"] and A["parseval"] and A["balanced"]
    print(f"  Boolean: {A['boolean']}; Parseval sum F_S^2 = 4^n: {A['parseval']}; balanced (F_empty = 0): "
          f"{A['balanced']}; nonzero coefficients: {A['nonzero']}")
    ok &= good
    good = A["I"] == row["I"]
    print(f"  I = {A['I']} exactly (the row: {row['I']}): {good}")
    ok &= good
    if extra is not None:
        same = all(int(v) == A["F"][int(k)] for k, v in extra["spectrum_nonzero"].items()) and \
            sum(1 for x in A["F"] if x) == len(extra["spectrum_nonzero"])
        print(f"  the forward transform of the rebuilt table returns the certificate's spectrum: {same}")
        ok &= same
    mono = monotone(f, n)
    print(f"  monotone ({'f' if mono == 1 else '-f' if mono == -1 else 'neither f nor -f'} non-decreasing in every "
          f"coordinate): {mono != 0}")
    ok &= mono != 0
    print(f"  H = {A['H'][0]:.40f}... (bits)")
    lo, hi = ratio(A)
    print(f"  H/(I - 1) in [{lo:.40f}, {hi:.40f}]")
    claim = D(row["claim"])
    good = lo > claim
    print(f"  the row's {row['claim']} < H/(I - 1): {good} (margin {float(lo - claim):.3e})")
    ok &= good
    if key == "num2026" and os.path.exists(os.path.join(certs, *ROWS["mi2026b"]["file"])):
        f17, _ = load("mi2026b", certs)
        B = analyse(f17, 17)
        dH = A["H"][0] - B["H"][0]
        good = abs(dH - D(1) / 128) < D("1e-45") and B["I"] == A["I"]
        lo17, _ = ratio(B)
        print(f"  against the [MI2026b] function: the same I, H larger by {dH:.40f} = 1/128 (to 1e-45), so the "
              f"gain is (1/128)/(133/128) = 1/133, as the row states: {good}")
        ok &= good
    if controls:
        print("Red controls (each must not certify the row's value):")
        g = list(f)
        i = next(x for x in range(1 << n) if g[x] == 1)
        g[i] = -1
        B = analyse(g, n)
        print(f"  one value of the table flipped: balanced {B['balanced']} (the rule needs balance): "
              f"{'refused, as expected' if not B['balanced'] else 'NOT REFUSED'}")
        ok &= not B["balanced"]
        j = next(x for x in range(1 << n) if g[x] == 1 and x != i)
        g = list(f)
        g[i], g[j] = -1, -1
        k1 = next(x for x in range(1 << n) if f[x] == -1)
        k2 = next(x for x in range(k1 + 1, 1 << n) if f[x] == -1)
        g[k1], g[k2] = 1, 1                           # a balanced change of four values
        B = analyse(g, n)
        lo2, hi2 = ratio(B) if B["I"] > 1 else (D(0), D(0))
        good = B["balanced"] and hi2 < claim
        print(f"  two +1 and two -1 values swapped (still balanced): H/(I - 1) <= {hi2:.12f}: "
              f"{'below, as expected' if good else 'NOT BELOW'}")
        ok &= good
    print(f"RESULT: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
