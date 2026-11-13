"""Program B (written separately): the same quantity by 17-bit limbs and exact int64 correlations in numpy.
Limb products are below 2^34 and at most 2^17 terms are summed, so every partial sum is below 2^51 (no wraparound).
Usage: c1c366_limbs.py FILE [--bump J]"""
import json, sys
import numpy as np
from fractions import Fraction
d = json.load(open(sys.argv[1])); K = d["K"]; w = [int(x) for x in d["weights_int"]]
if "--bump" in sys.argv: w[int(sys.argv[sys.argv.index("--bump") + 1])] += 1
n = len(w); assert max(w) < 2 ** 51 and n <= 2 ** 17
L = [np.array([(x >> (17 * a)) & 0x1FFFF for x in w], dtype=np.int64) for a in range(3)]
tot = [0] * (K + 1)
for a in range(3):
    for b in range(3):
        c = np.correlate(L[b], L[a], mode="full")[n - 1:n + K]   # c[k] = sum_i L[a][i] * L[b][i+k], k = 0..K
        assert c.dtype == np.int64 and int(c.max()) < 2 ** 51
        sh = 17 * (a + b)
        for k in range(K + 1): tot[k] += int(c[k]) << sh
S = sum(w); m = min(tot); ks = [k for k, v in enumerate(tot) if v == m]
f = Fraction(K * m, S * S)
print("B: n", n, "K", K, "S", S, "min A_k", m, "at k", ks)
print("B: bound", f.numerator, "/", f.denominator, "=", float(f))
print("B: equals stated fraction:", f == Fraction(d["score_num"], d["score_den"]))
