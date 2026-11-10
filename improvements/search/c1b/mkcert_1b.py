"""Round a float C_1b step function (heights h in [0,1] on n cells of [0,2]) to integers a_i in [0, D] with
sum a = n D / 2 exactly, and write a certificate JSON (the exact fraction is added after the exact checks)."""
import json, sys
import numpy as np
from fractions import Fraction

src, out = sys.argv[1], sys.argv[2]
D = int(sys.argv[3]) if len(sys.argv) > 3 else 1 << 48
d = json.load(open(src))
h = np.array(d['h'], dtype=np.float64)
n = len(h)
assert n % 2 == 0
x = np.clip(h, 0.0, 1.0) * D
a = np.rint(x).astype(object)
a = [int(min(max(v, 0), D)) for v in a]
need = n * D // 2 - sum(a)
resid = x - np.array([float(v) for v in a])          # positive: rounded down
if need > 0:
    order = np.argsort(-resid)                        # raise the cells rounded down the most
    for i in order:
        if need == 0:
            break
        if a[i] < D:
            a[i] += 1; need -= 1
elif need < 0:
    order = np.argsort(resid)
    for i in order:
        if need == 0:
            break
        if a[i] > 0:
            a[i] -= 1; need += 1
assert need == 0 and 2 * sum(a) == n * D and min(a) >= 0 and max(a) <= D
hf = np.array(a, dtype=np.float64) / D
from scipy.fft import rfft, irfft, next_fast_len
L = next_fast_len(2 * n, real=True)
c = irfft(np.conj(rfft(hf, L)) * rfft(1.0 - hf, L), L)
v = 2.0 * max(c[:n].max(), c[L - n + 1:].max()) / n
print(f'n={n} D={D} rounded float value {v:.16f} (source {d.get("float_score")})')
json.dump({'constant': 'C_1b', 'bound': 'upper',
           'definition': 'f = a_i / D on the n equal cells of [-1,1] (i = 0..n-1), g = 1 - f on [-1,1]; '
                         'C_1b <= (2/n) max_{|k|<n} sum_i a_i (D - a_{i+k}) / D^2',
           'n': n, 'D': D, 'a': a}, open(out, 'w'), separators=(',', ':'))
