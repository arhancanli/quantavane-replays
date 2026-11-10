"""Round our solved spike-plus-profile step function to integers (largest weight 2^60) and write a certificate JSON."""
import json, sys, numpy as np
from fp import weights, value
K = int(sys.argv[1]); out = sys.argv[2]
g = np.load(f"g_K{K}_m1.npy"); M = float(open(f"M_K{K}_m1.txt").read())
w = weights(K, 1, 0.5, M, g)
scale = 2.0 ** 60 / w.max()
wi = np.rint(w * scale)
assert wi.min() >= 0 and wi.max() <= 2 ** 60
ints = [int(x) for x in wi]
while ints and ints[-1] == 0: ints.pop()
n = len(ints)
if n <= K:                      # the checkers need n > K so that every A_k, k <= K, is a full sum
    ints += [0] * (K + 1 - n)
v_float, k_float, _ = value(K, w)
v_int, k_int, _ = value(K, np.array(ints, dtype=float))
print(f"K={K} n={len(ints)} nonzero={sum(1 for x in ints if x)} max={max(ints)} float value {v_float:.15f} (k={k_float}) rounded {v_int:.15f} (k={k_int})")
json.dump({"K": K, "weights_int": ints}, open(out, "w"))
