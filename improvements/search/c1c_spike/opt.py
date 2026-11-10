import sys, time, numpy as np
from fp import solve, weights, value
def best_M(K, m, r, lo, hi, tol=1e-6, g=None, alpha=0.3):
    cache = {}
    def f(M):
        nonlocal g
        if M not in cache:
            g, it, d = solve(K, m, r, M, g=g, iters=6000, alpha=0.5)
            v, k, rr = value(K, weights(K, m, r, M, g))
            cache[M] = (v, g.copy())
        return cache[M][0]
    gr = (5 ** 0.5 - 1) / 2
    a, b = lo, hi; c = b - gr * (b - a); d = a + gr * (b - a)
    while b - a > tol * (abs(a) + abs(b)):
        if f(c) > f(d): b, d = d, c; c = b - gr * (b - a)
        else: a, c = c, d; d = a + gr * (b - a)
    M = (a + b) / 2; v = f(M)
    return M, v, cache[M][1]
if __name__ == "__main__":
    K = int(sys.argv[1])
    for spec in sys.argv[2:]:
        m, r, lo, hi = spec.split(":"); m = int(m); r = float(r); lo = float(lo); hi = float(hi)
        t = time.time(); M, v, g = best_M(K, m, r, lo, hi)
        print(f"K={K} m={m} r={r}: best M={M:.6f} value {v:.12f} (law 0.410378-0.205/K = {0.410378-0.205/K:.9f}) {time.time()-t:.0f}s", flush=True)
