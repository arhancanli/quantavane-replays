import sys, time, numpy as np
from fp import solve, weights, value
from opt import best_M
Ks = [int(x) for x in sys.argv[1].split(",")]
M16 = 80.205888
for K in Ks:
    c = M16 * (K / 16384) ** 0.5
    t = time.time(); M, v, g = best_M(K, 1, 0.5, c * 0.995, c * 1.005, tol=2e-7)
    nz = np.nonzero(g)[0]
    print(f"K={K}: M={M:.6f} (M/sqrt(K) {M/K**0.5:.6f}) value {v:.13f}  law {0.410378-0.205/K:.13f}  (v + 0.205/K = {v+0.205/K:.10f}) j1/K {nz[0]/K:.5f} {time.time()-t:.0f}s", flush=True)
    np.save(f"g_K{K}_m1.npy", g); open(f"M_K{K}_m1.txt", "w").write(repr(M))
