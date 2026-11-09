"""CPU continuation for the C1d upper bound (numpy/scipy): p-norm continuation restarts from the current best
integer certificate, L-BFGS on log ||v*v||_p - 2 log S; every improvement of the float score is rounded to integers
(scale 2^56) and written as best_c1d.json. Exactness is checked separately by improvements/c1d/check.py."""
import json, os, sys, time
import numpy as np
from scipy.optimize import minimize
from scipy.fft import rfft, irfft, next_fast_len

warm, out, total = sys.argv[1], sys.argv[2], float(sys.argv[3])
restarts = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [8, 128, 24, 48, 12, 96]
t0 = time.time()
d = json.load(open(warm))
v = np.array([float(x) for x in d["v_int"]])
n = len(v)
L = next_fast_len(2 * n - 1, real=True)
logf = open(os.path.join(out, "log.txt"), "a")

def log(s):
    line = f"[{time.time() - t0:7.0f}s] {s}"
    print(line, flush=True); logf.write(line + "\n"); logf.flush()

def value(x):
    X = rfft(x, L)
    c = irfft(X * X, L)[: 2 * n - 1]
    return 2 * n * np.abs(c).max() / x.sum() ** 2

def fg(x, p):
    X = rfft(x, L)
    S = x.sum()
    c = irfft(X * X, L)[: 2 * n - 1]
    a = np.abs(c)
    M = a.max()
    z = a / M
    zp = z ** (p - 1)
    s = (zp * z).sum()
    f = np.log(M) + np.log(s) / p - 2 * np.log(S)
    gc = np.sign(c) * zp / (M * s)
    g = 2 * irfft(rfft(gc, L) * np.conj(X), L)[:n] - 2 / S   # d c_m/d v_i = 2 v_{m-i}: correlation of gc with v
    return f, g

best = {"val": value(v / v.sum()), "x": v / v.sum() * 2 * n}
log(f"start n={n} L={L} value {best['val']:.17f}")

def consider(x, tag):
    val = value(x)
    if val < best["val"]:
        xi = np.rint(x * (2.0 ** 56 / np.abs(x).max()))
        vi = value(xi)
        best["val"], best["x"] = val, x.copy()
        rec = {"problem": "C_1d upper bound witness: n equal cells on [-1/4,1/4], values v_int, bound 2n max|v*v| / S^2",
               "n": n, "v_int": [int(t) for t in xi], "float_score": vi, "tag": tag}
        p = os.path.join(out, "best_c1d.json")
        with open(p + ".tmp", "w") as fh:
            json.dump(rec, fh)
        os.replace(p + ".tmp", p)
        log(f"BEST {tag}: {val:.17f} (integer {vi:.17f})")
    return val

ri = 0
while time.time() - t0 < total - 600:
    p0 = restarts[ri % len(restarts)]; ri += 1
    x = best["x"].copy(); p = p0; its = 0; ts = time.time()
    while p <= 1.01e11 and time.time() - t0 < total - 300:
        res = minimize(lambda z: fg(z, p), x, jac=True, method="L-BFGS-B",
                       options={"maxiter": 4000, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-300, "maxfun": 5300})
        x = res.x / res.x.sum() * 2 * n
        its += res.nit
        p *= 2
    v1 = consider(x, f"restart p0={p0:g}")
    log(f"restart p0={p0:g} reached p={p / 2:.1e} its={its} value {v1:.17f} best {best['val']:.17f} {time.time() - ts:.0f}s")
log(f"done best {best['val']:.17f}")
