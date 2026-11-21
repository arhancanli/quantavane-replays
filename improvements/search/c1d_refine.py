"""Local refinement of a C_1d upper-bound step function (search only; exactness is checked afterwards by the two exact
checkers in quantavane-replays/improvements/c1d: check.py and ntt_check_conv.py).

The objective is the lab's own (quantavane-replays/improvements/search/c1ad_modal/c1ad_retake2.py, 10 Nov 2026):
L-BFGS on log ||v*v||_p - 2 log S with p doubled from a restart value p0, best vector kept.  Warm start: C. Weill's
public certificate of teorth/optimizationproblems PR #374 (cweill/c1d-upper-bound @ 7889e93, data only; his programs
were not read or run).

    python -u c1d_refine.py WARM.json OUT_PREFIX --up 1 --p0 3e7 --pmax 3e12 --iters 3000 --budget 3600
"""
import argparse
import json
import os
import time

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.optimize import minimize


def wait_load(limit=60.0):
    while os.getloadavg()[0] >= limit:
        time.sleep(30)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("warm"); ap.add_argument("out")
    ap.add_argument("--up", type=int, default=1)
    ap.add_argument("--p0", type=float, default=3e7)
    ap.add_argument("--pmax", type=float, default=3.3e12)
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--budget", type=float, default=3600)
    ap.add_argument("--nonneg", action="store_true", help="C_1a: cell values >= 0 (L-BFGS-B bounds)")
    a = ap.parse_args()
    t0 = time.time()
    d = json.load(open(a.warm))
    v = np.array([float(x) for x in (d.get("values") or d.get("v_int"))])
    if a.up > 1:
        v = np.repeat(v, a.up)
    n = len(v)
    L = next_fast_len(2 * n - 1, real=True)

    def value(x):
        X = rfft(x, L)
        c = irfft(X * X, L)[: 2 * n - 1]
        return 2 * n * np.abs(c).max() / x.sum() ** 2

    def fg(x, p):
        X = rfft(x, L)
        S = x.sum()
        c = irfft(X * X, L)[: 2 * n - 1]
        aa = np.abs(c)
        M = aa.max()
        z = aa / M
        zp = z ** (p - 1)
        s = (zp * z).sum()
        f = np.log(M) + np.log(s) / p - 2 * np.log(S)
        gc = np.sign(c) * zp / (M * s)
        g = 2 * irfft(rfft(gc, L) * np.conj(X), L)[:n] - 2 / S
        return f, g

    x0 = v / v.sum() * (2 * n)
    start = value(x0)
    best = {"val": start, "x": x0}
    print(f"start n={n} value {start:.17f}", flush=True)
    chain = 0
    while time.time() - t0 < a.budget - 240:
        q0 = a.p0 * (3.0 ** (chain % 3))
        chain += 1
        x = best["x"].copy()
        p = q0
        while p <= a.pmax and time.time() - t0 < a.budget - 90:
            wait_load()
            res = minimize(lambda z: fg(z, p), x, jac=True, method="L-BFGS-B",
                           bounds=[(0.0, None)] * n if a.nonneg else None,
                           options={"maxiter": a.iters, "maxcor": 20, "ftol": 1e-15, "gtol": 1e-300,
                                    "maxfun": int(1.3 * a.iters) + 10})
            x = (np.maximum(res.x, 0.0) if a.nonneg else res.x)
            x = x / x.sum() * (2 * n)
            val = value(x)
            print(f"chain {chain} p={p:.2e} it={res.nit} value {val:.17f} best {min(val, best['val']):.17f} "
                  f"t={time.time() - t0:.0f}s", flush=True)
            if val < best["val"]:
                best = {"val": val, "x": x.copy()}
                xi = np.rint(x * (2.0 ** 56 / np.abs(x).max()))
                if a.nonneg:
                    xi = np.maximum(xi, 0)
                json.dump({"n": n, "start": start, "best_float": val, "best_int_float": value(xi),
                           "warm": os.path.basename(a.warm), "up": a.up,
                           "v_int": [int(t) for t in xi]}, open(a.out + "_best.json", "w"))
            p *= 2
    print(f"done start {start:.17f} best {best['val']:.17f}", flush=True)


if __name__ == "__main__":
    main()
