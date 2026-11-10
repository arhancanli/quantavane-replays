"""Modal run (lane tao-c1), robust version: 16 independent one-core calls (.map), each its own input with its own
timeout, each writing its best certificate to a modal.Volume after every improving stage, so a lost input or a dropped
client loses at most one worker's unsaved stage.  Problems: C1a and C1d upper-bound step functions of
teorth/optimizationproblems, refined from vidigoat's 10 Nov certificates (#270; github.com/vidigoat/optimization-records
@ a559680e, published data only) by L-BFGS(-B) on log ||v*v||_p - 2 log S with p doubled from a restart value p0.
Exactness is checked locally afterwards by our two exact checkers.  No secrets in the image: numpy, scipy, this file and
the two public certificates.

  ~/kaggle/modal-venv/bin/modal run c1ad_retake2.py       (after 16 modal_budget.py reservations of 1 core, 1 GiB, 3600 s)
"""
import json
import os
import time

for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_k] = "1"

CPU, MEM_GIB, TIMEOUT_S = 1, 1, 3600          # each call; must equal each reservation
WORK_S = 3300
P0S = [1e8, 5e7, 2e8, 3e7, 3e8, 1.5e8, 7e7, 1e9]
HERE = os.path.dirname(os.path.abspath(__file__))
VOL_PATH = "/vol"


def refine(job, save):
    import numpy as np
    from scipy.optimize import minimize
    from scipy.fft import rfft, irfft, next_fast_len
    t0 = time.time()
    kind, p0, pmax, iters = job["kind"], job["p0"], job["pmax"], job["iters"]
    NN = kind == "1a"
    v = np.array([float(x) for x in json.load(open(job["warm"]))["v_int"]])
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
        a = np.abs(c)
        M = a.max()
        z = a / M
        zp = z ** (p - 1)
        s = (zp * z).sum()
        f = np.log(M) + np.log(s) / p - 2 * np.log(S)
        gc = np.sign(c) * zp / (M * s)
        g = 2 * irfft(rfft(gc, L) * np.conj(X), L)[:n] - 2 / S
        return f, g

    def norm(x):
        if NN:
            x = np.maximum(x, 0.0)
        return x / x.sum() * (2 * n)

    x0 = norm(v)
    start = value(x0)
    best = {"val": start, "x": x0}
    log = [f"start {kind} n={n} p0={p0:g} value {start:.17f}"]
    chain = 0
    while time.time() - t0 < job["budget"] - 240:
        q0 = p0 * (3.0 ** (chain % 3))
        chain += 1
        x = best["x"].copy()
        p = q0
        while p <= pmax and time.time() - t0 < job["budget"] - 90:
            res = minimize(lambda z: fg(z, p), x, jac=True, method="L-BFGS-B",
                           bounds=[(0.0, None)] * n if NN else None,
                           options={"maxiter": iters, "maxcor": 20, "ftol": 1e-15, "gtol": 1e-300,
                                    "maxfun": int(1.3 * iters) + 10})
            x = norm(res.x)
            val = value(x)
            log.append(f"chain {chain} p={p:.2e} it={res.nit} value {val:.17f} best {min(val, best['val']):.17f} "
                       f"t={time.time() - t0:.0f}s")
            if val < best["val"]:
                best = {"val": val, "x": x.copy()}
                xi = np.rint(x * (2.0 ** 56 / np.abs(x).max()))
                if NN:
                    xi = np.maximum(xi, 0)
                save({"job": job["name"], "kind": kind, "n": n, "start": start, "best": val,
                      "best_int_float": value(xi), "v_int": [int(t) for t in xi], "log": log})
            p *= 2
    return {"job": job["name"], "kind": kind, "start": start, "best": best["val"], "log": log}


def jobs(data_dir, budget=WORK_S, iters=3000):
    J = []
    for kind, warm, pmax in (("1d", "vidigoat_c1d_n163840.json", 3.3e12), ("1a", "vidigoat_c1a_n262144_v3.json", 1.1e12)):
        for i, p0 in enumerate(P0S):
            J.append({"name": f"{kind}_{i}", "kind": kind, "p0": p0, "pmax": pmax, "iters": iters, "budget": budget,
                      "warm": os.path.join(data_dir, warm)})
    return J


try:
    import modal
except ImportError:
    modal = None

if modal is not None:
    image = (modal.Image.debian_slim(python_version="3.12")
             .pip_install("numpy==2.1.3", "scipy==1.14.1")
             .add_local_file(os.path.join(HERE, "data", "vidigoat_c1d_n163840.json"), "/data/vidigoat_c1d_n163840.json")
             .add_local_file(os.path.join(HERE, "data", "vidigoat_c1a_n262144_v3.json"), "/data/vidigoat_c1a_n262144_v3.json"))
    vol = modal.Volume.from_name("qv-tao-c1", create_if_missing=True)
    app = modal.App("qv-tao-c1-retake2", image=image)

    @app.function(cpu=CPU, memory=MEM_GIB * 1024, timeout=TIMEOUT_S, volumes={VOL_PATH: vol})
    def worker(job):
        def save(rec):
            path = os.path.join(VOL_PATH, f"best_{rec['job']}.json")
            with open(path + ".tmp", "w") as fh:
                json.dump(rec, fh)
            os.replace(path + ".tmp", path)
            vol.commit()
        return refine(job, save)

    @app.local_entrypoint()
    def main():
        t = time.time()
        os.makedirs(os.path.join(HERE, "out2"), exist_ok=True)
        for r in worker.map(jobs("/data"), return_exceptions=True):
            if isinstance(r, Exception):
                print("worker failed:", repr(r)[:200], flush=True)
                continue
            json.dump(r, open(os.path.join(HERE, "out2", f"summary_{r['job']}.json"), "w"))
            print(r["job"], f"{r['start']:.15f} -> {r['best']:.15f}", f"wall {int(time.time() - t)} s", flush=True)
        print("all calls returned, wall", int(time.time() - t), "s")
