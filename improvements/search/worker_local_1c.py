#!/usr/bin/env python3
'''Worker for the Kaggle GPU job: refine step-function certificates for C_1c (lower bound) and C_1d (upper bound).

Smooth surrogates minimised by L-BFGS(-B) with float64 FFTs on the GPU:
  1c: soft-min_beta of r_k = K A_k / S^2 over 0 <= k <= K (w >= 0), beta raised in steps;
  1d: log ||v*v||_p - 2 log S, p doubled from p0 up to 1e11, restarted from the best vector.
Warm starts: Julsgaard's certificates (teorth/optimizationproblems #318 for 1c, #312 for 1d), from the dataset.
Every improvement of the float64 score is rounded to integers and written to /kaggle/working/best_<job>.json.
'''
import glob
import json
import os
import subprocess
import sys
import time

TOTAL = float(os.environ.get("C1_TOTAL", "9300"))
SMOKE = os.environ.get("C1_SMOKE") == "1"
IT = (lambda k: max(5, k // 100)) if SMOKE else (lambda k: k)
JOBS = {
    "Lc1": {"kind": "1c", "gpu": 0, "warm": "best_c1_n98368.json", "factor": 1,
            "cycle": [1e8, 3e8, 1e9, 3e9, 1e10, 3e7], "stage_iters": 20000, "chunk": 2500},
}
OUT = "/kaggle/working" if os.path.isdir("/kaggle/working") else os.path.abspath("out")
os.makedirs(OUT, exist_ok=True)


def find(name):
    hits = glob.glob("/kaggle/input/**/" + name, recursive=True) + glob.glob(name)
    if not hits:
        raise SystemExit("missing input " + name)
    return hits[0]


def smooth_len(N):
    best, a = None, 1
    while a < 4 * N:
        b = a
        while b < 4 * N:
            c = b
            while c < 4 * N:
                if c >= N and (best is None or c < best):
                    best = c
                c *= 5
            b *= 3
        a *= 2
    return best


def worker(name):
    import numpy as np
    import torch
    from scipy.optimize import minimize
    job = JOBS[name]
    kind = job["kind"]
    ng = torch.cuda.device_count()
    dev = torch.device(f"cuda:{job['gpu'] % ng}" if ng else "cpu")
    if dev.type == "cuda":
        torch.cuda.set_device(dev)
    t_start = time.time()
    logf = open(os.path.join(OUT, f"log_{name}.txt"), "a")

    def log(s):
        line = f"[{time.time() - t_start:7.0f}s] {s}"
        print(f"{name} {line}", flush=True)
        logf.write(line + "\n")
        logf.flush()

    def left():
        return TOTAL - (time.time() - t_start)

    log(f"device {dev} {torch.cuda.get_device_name(dev) if dev.type == 'cuda' else 'cpu'} torch {torch.__version__}")

    class P:
        def __init__(self, n, K):
            self.n, self.K = n, K
            self.L = smooth_len(2 * n - 1)

        def value(self, x):
            X = torch.fft.rfft(x, self.L)
            S = x.sum()
            if kind == "1d":
                c = torch.fft.irfft(X * X, self.L)[: 2 * self.n - 1]
                return float(2 * self.n * c.abs().max() / S ** 2)
            A = torch.fft.irfft(X * X.conj(), self.L)[: self.K + 1]
            return float(self.K * A.min() / S ** 2)

        def fg(self, x, t):
            n, L = self.n, self.L
            X = torch.fft.rfft(x, L)
            S = x.sum()
            if kind == "1d":
                p = t
                c = torch.fft.irfft(X * X, L)[: 2 * n - 1]
                a = c.abs()
                M = a.max()
                z = a / M
                s = (z ** p).sum()
                f = torch.log(M) + torch.log(s) / p - 2 * torch.log(S)
                gc = torch.sign(c) * z ** (p - 1) / (M * s)
                g = 2 * torch.fft.irfft(torch.fft.rfft(gc, L) * X.conj(), L)[:n] - 2 / S
                return float(f), g
            beta, K = t, self.K
            A = torch.fft.irfft(X * X.conj(), L)[: K + 1]
            r = K * A / S ** 2
            m = r.min()
            e = torch.exp(-beta * (r - m))
            Z = e.sum()
            F = m - torch.log(Z) / beta
            pk = e / Z
            Pf = torch.fft.rfft(pk, L)
            s1 = torch.fft.irfft(X * Pf.conj(), L)[:n] + torch.fft.irfft(Pf * X, L)[:n]
            g = (K / S ** 2) * s1 - (2 / S) * (pk * r).sum()
            return -float(F), -g

    better = (lambda a, b: a > b) if kind == "1c" else (lambda a, b: a < b)
    st = {"best": None, "x": None, "prob": None}

    def norm(prob, x):
        return x / x.sum() * (prob.K if kind == "1c" else 2 * prob.n)

    def consider(prob, x, tag):
        val = prob.value(x)
        if st["best"] is None or st["prob"].n != prob.n or better(val, st["best"]):
            xa = x.detach().cpu().numpy().astype(np.float64)
            xi = np.rint(xa * (2.0 ** 56 / np.abs(xa).max()))
            if kind == "1c":
                xi = np.maximum(xi, 0)
            vi = prob.value(torch.from_numpy(xi).to(dev))
            st["best"], st["x"], st["prob"] = val, x.clone(), prob
            ints = [int(v) for v in xi]
            if kind == "1c":
                rec = {"problem": "C_1c lower bound witness: f = sum_i w_i 1_[i/K,(i+1)/K), "
                                  "bound K min_{0<=k<=K} A_k / S^2", "K": prob.K, "weights_int": ints}
            else:
                rec = {"problem": "C_1d upper bound witness: n equal cells on [-1/4,1/4], values v_int, "
                                  "bound 2n max|v*v| / S^2", "n": prob.n, "v_int": ints}
            rec.update({"float_score": vi, "tag": tag, "job": name})
            path = os.path.join(OUT, f"best_{name}_n{prob.n}.json")
            with open(path + ".tmp", "w") as fh:
                json.dump(rec, fh)
            os.replace(path + ".tmp", path)
            log(f"BEST {tag}: {val:.17f} (integer {vi:.17f}) n={prob.n}")
        return val

    def lbfgs(prob, x, t, maxiter, nonneg, ftol=1e-16, m=20):
        def fg(z):
            f, g = prob.fg(torch.from_numpy(z).to(dev), t)
            return f, g.detach().cpu().numpy()
        bounds = [(0, None)] * prob.n if nonneg else None
        res = minimize(fg, x.detach().cpu().numpy(), jac=True, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": maxiter, "maxcor": m, "ftol": ftol, "gtol": 1e-300,
                                "maxfun": int(1.3 * maxiter) + 10})
        return norm(prob, torch.from_numpy(res.x).to(dev)), res

    if kind == "1c":
        d = json.load(open(find(job["warm"])))
        x0 = np.array(d["weights_int"], dtype=np.float64)
        prob = P(len(x0), int(d["K"]))
        x = norm(prob, torch.from_numpy(x0).to(dev))
        log(f"start {job['warm']} K={prob.K} n={prob.n} L={prob.L} value {prob.value(x):.17f}")
        if job["factor"] > 1:
            consider(prob, x, "start")
            xs = np.repeat(x.detach().cpu().numpy(), job["factor"])
            prob = P(prob.n * job["factor"], prob.K * job["factor"])
            x = norm(prob, torch.from_numpy(xs).to(dev))
            log(f"split x{job['factor']}: K={prob.K} n={prob.n} L={prob.L} value {prob.value(x):.17f}")
        consider(prob, x, "start")
        cyc = 0
        while left() > 60:
            cyc += 1
            x = st["x"].clone()
            for beta in job["cycle"]:
                its, t0, v_start = 0, time.time(), prob.value(x)
                while its < IT(job["stage_iters"]) and left() > 30:
                    x, res = lbfgs(prob, x, beta, IT(job["chunk"]), True)
                    its += res.nit
                    v1 = consider(prob, x, f"cycle {cyc} beta={beta:.0e}")
                    if res.nit < 0.5 * IT(job["chunk"]):
                        break
                log(f"cycle {cyc} beta={beta:.1e} its={its} value {prob.value(x):.17f} (stage gain "
                    f"{prob.value(x) - v_start:.3e}) best {st['best']:.17f} {time.time() - t0:.0f}s")
                if left() <= 30:
                    break
        log(f"done 1c best {st['best']:.17f} K={st['prob'].K}")
        return

    d = json.load(open(find("julsgaard_certificate.json")))
    x0 = np.array(d["v_int"], dtype=np.float64)
    prob = P(len(x0) * job["factor"], None)
    x = norm(prob, torch.from_numpy(np.repeat(x0, job["factor"])).to(dev))
    log(f"start n={prob.n} (factor {job['factor']}) L={prob.L} value {prob.value(x):.17f}")
    consider(prob, x, "start")
    if job["polish"]:
        end = time.time() + 0.15 * left()
        for p in [1e8, 1e9, 1e10, 1e11]:
            t0 = time.time()
            v0 = prob.value(x)
            x, res = lbfgs(prob, x, p, IT(3000), False, ftol=1e-15)
            v1 = consider(prob, x, f"polish p={p:.0e}")
            log(f"polish p={p:.0e} it={res.nit} value {v1:.17f} gain {v0 - v1:.3e} {time.time() - t0:.0f}s")
            if time.time() > end:
                break
    ri = 0
    while left() > 240:
        p0 = job["restarts"][ri % len(job["restarts"])]
        ri += 1
        xr = st["x"].clone()
        p, t0, its = p0, time.time(), 0
        end = time.time() + left() - 60
        while p <= 1.01e11 and time.time() < end:
            xr, res = lbfgs(prob, xr, p, IT(4000), False, ftol=1e-14)
            its += res.nit
            p *= 2
        v = consider(prob, xr, f"restart p0={p0:g}")
        log(f"restart p0={p0:g} reached p={p / 2:.1e} its={its} value {v:.17f} best {st['best']:.17f} "
            f"{time.time() - t0:.0f}s")
    log(f"done 1d best {st['best']:.17f}")


if __name__ == "__main__":
    worker(sys.argv[1])
