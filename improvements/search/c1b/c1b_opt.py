"""Our own search for upper bounds on C_1b (Erdos minimum overlap), teorth/optimizationproblems constants/1b.

A step function h on n equal cells of [0,2], 0 <= h <= 1, sum(h) = n/2 (mass 1), g = 1 - h on [0,2], zero outside.
Score = (2/n) max_k C_k, C_k = sum_i h_i g_{i+k}, k in (-n, n) (the overlap is linear between grid shifts).
We minimise log ||C||_p (p-norm continuation, L-BFGS-B, box [0,1]) with a quadratic penalty on sum(h) - n/2;
the exact sum is restored afterwards. Exactness is the job of the checkers, never of this search.
"""
import json, os, sys, time, math
import numpy as np
from scipy.fft import rfft, irfft, next_fast_len
from scipy.optimize import minimize


def load_heights(path):
    d = json.load(open(path))
    if isinstance(d, dict) and 'values' in d:
        v = d['values']
        if isinstance(v[0], str):
            return np.array([int(s.split('/')[0]) / int(s.split('/')[1]) for s in v], dtype=np.float64)
        return np.array(v, dtype=np.float64)
    if isinstance(d, dict) and 'h' in d:
        return np.array(d['h'], dtype=np.float64)
    if isinstance(d, dict) and 'a' in d:
        D = int(d['D'])
        return np.array([int(x) for x in d['a']], dtype=np.float64) / D
    return np.array(d, dtype=np.float64)


def corr_all(h):
    """C_k for k = -(n-1)..(n-1), returned as an array of length 2n-1 (index k + n - 1)."""
    n = len(h)
    L = next_fast_len(2 * n, real=True)
    H = rfft(h, L)
    G = rfft(1.0 - h, L)
    c = irfft(np.conj(H) * G, L)
    return np.concatenate([c[L - n + 1:], c[:n]])


def score(h):
    n = len(h)
    return 2.0 * corr_all(h).max() / n


class Obj:
    def __init__(self, n, p, rho, sym):
        self.n = n; self.p = p; self.rho = rho; self.sym = sym
        self.L = next_fast_len(2 * n, real=True)
        self.s = n / 2.0
        self.best = None

    def full(self, y):
        return np.concatenate([y, y[::-1]]) if self.sym else y

    def __call__(self, y):
        n, L, p = self.n, self.L, self.p
        h = self.full(y)
        g = 1.0 - h
        H = rfft(h, L); G = rfft(g, L)
        c = irfft(np.conj(H) * G, L)          # c[k mod L] = C_k
        m = max(c[:n].max(), c[L - n + 1:].max())
        u = np.zeros(L)
        u[:n] = c[:n] / m
        u[L - n + 1:] = c[L - n + 1:] / m
        np.clip(u, 1e-300, None, out=u)
        lu = np.log(u)
        e = np.exp(p * lu)
        e[n:L - n + 1] = 0.0
        S = e.sum()
        F = math.log(m) + math.log(S) / p
        w = e / (S * u) / m                      # dF/dC_k
        w[n:L - n + 1] = 0.0
        W = rfft(w, L)
        a = irfft(np.conj(W) * G, L)[:n]         # sum_k w_k g_{j+k}
        b = irfft(W * H, L)[:n]                  # sum_k w_k h_{j-k}
        grad = a - b
        dev = h.sum() - self.s
        F += 0.5 * self.rho * dev * dev / n
        grad = grad + self.rho * dev / n
        if self.sym:
            half = len(y)
            grad = grad[:half] + grad[half:][::-1]
        sc = 2.0 * m / n
        if self.best is None or sc < self.best[0]:
            self.best = (sc, dev)
        return F, grad


def fix_sum(h, iters=60):
    """Restore sum(h) = n/2 by a uniform shift of the cells strictly inside (0,1), clipped."""
    n = len(h); s = n / 2.0
    h = h.copy()
    for _ in range(iters):
        dev = s - h.sum()
        if abs(dev) < 1e-13 * n:
            break
        free = (h > 0) & (h < 1)
        if free.sum() == 0:
            break
        h[free] += dev / free.sum()
        np.clip(h, 0.0, 1.0, out=h)
    return h


def wait_load(limit=48.0):
    while True:
        try:
            l1 = os.getloadavg()[0]
        except OSError:
            return
        if l1 <= limit:
            return
        time.sleep(30)


def polish(h, ps, iters, rho=1e4, sym=False, log=None, tag=''):
    n = len(h)
    if sym:
        assert n % 2 == 0
        y = 0.5 * (h[: n // 2] + h[n // 2:][::-1])
    else:
        y = h.copy()
    bestp = None
    for p in ps:
        wait_load()
        ob = Obj(n, p, rho, sym)
        t0 = time.time()
        r = minimize(ob, y, jac=True, method='L-BFGS-B', bounds=[(0.0, 1.0)] * len(y),
                     options={'maxiter': iters, 'maxcor': 30, 'ftol': 0.0, 'gtol': 0.0, 'maxfun': iters * 2})
        y = r.x
        hf = fix_sum(ob.full(y))
        sc = score(hf)
        msg = f'{tag} p={p:.3g} it={r.nit} {str(r.message)[:28]} score_fixed={sc:.13f} dev={hf.sum()-n/2:.2e} t={time.time()-t0:.0f}s'
        print(msg, flush=True)
        if log:
            with open(log, 'a') as fh:
                fh.write(msg + '\n')
        if bestp is None or sc < bestp[0]:
            bestp = (sc, hf.copy())
            if os.environ.get('C1B_SAVE'):
                save(hf, os.environ['C1B_SAVE'], note=f'{tag} p={p:.3g}')
    return bestp[1]


def upsample(h, f):
    return np.repeat(h, f)


def save(h, path, note=''):
    json.dump({'n': len(h), 'score_float': score(h), 'note': note, 'h': [float(x) for x in h]}, open(path, 'w'))


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('src'); ap.add_argument('out')
    ap.add_argument('--up', type=int, default=1)
    ap.add_argument('--ps', default='1e4,1e5,1e6,1e7,1e8,1e9')
    ap.add_argument('--iters', type=int, default=2000)
    ap.add_argument('--rho', type=float, default=1e4)
    ap.add_argument('--sym', action='store_true')
    ap.add_argument('--cycles', type=int, default=1)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    h = load_heights(a.src)
    if a.up > 1:
        h = upsample(h, a.up)
    print('start n', len(h), 'score', f'{score(h):.13f}', flush=True)
    ps = [float(x) for x in a.ps.split(',')]
    best = (score(h), h)
    for cyc in range(a.cycles):
        h = polish(h, ps, a.iters, a.rho, a.sym, a.log, tag=f'c{cyc}')
        sc = score(h)
        if sc < best[0]:
            best = (sc, h.copy())
            save(h, a.out, note=f'cycle {cyc}')
        print(f'cycle {cyc} score {sc:.13f} best {best[0]:.13f}', flush=True)
