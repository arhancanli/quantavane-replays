"""Our own construction for C_1c lower bounds in the spike-plus-profile family described in teorth/optimizationproblems
PR #354 (SeanSagait): a geometric spike s_i = M r^i on cells 0..m-1, then f = 0, then a profile g on [j1, K] with every
lag k in [j1, K] active: sum_{i<m} s_i g_{k+i} + (g*g)_k = T.  Written from the PR's description only (its programs
were not read).  The equations are solved by a damped fixed point: H_k = (T - c_k)/M with c = g autocorrelation,
g_k = H_k - r H_{k+1} + r^m g_{k+m} (backward), negatives clipped."""
import numpy as np
from scipy.fft import rfft, irfft, next_fast_len

def autocorr(x, upto):
    L = next_fast_len(2 * len(x), real=True)
    X = rfft(x, L)
    return irfft(X * np.conj(X), L)[:upto]

def recursion(b, q, m):
    # g_k = b_k + q g_{k+m}, g_k = 0 for k >= len(b)
    if m == 1:
        out = np.empty_like(b); acc = 0.0
        # g_k = b_k + q g_{k+1}: a first-order filter run backwards
        from scipy.signal import lfilter
        return lfilter([1.0], [1.0, -q], b[::-1])[::-1]
    n = len(b); nb = -(-n // m); pad = nb * m - n
    B = np.concatenate([b, np.zeros(pad)]).reshape(nb, m)
    G = np.empty_like(B); G[-1] = B[-1]
    for j in range(nb - 2, -1, -1):
        G[j] = B[j] + q * G[j + 1]
    return G.reshape(-1)[:n]

def solve(K, m, r, M, T=1.0, iters=400, alpha=0.5, g=None, tol=1e-15, verbose=False):
    n = K + 1
    if g is None:
        g = np.zeros(n); g[n // 6:] = T * (1 - r) / M
    q = r ** m
    prev = None; a = alpha
    for it in range(iters):
        c = autocorr(g, n)
        H = (T - c) / M
        b = H - r * np.concatenate([H[1:], [0.0]])
        gn = recursion(b, q, m)
        gn[:m + 1] = 0.0
        np.maximum(gn, 0.0, out=gn)
        d = np.abs(gn - g).max() / max(gn.max(), 1e-300)
        if prev is not None and d > prev:
            a = max(a * 0.7, 0.01)
        else:
            a = min(a * 1.03, alpha)
        prev = d
        g = (1 - a) * g + a * gn
        if verbose and it % 50 == 0:
            print(it, d, flush=True)
        if d < tol:
            break
    return g, it, d

def weights(K, m, r, M, g):
    w = g.copy(); w[:m] += M * r ** np.arange(m)
    return w

def value(K, w):
    A = autocorr(w, K + 1)
    S = w.sum()
    rr = K * A / S ** 2
    return rr.min(), int(rr.argmin()), rr
