"""Search (floating point, information only): Delsarte's LP bound for the kissing number in dimension 5.

minimize f(1) = 1 + sum_k f_k over f = 1 + sum_{k=1..d} f_k G_k, f_k >= 0, f(t) <= 0 on a grid of [-1, 1/2],
G_k the Gegenbauer polynomials of S^4 normalised by G_k(1) = 1: (k + 3) G_{k+1} = (2k + 3) t G_k - k G_{k-1}.
Solved through its dual  max sum_i y_i  s.t.  sum_i y_i (-G_k(t_i)) <= 1 (k = 1..d), y >= 0  by a dense simplex;
the primal f_k are the dual prices. Usage: python3 -I lp_search.py D M
"""
import sys


def gegen(d, t):
    g = [1.0, t]
    for k in range(1, d):
        g.append(((2 * k + 3) * t * g[k] - k * g[k - 1]) / (k + 3))
    return g[: d + 1]


def simplex_max(A, b, c, iters=100000):
    """max c.y s.t. A y <= b, y >= 0, b >= 0 (dense tableau, Bland's rule on ties). Returns (value, y, prices)."""
    m, n = len(A), len(A[0])
    T = [A[i][:] + [1.0 if j == i else 0.0 for j in range(m)] + [b[i]] for i in range(m)]
    z = [-x for x in c] + [0.0] * m + [0.0]
    basis = [n + i for i in range(m)]
    eps = 1e-12
    for _ in range(iters):
        # entering: most negative reduced cost
        e = min(range(n + m), key=lambda j: z[j])
        if z[e] >= -eps:
            break
        best, r = None, None
        for i in range(m):
            a = T[i][e]
            if a > eps:
                q = T[i][-1] / a
                if best is None or q < best - 1e-15:
                    best, r = q, i
        if r is None:
            raise RuntimeError("unbounded")
        piv = T[r][e]
        T[r] = [x / piv for x in T[r]]
        for i in range(m):
            if i != r and T[i][e] != 0.0:
                f = T[i][e]
                Ti, Tr = T[i], T[r]
                T[i] = [Ti[j] - f * Tr[j] for j in range(len(Ti))]
        f = z[e]
        z = [z[j] - f * T[r][j] for j in range(len(z))]
        basis[r] = e
    y = [0.0] * n
    for i, bi in enumerate(basis):
        if bi < n:
            y[bi] = T[i][-1]
    prices = z[n:n + m]
    return z[-1], y, prices


def main():
    d, M = int(sys.argv[1]), int(sys.argv[2])
    ts = [-1 + 1.5 * i / (M - 1) for i in range(M)]
    G = [gegen(d, t) for t in ts]
    A = [[-G[i][k] for i in range(M)] for k in range(1, d + 1)]
    b = [1.0] * d
    c = [1.0] * M
    val, y, prices = simplex_max(A, b, c)
    f = [1.0] + prices
    print("degree", d, "grid", M, "LP value (bound f(1)/f_0) =", 1 + val, "check", sum(f))
    print("f_k:", " ".join(repr(x) for x in f))
    mx = max(sum(f[k] * gk for k, gk in enumerate(gegen(d, -1 + 1.5 * i / 20000))) for i in range(20001))
    print("max of f on a fine grid of [-1, 1/2]:", mx)
    sup = [(ts[i], y[i]) for i in range(M) if y[i] > 1e-12]
    print("dual support:", " ".join("%.6f:%r" % s for s in sup))


if __name__ == "__main__":
    main()
