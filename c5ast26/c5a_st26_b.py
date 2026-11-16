#!/usr/bin/env python3
"""5a [St26] C_5a <= 0.98173* and [CHO25] C_5a <= 0.98183, program B (written separately from program A).

Source read as text: Carter, Hunter and O'Bryant, arXiv 2310.20032 (Acta Math. Hungar. 175 (2025)), Section 3:
Lemma 3.1 (one bound, affine in w), Lemma 3.2 (one bound per window ratio c, piecewise affine), Section 3.2.1
(cells from the interlacings s of q_j = w_j + c and r_j = w_j; p, zeta, eta) and Proposition 1.3 (C_5a <= b/2).

Route of this program, different from program A at every step:
- cells: for each window c, the interlacing tuple s (s_j = the largest index i with r_i <= q_j), chosen jointly for
  the three windows j by j; strict feasibility of the difference constraints by Bellman-Ford from a virtual source
  (integers at scale 16 * 10^6, a strict bound tightened by one unit; a cycle has at most 8 edges);
- cell functions: the paper's p_0 < ... < p_L, zeta_j = min{j' : q_j' >= p_j}, eta_j = min{j' : r_j' >= p_j} and
  the sum over j of (p_j - p_{j-1}) min (z - c)^2 over [alpha_{eta_j - 1} - alpha_{zeta_j}, alpha_{eta_j} - alpha_{zeta_j - 1}];
- pruning: an upper bound of each bound over the cell from the ranges of the increments w_{j+1} - w_j (greedy over
  sum = 1);
- the linear program in its DUAL form, min F.lambda + h.mu subject to sum lambda = 1, A^T lambda = G^T mu,
  lambda, mu >= 0, by a revised-free two-phase simplex on the equality system (7 rows), exact Fractions, Bland's rule.
Standard library only.
"""
import sys, time
from fractions import Fraction as Q

PAR = {
    "st26": ("1.128028", "0.708402 0.785148 0.869892 1.124418 1.182121 1.243221", "0.664209 0.664134 0.855601",
             "0.98173", "0.0943872 0.335791 0.3880977 0.8328982 1 1"),
    "cho25": ("1.12733", "0.70749 0.78822 0.87175 1.12464 1.18020 1.24610", "0.66461 0.67780 0.71884",
              "0.98183", "0.13398 0.30015 0.46220 0.96476 0.97795 1"),
}
U = 10 ** 6
S = 16 * U


def sqmin(lo, hi, c):
    """min (z - c)^2 over lo <= z <= hi; None stands for -infinity (lo) or +infinity (hi)."""
    if hi is not None and hi < c:
        return (c - hi) * (c - hi)
    if lo is not None and lo > c:
        return (lo - c) * (lo - c)
    return Q(0)


class Bounds:
    def __init__(self, tau, al, cs, factor=2):
        self.tau, self.al, self.cs, self.factor = tau, al, cs, factor
        self.K = len(al)

    def A(self, i):
        if i < 1:
            return Q(0)
        return None if i > self.K else self.al[i - 1]

    def lemma31(self):
        """coefficients over the increments delta_j = w_{j+1} - w_j, j = 0..K, and the constant"""
        t = self.tau
        e = []
        for j in range(self.K + 1):
            e.append(-2 * t * sqmin(self.A(j), self.A(j + 1), Q(1)))
        return t + 1 / t, e

    def lemma32(self, c, order):
        """order: the symbolic increasing list of the points in P, each ('r', i) or ('q', i), from r_0 = 0 to
        r_{K+1} = 1. Returns the bound as constant + coefficients over w_1..w_K."""
        K = self.K
        t = self.tau
        L = len(order) - 1
        const = c * t + 1 / (c * t)
        coef = [Q(0)] * (K + 2)
        k = -self.factor * t / (c * c)
        pos = {pt: n for n, pt in enumerate(order)}
        for j in range(1, L + 1):
            # zeta_j: smallest j' with q_j' >= p_j ; eta_j: smallest j' with r_j' >= p_j (ties: the order itself)
            z = next(i for i in range(K + 2) if ('q', i) not in pos or pos[('q', i)] >= j)
            e = next(i for i in range(K + 2) if pos[('r', i)] >= j)
            Aem1 = self.A(e - 1)
            Az = self.A(z)
            lo = None if Az is None else Aem1 - Az
            Ae = self.A(e)
            hi = None if Ae is None else Ae - self.A(z - 1)
            m = sqmin(lo, hi, c)
            if m == 0:
                continue
            for (kind, i), sg in ((order[j], 1), (order[j - 1], -1)):
                coef[i] += k * m * sg
                if kind == 'q':
                    const += k * m * sg * c
        const += coef[K + 1]            # w_{K+1} = 1
        return const, coef[1:K + 1]

    def direct(self, c, w):
        """Lemma 3.2's bound at a point w = (w_0, ..., w_{K+1}), from x_i, y_i of its proof (fine grid in t)."""
        K = self.K
        brk = sorted(set(list(w) + [x + c for x in w]))
        brk = [x for x in brk if x <= 1]
        tot = Q(0)
        for a, b in zip(brk, brk[1:]):
            tm = (a + b) / 2
            z = min(i for i in range(K + 2) if w[i] + c >= tm)
            e = min(i for i in range(K + 2) if w[i] >= tm)
            Az = self.A(z)
            lo = None if Az is None else self.A(e - 1) - Az
            hi = None if self.A(e) is None else self.A(e) - self.A(z - 1)
            tot += (b - a) * sqmin(lo, hi, c)
        return c * self.tau + 1 / (c * self.tau) - self.factor * self.tau / (c * c) * tot


def bellman_ford(n, edges, start=None):
    """edges (u, v, W): x_v - x_u <= W. Returns a feasible potential, or None if there is a negative cycle.
    Starting from any potential is valid (it acts as a virtual source with those edge weights); n + 1 passes."""
    d = list(start) if start is not None else [0] * n
    for _ in range(n + 1):
        ch = False
        for u, v, W in edges:
            if d[u] + W < d[v]:
                d[v] = d[u] + W; ch = True
        if not ch:
            return d
    return None


def dist_from(n, edges, src):
    INF = None
    d = [INF] * n
    d[src] = 0
    for _ in range(n):
        for u, v, W in edges:
            if d[u] is not None and (d[v] is None or d[u] + W < d[v]):
                d[v] = d[u] + W
    return d


def cells(K, cint):
    """Joint interlacings for the windows (c's as integers at scale U, any order). Yields (s tuples, edges)."""
    n = K + 2
    nw = len(cint)
    base = [(n - 1, 0, -S), (0, n - 1, S)]                 # w_{K+1} - w_0 = 1
    base += [(j + 1, j, -1) for j in range(n - 1)]          # w_j < w_{j+1}  (x_j - x_{j+1} <= -1 unit)
    out = []

    def cons(j, s, c):
        """q_j = w_j + c lies after r_s and before r_{s+1} (strictly)."""
        e = []
        # w_s - w_j < c   (if s > j; for s == j it is automatic)
        if s > j:
            e.append((j, s, 16 * c - 1))
        if s + 1 <= n - 1:
            # w_{s+1} - w_j > c  ->  w_j - w_{s+1} < -c
            e.append((s + 1, j, -16 * c - 1))
        return e

    def rec(j, prev, edges, acc, pot):
        if j == n:
            out.append((acc[:], edges))
            return
        def rec_w(wi, cur, edges2, pot2):
            if wi == nw:
                rec(j + 1, cur, edges2, acc + [tuple(cur)], pot2)
                return
            lo = max(j, prev[wi] if prev else j)
            for s in range(lo, n):
                e2 = edges2 + cons(j, s, cint[wi])
                p2 = bellman_ford(n, e2, pot2)
                if p2 is not None:
                    rec_w(wi + 1, cur + [s], e2, p2)
        rec_w(0, [], edges, pot)

    rec(0, None, base, [], bellman_ford(n, base))
    return out


def order_for(K, sj):
    """symbolic order of P for one window: r's, and q_j right after r_{s_j} (q_j <= 1 iff s_j < K+1)."""
    n = K + 2
    o = []
    for i in range(n):
        o.append(('r', i))
        if i == n - 1:
            break
        o.extend(('q', j) for j in range(n) if sj[j] == i)
    return o


def dual_lp(forms, G, h):
    """min over lambda (one per form), mu (one per row of G w <= h) of sum lambda_k F_k + sum mu_i h_i
    subject to sum lambda = 1 and sum_k lambda_k a_k - sum_i mu_i g_i = 0 (one row per w-coordinate)."""
    nk = len(forms)
    K = len(forms[0][1])
    cols = []      # each column: (cost, [entries of the 1 + K equality rows])
    for F0, a in forms:
        cols.append((F0, [Q(1)] + list(a)))
    for g, hh in zip(G, h):
        cols.append((hh, [Q(0)] + [-x for x in g]))
    m = 1 + K
    rhs = [Q(1)] + [Q(0)] * K
    ncol = len(cols)
    # tableau with artificials (one per row); make rhs >= 0 (already)
    T = []
    for r in range(m):
        T.append([cols[c][1][r] for c in range(ncol)] + [Q(1) if rr == r else Q(0) for rr in range(m)] + [rhs[r]])
    basis = [ncol + r for r in range(m)]
    tot = ncol + m

    def piv(r, c):
        p = T[r][c]
        T[r] = [x / p for x in T[r]]
        for i in range(m):
            if i != r and T[i][c] != 0:
                f = T[i][c]
                T[i] = [x - f * y for x, y in zip(T[i], T[r])]
        basis[r] = c

    def solve(cost, banned):
        # minimise cost.x
        while True:
            enter = None
            for c in range(tot):
                if c in banned or c in basis:
                    continue
                red = cost[c] - sum(cost[basis[i]] * T[i][c] for i in range(m))
                if red < 0:
                    enter = c; break
            if enter is None:
                return
            cand = [(T[i][-1] / T[i][enter], basis[i], i) for i in range(m) if T[i][enter] > 0]
            if not cand:
                raise RuntimeError("dual unbounded: the cell would be empty")
            piv(min(cand)[2], enter)

    c1 = [Q(0)] * ncol + [Q(1)] * m
    solve(c1, set())
    if sum(T[i][-1] for i in range(m) if basis[i] >= ncol) != 0:
        raise RuntimeError("dual infeasible")
    for i in range(m):
        if basis[i] >= ncol:
            for c in range(ncol):
                if T[i][c] != 0 and c not in basis:
                    piv(i, c); break
    c2 = [cols[c][0] for c in range(ncol)] + [Q(0)] * m
    solve(c2, set(range(ncol, tot)))
    val = sum(c2[basis[i]] * T[i][-1] for i in range(m))
    lam = [Q(0)] * nk
    for i in range(m):
        if basis[i] < nk:
            lam[basis[i]] = T[i][-1]
    return val, lam


def main():
    a = sys.argv[1:]
    which = a[0] if a else "st26"
    tau, al, cs, claim, pt = PAR[which]
    tau = Q(tau); al = [Q(x) for x in al.split()]; cs = [Q(x) for x in cs.split()]; claim = Q(claim)
    pt = [Q(0)] + [Q(x) for x in pt.split()] + [Q(1)]
    factor = 2; tag = "main"
    if "--claim" in a:
        claim = Q(a[a.index("--claim") + 1]); tag = "CONTROL claim %s" % claim
    if "--factor1" in a:
        factor = 1; tag = "CONTROL Lemma 3.2 with tau/c^2 in place of 2 tau/c^2"
    if "--drop" in a:
        gone = set(int(x) for x in a[a.index("--drop") + 1].split(","))
        cs = [c for i, c in enumerate(cs) if i not in gone]
        tag = "CONTROL windows %s removed" % ", ".join("c_%d" % (i + 1) for i in sorted(gone))
    if "--alpha4" in a:
        al[3] = Q(a[a.index("--alpha4") + 1]); tag = "CONTROL alpha_4 = %s" % al[3]
    t0 = time.time()
    B = Bounds(tau, al, cs, factor)
    K = B.K
    print("program B, parameters %s (%s)" % (which, tag))
    F1c, F1e = B.lemma31()
    # Lemma 3.1 over w_1..w_K: sum e_j (w_{j+1} - w_j)
    f1_const = F1c + F1e[K]                       # e_K * w_{K+1} with w_{K+1} = 1
    f1_lin = [F1e[j - 1] - F1e[j] for j in range(1, K + 1)]

    def f1_at(w):
        return F1c + sum(F1e[j] * (w[j + 1] - w[j]) for j in range(K + 1))

    vals = [f1_at(pt)] + [B.direct(c, pt) for c in cs]
    LB = min(vals)
    print("bounds at the given point: " + ", ".join("%.9f" % float(v) for v in vals))
    cint = [int(c * U) for c in cs]
    assert all(Q(ci, U) == c for ci, c in zip(cint, cs))
    allc = cells(K, cint)
    print("full-dimensional cells (joint interlacings): %d" % len(allc))
    # self-test: the cell function against the direct evaluation at the cell's own interior point
    import random
    rr = random.Random(7)
    n = K + 2
    best = (LB, None)
    nlp = 0
    tested = 0

    def incr_max(form, up):
        """max of the form over {delta_j in [0, up_j], sum delta = 1} (a relaxation of the closed cell), greedy."""
        c0, lin = form
        ed = [sum(lin[i - 1] for i in range(j + 1, K + 1)) for j in range(n - 1)]
        rem = Q(1); val = c0
        for e, j in sorted(((ed[j], j) for j in range(n - 1)), reverse=True):
            take = min(rem, up[j])
            val += e * take; rem -= take
            if rem == 0:
                break
        return val

    prepared = []
    for sjs, edges in allc:
        forms = [(f1_const, f1_lin)]
        for wi, c in enumerate(cs):
            sj = [s[wi] for s in sjs]
            forms.append(B.lemma32(c, order_for(K, sj)))
        # interior point from shortest paths (a feasible point of the strict system), used as self-test
        if tested < 300 and rr.random() < 0.05:
            d = dist_from(n + 1, edges + [(n, i, 0) for i in range(n)], n)
            w = [Q(d[i] - d[0], S) for i in range(n)]
            if w[n - 1] == 1:
                for wi, c in enumerate(cs):
                    lhs = forms[wi + 1][0] + sum(x * y for x, y in zip(forms[wi + 1][1], w[1:n - 1]))
                    if lhs != B.direct(c, w):
                        print("SELF-TEST FAILED at cell", sjs); sys.exit(2)
                tested += 1
        # the closed cell: strict bounds back to their exact values (16 X - 1 -> 16 X)
        closed = [(u, v, -(-W // 16) * 16) for u, v, W in edges]
        # ranges of the increments delta_j = w_{j+1} - w_j over the closed cell (shortest paths; upper ends)
        up = []
        for j in range(n - 1):
            d = dist_from(n, closed, j)
            up.append(Q(d[j + 1], S))
        prepared.append((min(incr_max(f, up) for f in forms), sjs, forms, closed, up))
    # largest relaxation bound first; skip a cell when its bound, or the relaxation bound of sum lam_k f_k for the
    # dual weights lam of the best cell so far (min_k f_k <= sum lam_k f_k when lam >= 0, sum lam = 1), is <= best
    prepared.sort(key=lambda r: r[0], reverse=True)
    lam_best = None
    skipped = 0
    for idx, (u0, sjs, forms, closed, up) in enumerate(prepared):
        if u0 <= best[0]:
            skipped += len(prepared) - idx
            break
        if lam_best is not None:
            comb = (sum(l * f[0] for l, f in zip(lam_best, forms)),
                    [sum(l * f[1][i] for l, f in zip(lam_best, forms)) for i in range(K)])
            if incr_max(comb, up) <= best[0]:
                skipped += 1
                continue
        # constraints G w <= h of the closed cell, from the edges: x_v - x_u <= W/S
        G = []; h = []
        for u, v, W in closed:
            if (u, v) in ((n - 1, 0), (0, n - 1)):
                continue
            Wv = Q(W, S)
            g = [Q(0)] * K; rhs = Wv
            if 1 <= v <= K: g[v - 1] += 1
            elif v == n - 1: rhs -= 1
            if 1 <= u <= K: g[u - 1] -= 1
            elif u == n - 1: rhs += 1
            G.append(g); h.append(rhs)
        val, lam = dual_lp(forms, G, h)
        nlp += 1
        if val > best[0]:
            best = (val, (sjs, lam))
            if sum(lam) == 1 and min(lam) >= 0:
                lam_best = lam
    print("cells skipped by the relaxation bounds: %d" % skipped)
    b = best[0]
    print("self-tests of cell functions at interior points: %d cells, all equal" % tested)
    print("dual linear programs solved exactly: %d" % nlp)
    print("b* = %s" % b)
    print("b* = " + dstr(b, 20))
    if best[1]:
        print("dual weights on (Lemma 3.1, windows): " + ", ".join(dstr(x, 12) for x in best[1][1]))
    C = b / 2
    print("C_5a <= b*/2 = " + dstr(C, 20))
    ok = C <= claim
    print("row %s: C_5a <= %s  %s" % (which, dstr(claim, 6), "holds" if ok else "NOT implied"))
    print("VERDICT " + (("PASS" if ok else "FAILED") if tag == "main" else ("REFUSED" if not ok else "NOT REFUSED")))
    print("run time %.1f s" % (time.time() - t0))


def dstr(x, k):
    x = Q(x); s = "-" if x < 0 else ""; x = abs(x)
    q = x.numerator * 10 ** k // x.denominator
    t = str(q).rjust(k + 1, "0")
    return s + t[:-k] + "." + t[-k:]


if __name__ == "__main__":
    main()
