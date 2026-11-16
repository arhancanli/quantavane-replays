#!/usr/bin/env python3
"""5a [St26] C_5a <= 0.98173* and [CHO25] C_5a <= 0.98183, program A.

Written from Carter, Hunter and O'Bryant, "On the diameter of finite Sidon sets" (arXiv 2310.20032; Acta Math.
Hungar. 175 (2025)), Lemmas 3.1 and 3.2 and Proposition 1.3, read as text. With levels alpha_1 < ... < alpha_K,
tau and window ratios c_1, c_2, c_3, every Sidon set gives b <= max_w min(f_1(w), f_c1(w), f_c2(w), f_c3(w)) over
0 = w_0 <= w_1 <= ... <= w_K <= w_{K+1} = 1, and C_5a <= b/2. Here

  f_1(w) = tau + 1/tau - 2 tau sum_{j=0..K} (w_{j+1} - w_j) min_{alpha_j <= z <= alpha_{j+1}} (z - 1)^2,
  f_c(w) = c tau + 1/(c tau) - (2 tau / c^2) int_0^1 G(zeta(t), eta(t)) dt,

zeta(t) = #{j : w_j + c < t}, eta(t) = #{j : w_j < t}, G(z, e) = min (x - c)^2 over
alpha_{e-1} - alpha_z <= x <= alpha_e - alpha_{z-1} (alpha_{-1} = alpha_0 = 0, alpha_{K+1} = infinity).

Method of this program: the w-simplex is cut into the full-dimensional cells on which every difference
w_b - w_a lies strictly between two consecutive thresholds of {0, c's, 1}; cells are enumerated depth-first with an
incremental Floyd-Warshall test of the difference constraints in integers (strictness carried in the low bits). On
a cell all four bounds are affine; f_c's coefficients come from counting q's and r's in the merged order. A cell is
skipped when a box bound (from the shortest-path distances) of one of the four bounds is already <= a value attained
at a known point; on every other cell the linear program  max t  s.t.  t <= f_k(w), w in the closed cell  is solved
exactly by a two-phase tableau simplex over Fractions with Bland's rule. Standard library only.
"""
import sys, time
from fractions import Fraction as F

SCALE = 10 ** 6          # every parameter has six decimals
PARAMS = {
    # St26: records/certified/record_1771350836_46537e4f6f8d/certified_record.json of
    # rynald0cst0ltziam/Erdos-Sidon-Records-30 @ aa40ce9 (the six-decimal parameters of the certified record)
    "st26": dict(tau="1.128028", alphas=["0.708402", "0.785148", "0.869892", "1.124418", "1.182121", "1.243221"],
                 cs=["0.664209", "0.664134", "0.855601"], claim="0.98173",
                 point=["0.0943872", "0.335791", "0.3880977", "0.8328982", "1", "1"]),
    # CHO25: the proof of Theorem 3.3
    "cho25": dict(tau="1.12733", alphas=["0.70749", "0.78822", "0.87175", "1.12464", "1.18020", "1.24610"],
                  cs=["0.66461", "0.67780", "0.71884"], claim="0.98183",
                  point=["0.13398", "0.30015", "0.46220", "0.96476", "0.97795", "1"]),
}


def dec(s):
    return F(s)


class Problem:
    def __init__(self, tau, alphas, cs, lemma32_factor=2, drop_window=None):
        self.tau = tau
        self.K = len(alphas)
        self.al = alphas                       # alpha_1..alpha_K
        self.cs = [c for i, c in enumerate(cs) if not (drop_window and i in drop_window)]
        self.fac = lemma32_factor
        K = self.K
        # f_1 as an affine form [const, w_1..w_K] (w_0 = 0, w_{K+1} = 1)
        A = [F(0)] + list(alphas) + [None]     # alpha_0..alpha_{K+1}
        m = []
        for j in range(K + 1):
            lo, hi = A[j], A[j + 1]
            if hi is not None and hi < 1:
                m.append((1 - hi) ** 2)
            elif lo > 1:
                m.append((lo - 1) ** 2)
            else:
                m.append(F(0))
        form = [tau + 1 / tau] + [F(0)] * K
        for j in range(K + 1):           # - 2 tau m_j (w_{j+1} - w_j)
            k = -2 * tau * m[j]
            self._addw(form, j + 1, k)
            self._addw(form, j, -k)
        self.f1 = form

    def _addw(self, form, idx, coef):
        K = self.K
        if idx == 0:
            return
        if idx == K + 1:
            form[0] += coef
        else:
            form[idx] += coef

    def alpha(self, i):                  # alpha_i with alpha_{-1} = alpha_0 = 0, alpha_{K+1} = None (infinity)
        if i <= 0:
            return F(0)
        if i > self.K:
            return None
        return self.al[i - 1]

    def G(self, z, e, c):
        a = self.alpha(e - 1)
        b = self.alpha(z)
        lo = None if b is None else (a - b)            # a is finite since e - 1 <= K
        top = self.alpha(e)
        hi = None if top is None else top - self.alpha(z - 1)
        if hi is not None and hi < c:
            return (hi - c) ** 2
        if lo is not None and lo > c:
            return (lo - c) ** 2
        return F(0)

    def fc_form(self, c, less):
        """Affine form of f_c on a cell. less(a, b) is True when w_b < w_a + c on the cell (a < b)."""
        K = self.K
        n = K + 2
        # merge r_0 < ... < r_{K+1} (r_j = w_j) and q_0 < q_1 < ... (q_j = w_j + c) up to r_{K+1} = 1
        seq = []
        i = j = 0
        while i < n:
            # compare r_i with q_j: r_i < q_j iff (j >= i) or less(j, i) when j < i
            if j < n and not (j >= i or less(j, i)):
                seq.append(("q", j)); j += 1
            else:
                seq.append(("r", i)); i += 1
        form = [c * self.tau + 1 / (c * self.tau)] + [F(0)] * K
        k = -self.fac * self.tau / (c * c)
        nq = nr = 0
        for t in range(1, len(seq)):
            kind, idx = seq[t - 1]
            if kind == "q":
                nq += 1
            else:
                nr += 1
            g = self.G(nq, nr, c)
            if g == 0:
                continue
            # + k g (P_t - P_{t-1})
            for (kd, ix), sgn in ((seq[t], 1), (seq[t - 1], -1)):
                coef = k * g * sgn
                self._addw(form, ix, coef)
                if kd == "q":
                    form[0] += coef * c
        return form

    def fc_direct(self, c, w):
        """f_c at a point w = [w_0, ..., w_{K+1}] straight from the definition (for self-tests)."""
        pts = sorted(set([x for x in w] + [x + c for x in w if x + c <= 1]))
        tot = F(0)
        for a, b in zip(pts, pts[1:]):
            t = (a + b) / 2
            z = sum(1 for x in w if x + c < t)
            e = sum(1 for x in w if x < t)
            tot += (b - a) * self.G(z, e, c)
        return c * self.tau + 1 / (c * self.tau) - self.fac * self.tau / (c * c) * tot

    def f1_direct(self, w):
        A = [F(0)] + list(self.al) + [None]
        tot = F(0)
        for j in range(self.K + 1):
            lo, hi = A[j], A[j + 1]
            mm = (1 - hi) ** 2 if (hi is not None and hi < 1) else ((lo - 1) ** 2 if lo > 1 else F(0))
            tot += (w[j + 1] - w[j]) * mm
        return self.tau + 1 / self.tau - 2 * self.tau * tot


def ev(form, w):    # w = [w_1..w_K]
    return form[0] + sum(a * x for a, x in zip(form[1:], w))


# ---------------------------------------------------------------- exact LP: two-phase tableau simplex, Bland's rule
def simplex_max(nvar, rows, obj, want_duals=False):
    """maximize obj.x subject to rows (coeffs, sense, rhs), x >= 0. Returns (value, x) or (None, None) if infeasible."""
    # normalise to rhs >= 0
    R = []
    for a, s, b in rows:
        a = list(a)
        if b < 0:
            a = [-x for x in a]; b = -b; s = {"<=": ">=", ">=": "<=", "==": "=="}[s]
        R.append((a, s, b))
    m = len(R)
    nsl = sum(1 for _, s, _ in R if s != "==")
    nart = sum(1 for _, s, b in R if s == ">=" or s == "==")
    ncol = nvar + nsl + nart
    T = []
    basis = []
    si = nvar
    ai = nvar + nsl
    arts = []
    for a, s, b in R:
        row = a + [F(0)] * (nsl + nart) + [b]
        if s == "<=":
            row[si] = F(1); basis.append(si); si += 1
        elif s == ">=":
            row[si] = F(-1); si += 1
            row[ai] = F(1); basis.append(ai); arts.append(ai); ai += 1
        else:
            row[ai] = F(1); basis.append(ai); arts.append(ai); ai += 1
        T.append(row)

    def pivot(r, col):
        pr = T[r]
        p = pr[col]
        if p != 1:
            pr = [x / p for x in pr]; T[r] = pr
        for i in range(m):
            if i != r:
                f = T[i][col]
                if f:
                    Ti = T[i]
                    T[i] = [x - f * y for x, y in zip(Ti, pr)]
        basis[r] = col

    def run(cost, allowed):
        # maximise cost.x ; reduced costs computed each iteration (small tableaux)
        while True:
            cb = [cost[basis[i]] for i in range(m)]
            enter = None
            for j in range(ncol):
                if not allowed[j] or j in basis:
                    continue
                rc = cost[j] - sum(cb[i] * T[i][j] for i in range(m) if T[i][j])
                if rc > 0:
                    enter = j; break           # Bland: least index
            if enter is None:
                return True
            best = None
            for i in range(m):
                a = T[i][enter]
                if a > 0:
                    ratio = T[i][-1] / a
                    if best is None or ratio < best[0] or (ratio == best[0] and basis[i] < basis[best[1]]):
                        best = (ratio, i)
            if best is None:
                return False                   # unbounded
            pivot(best[1], enter)

    allowed = [True] * ncol
    if arts:
        cost1 = [F(0)] * ncol
        for j in arts:
            cost1[j] = F(-1)
        run(cost1, allowed)
        if sum(T[i][-1] for i in range(m) if basis[i] in arts) != 0:
            return (None, None, None) if want_duals else (None, None)
        # drive remaining (zero-level) artificials out of the basis where possible
        for i in range(m):
            if basis[i] in arts:
                for j in range(nvar + nsl):
                    if T[i][j] != 0 and j not in basis:
                        pivot(i, j); break
        for j in arts:
            allowed[j] = False
    cost2 = list(obj) + [F(0)] * (ncol - nvar)
    if not run(cost2, allowed):
        raise RuntimeError("unbounded LP")
    x = [F(0)] * nvar
    for i in range(m):
        if basis[i] < nvar:
            x[basis[i]] = T[i][-1]
    if want_duals:
        duals = {}
        col = nvar
        for r, (a0, s0, b0) in enumerate(rows):
            if s0 == "==":
                continue
            if s0 == "<=" and b0 >= 0:
                duals[r] = sum(cost2[basis[i]] * T[i][col] for i in range(m))
            col += 1
        return sum(o * v for o, v in zip(obj, x)), x, duals
    return sum(o * v for o, v in zip(obj, x)), x


# ---------------------------------------------------------------- cells
def enumerate_cells(K, thr_int):
    """Depth-first over levels of every pair difference d_ab = w_b - w_a (0 <= a < b <= K+1) among the thresholds
    0 = t_0 < t_1 < ... < t_m < 1; yields (levels, dist) for every full-dimensional cell. Integers: value*64 - strict."""
    n = K + 2
    pairs = sorted(((a, b) for a in range(n) for b in range(a + 1, n)), key=lambda p: (p[1] - p[0], p[0]))
    INF = 1 << 62
    D0 = [[0 if i == j else INF for j in range(n)] for i in range(n)]

    def add(D, u, v, W):     # constraint w_v - w_u <= W (encoded); returns new matrix or None
        if D[v][u] < INF and D[v][u] + W < 0:
            return None
        if D[u][v] <= W:
            return D
        N = [row[:] for row in D]
        for i in range(n):
            if D[i][u] >= INF:
                continue
            diu = D[i][u] + W
            Ni = N[i]
            for j in range(n):
                if D[v][j] < INF:
                    x = diu + D[v][j]
                    if x < Ni[j]:
                        Ni[j] = x
        return N

    one = SCALE * 64
    D0 = add(D0, 0, n - 1, one)
    D0 = add(D0, n - 1, 0, -one)
    T = [0] + thr_int + [None]      # level L: T[L] < d < T[L+1]  (T[m+1] = None: no upper bound below 1)
    nlev = len(thr_int) + 1
    levels = {}
    out = []

    def rec(k, D):
        if k == len(pairs):
            out.append((dict(levels), D))
            return
        a, b = pairs[k]
        lo = 0
        if b - a > 1:
            lo = max(levels[(a + 1, b)], levels[(a, b - 1)])
        for L in range(lo, nlev):
            D1 = add(D, b, a, -T[L] * 64 - 1)              # d_ab > T[L]
            if D1 is None:
                break                                       # higher levels are infeasible too
            if T[L + 1] is not None:
                D1 = add(D1, a, b, T[L + 1] * 64 - 1)     # d_ab < T[L+1]
                if D1 is None:
                    continue
            levels[(a, b)] = L
            rec(k + 1, D1)
        levels.pop((a, b), None)

    rec(0, D0)
    return out


def main():
    args = sys.argv[1:]
    which = args[0] if args else "st26"
    P = PARAMS[which]
    tau = dec(P["tau"]); alphas = [dec(x) for x in P["alphas"]]; cs = [dec(x) for x in P["cs"]]
    claim = dec(P["claim"])
    fac = 2; drop = None; label = "main"
    if "--claim" in args:
        claim = dec(args[args.index("--claim") + 1]); label = "CONTROL claim %s" % claim
    if "--factor1" in args:
        fac = 1; label = "CONTROL Lemma 3.2 with tau/c^2 in place of 2 tau/c^2"
    if "--drop" in args:
        drop = [int(x) for x in args[args.index("--drop") + 1].split(",")]
        label = "CONTROL windows %s removed" % ", ".join("c_%d" % (i + 1) for i in drop)
    if "--alpha4" in args:
        alphas[3] = dec(args[args.index("--alpha4") + 1]); label = "CONTROL alpha_4 = %s" % alphas[3]
    t0 = time.time()
    pr = Problem(tau, alphas, cs, fac, drop)
    K = pr.K
    print("program A, parameters %s (%s)" % (which, label))
    print("tau = %s, alphas = %s, cs = %s" % (P["tau"], ", ".join(P["alphas"]), ", ".join(str(float(c)) for c in pr.cs)))

    # self-tests: the cell formula against the definition at random points, and the paper's own worst case
    import random
    rnd = random.Random(12345)
    pt = [F(0)] + [dec(x) for x in P["point"]] + [F(1)]
    for c in pr.cs:
        for _ in range(200):
            w = [F(0)] + sorted(F(rnd.randrange(0, SCALE + 1), SCALE) for _ in range(K)) + [F(1)]
            less = lambda a, b, w=w, c=c: w[b] < w[a] + c
            if ev(pr.fc_form(c, less), w[1:-1]) != pr.fc_direct(c, w):
                print("SELF-TEST FAILED: cell formula of f_c at", w); sys.exit(2)
    for _ in range(50):
        w = [F(0)] + sorted(F(rnd.randrange(0, SCALE + 1), SCALE) for _ in range(K)) + [F(1)]
        if ev(pr.f1, w[1:-1]) != pr.f1_direct(w):
            print("SELF-TEST FAILED: f_1"); sys.exit(2)
    vals = [pr.f1_direct(pt)] + [pr.fc_direct(c, pt) for c in pr.cs]
    LB = min(vals)
    print("self-test: cell formulas equal the definition at 200 random points per window: ok")
    print("the four bounds at w = (%s): %s" % (", ".join(P["point"]), ", ".join("%.9f" % float(v) for v in vals)))

    thr = sorted(set(int(c * SCALE) for c in pr.cs))
    assert all(F(t, SCALE) in pr.cs for t in thr)
    cells = enumerate_cells(K, thr)
    print("full-dimensional cells: %d (%.1f s)" % (len(cells), time.time() - t0), flush=True)
    thrF = [F(0)] + [F(t, SCALE) for t in thr] + [None]
    best = (LB, None, None)
    nlp = 0

    def boxmax(f, lo, hi):
        return f[0] + sum(max(a * lo[j + 1], a * hi[j + 1]) for j, a in enumerate(f[1:]))

    scored = []
    for levels, D in cells:
        def less_for(c, levels=levels):
            return lambda a, b: (thrF[levels[(a, b)] + 1] is not None and thrF[levels[(a, b)] + 1] <= c)
        forms = [pr.f1] + [pr.fc_form(c, less_for(c)) for c in pr.cs]
        # box of the cell from the shortest-path distances: lo_j <= w_j <= hi_j
        hi = [F(0)] + [F((D[0][j] + 63) // 64, SCALE) for j in range(1, K + 1)]
        lo = [F(0)] + [F(-((D[j][0] + 63) // 64), SCALE) for j in range(1, K + 1)]
        scored.append((min(boxmax(f, lo, hi) for f in forms), levels, forms, lo, hi))
    # most promising cells first; a cell is skipped when a valid upper bound of its max-min is <= the best value so
    # far: min_k (box max of f_k), or the box max of sum_k lam_k f_k for the dual weights lam of the best cell so far
    scored.sort(key=lambda r: r[0], reverse=True)
    lam = None
    nskip = 0
    for ncell, (ub, levels, forms, lo, hi) in enumerate(scored):
        if ncell % 2000 == 0 and "--progress" in args:
            print("  cell %d, LPs %d, best %.12f, %.0f s" % (ncell, nlp, float(best[0]), time.time() - t0), flush=True)
        if ub <= best[0]:
            nskip += len(scored) - ncell
            break
        if lam is not None:
            comb = [sum(l * f[i] for l, f in zip(lam, forms)) for i in range(K + 1)]
            if boxmax(comb, lo, hi) <= best[0]:
                nskip += 1
                continue
        # LP: variables w_1..w_K, t (all >= 0); maximise t
        rows = []
        nv = K + 1
        def vec(d):
            v = [F(0)] * nv
            for i, x in d.items():
                v[i] = x
            return v
        n = K + 2
        for (a, b), L in levels.items():
            for bound, sense in ((thrF[L], ">="), (thrF[L + 1], "<=")):
                if bound is None or (bound == 0 and b - a > 1):
                    continue
                # implied by a pair inside (lower) or around (upper) with the same level: d is monotone in the pair
                if sense == ">=" and b - a > 1 and (levels[(a + 1, b)] == L or levels[(a, b - 1)] == L):
                    continue
                if sense == "<=" and ((a > 0 and levels[(a - 1, b)] == L) or (b < n - 1 and levels[(a, b + 1)] == L)):
                    continue
                d = {}; rhs = bound
                if b != n - 1: d[b - 1] = d.get(b - 1, 0) + 1
                else: rhs -= 1
                if a != 0: d[a - 1] = d.get(a - 1, 0) - 1
                if b == n - 1 and a == 0:
                    continue
                rows.append((vec(d), sense, rhs))
        frow = len(rows)
        for f in forms:      # t - f_lin.w <= f_const
            d = {i: -x for i, x in enumerate(f[1:]) if x}
            d[K] = F(1)
            rows.append((vec(d), "<=", f[0]))
        val, x, duals = simplex_max(nv, rows, [F(0)] * K + [F(1)], want_duals=True)
        nlp += 1
        if val is not None and val > best[0]:
            best = (val, x, [ev(f, x[:K]) for f in forms])
            if all(frow + k in duals for k in range(len(forms))):
                cand = [duals[frow + k] for k in range(len(forms))]
                if sum(cand) == 1 and min(cand) >= 0:      # then min_k f_k <= sum_k cand_k f_k everywhere
                    lam = cand
    print("cells skipped by the upper bounds: %d" % nskip)
    if lam is not None:
        print("dual weights of the four bounds at the maximum: " + ", ".join(decstr(v, 12) for v in lam))
    b = best[0]
    print("linear programs solved exactly: %d" % nlp)
    print("max over w of min of the four bounds (exact):")
    print("b* = %s" % b)
    print("b* = " + decstr(b, 20))
    if best[1] is not None:
        print("attained at w = (%s)" % ", ".join(decstr(v, 12) for v in best[1][:K]))
        print("the four bounds there: %s" % ", ".join(decstr(v, 15) for v in best[2]))
    else:
        print("attained at the given point (no cell beats it)")
    Cst = b / 2
    print("C_5a <= b*/2 = " + decstr(Cst, 20))
    print("row %s: C_5a <= %s  %s" % (which, decstr(claim, 6), "holds" if Cst <= claim else "NOT implied"))
    print("VERDICT %s" % (("PASS" if Cst <= claim else "FAILED") if label == "main" else
                          ("REFUSED" if Cst > claim else "NOT REFUSED")))
    print("run time %.1f s" % (time.time() - t0))


def decstr(x, digits):
    x = F(x)
    neg = x < 0
    x = abs(x)
    q = x.numerator * 10 ** digits // x.denominator
    s = str(q).rjust(digits + 1, "0")
    return ("-" if neg else "") + s[:-digits] + "." + s[-digits:]


if __name__ == "__main__":
    main()
