"""37a: independent replay of the statement of PR #283 (DennisKoshta), 'no Boolean function of degree 4 has
sensitivity 10', from the definition on the 37a page and the certificate read as data (Quantavane replays).
Standard library only; run with python3 -I.  The packet's programs (checker/, generator/, verifier/, sharpness/,
tools/) were neither read nor run; only the three certificate files are read, in the format stated in WRITEUP.md,
section 6: a header 'p cert-degsens n d', then the tree in preorder, 'b M' branching on g_M (the g_M = 0 subtree
first) and 'l k r_1 y_1 ... r_k y_k' a leaf, where row r is the r-th set S with |S| > d in increasing bitmask order
(bit i-1 of a mask is coordinate i).

The statement.  For g : {0,1}^n -> {0,1} write g_T = g(1_T).  Its multilinear polynomial has the Moebius
coefficients a_S = sum_{T subset S} (-1)^{|S|-|T|} g_T, so deg g <= d iff a_S = 0 for every |S| > d.  A g with
deg g <= d, g(0) = 0 and g(e_i) = 1 for all i exists iff the system
    A g = 0,  g in {0,1}^(2^n),  g_0 = 0,  g_{i} = 1 (1 <= i <= n)                                       (*)
has a solution (A: one row per |S| > d, entries (-1)^{|S|-|T|} for T subset S).  Re-derived reduction (the PR's
Lemma 1): if deg f <= d and s(f) >= n, restrict f to n coordinates sensitive at a point x, shift x to 0 by y -> 1 - y
where x_i = 1, and complement if g(0) = 1; the degree does not rise.  Restricting to x_11 = 0 maps an n = 11
solution to an n = 10 one, so infeasibility of (*) for (10, 4) gives s(f) <= 9 for every f with deg f <= 4.

Our check of a certificate.  We rebuild A from (n, d) alone.  For each node we keep the box B: 0 <= g <= 1 with the
root fixings and the fixings on the path.  At a leaf with multipliers y, the vector c = y^T A is computed exactly
(c_T = (-1)^{|T|} sum_{S superset T, |S| > d} (-1)^{|S|} y_S, by a superset-sum transform in integers) and we check
max_{g in B} c.g < 0, a sum over coordinates.  Every solution of (*) has c.g = 0, so none lies in B.  The two children
of a branch fix g_M = 0 and g_M = 1, so the leaf boxes cover every 0/1 point with the root fixings: (*) has no
solution.

Also: Proskurin's symmetric point g_T = p(|T|), p = 0, 1, 1, 7/12, 1/6, 0, 1/6, 7/12, 1, 1, 0, solves the LP
relaxation of (*) for (10, 4) exactly (so integrality is needed), and NAE_3 o NAE_3 on 9 variables has degree 4 and
sensitivity 9 (the statement is sharp), both from our own Moebius transform and sensitivity count.

usage: python3 -I c37_check.py CERTS_DIR [--controls]
"""
import hashlib
import lzma
import os
import sys
import time
from fractions import Fraction as F

SHA = {"d10_4.cert.xz": "3be6fa8b03d544367b3b61de058de4784885aa20e7e26b8b95813e35dfe7e9f4",
       "d10_4.cert": "7dc38b37e9af02340fe0282b523452968e6d0e72320eff4c8840d9a67b4e7e6c",
       "d7_3.cert": "f968ea647a26deef6cbcc10c8ec87a301742c4fa9659047ea8a50723381c2ab6",
       "d4_2.cert": "f58ba4cf223546bcd806d734ac3b2385153728e1959021fb7335b1e7e38cf123"}


def popcount(x):
    return bin(x).count("1")


KERNEL = {}   # (n, d) -> vectors h with A h = 0, used to test every leaf vector c (c.h must be 0)


def rows_of(n, d):
    return [S for S in range(1 << n) if popcount(S) > d]


def leaf_vector(n, rows, mult):
    """c = y^T A: c_T = (-1)^{|T|} sum_{S >= T} (-1)^{|S|} y_S (integers)."""
    N = 1 << n
    w = [0] * N
    for r, y in mult:
        S = rows[r]
        w[S] += -y if popcount(S) & 1 else y
    for i in range(n):          # superset sums: w[T] <- sum over S superset of T
        bit = 1 << i
        for T in range(N):
            if not T & bit:
                w[T] += w[T | bit]
    return [(-v if popcount(T) & 1 else v) for T, v in enumerate(w)]


def direct_row_combo(n, rows, mult, T):
    """the same coefficient by the definition, for spot checks."""
    s = 0
    for r, y in mult:
        S = rows[r]
        if S & T == T:
            s += y * (-1) ** (popcount(S) - popcount(T))
    return s


def parse(text):
    lines = [ln.split() for ln in text.strip().split("\n") if ln.strip()]
    if lines[0][:2] != ["p", "cert-degsens"]:
        raise ValueError("bad header")
    n, d = int(lines[0][2]), int(lines[0][3])
    return n, d, lines[1:]


def check_tree(n, d, body, run_spot=0):
    """walks the preorder tree; returns (ok, leaves, depth, message)."""
    rows = rows_of(n, d)
    N = 1 << n
    lo = [0] * N
    hi = [1] * N
    lo[0] = hi[0] = 0                     # g(0) = 0
    for i in range(n):
        lo[1 << i] = hi[1 << i] = 1       # g(e_i) = 1
    pos = 0
    leaves = 0
    maxdepth = 0
    spot = [0]

    def node(depth):
        nonlocal pos, leaves, maxdepth
        if pos >= len(body):
            raise ValueError("the tree ends early")
        tok = body[pos]
        pos += 1
        maxdepth = max(maxdepth, depth)
        if tok[0] == "b":
            M = int(tok[1])
            if not (0 <= M < N) or lo[M] != 0 or hi[M] != 1:
                raise ValueError("branch on a fixed or invalid variable %d" % M)
            hi[M] = 0
            ok0 = node(depth + 1)
            hi[M] = 1
            lo[M] = 1
            ok1 = node(depth + 1)
            lo[M] = 0
            return ok0 and ok1
        if tok[0] == "l":
            k = int(tok[1])
            vals = list(map(int, tok[2:]))
            if len(vals) != 2 * k:
                raise ValueError("leaf length")
            mult = list(zip(vals[0::2], vals[1::2]))
            if any(not (0 <= r < len(rows)) for r, _ in mult) or len({r for r, _ in mult}) != k:
                raise ValueError("leaf rows")
            c = leaf_vector(n, rows, mult)
            if spot[0] < run_spot:      # spot check the transform against the definition
                for T in (0, 1, N - 1, (N - 1) // 3, N // 2 + 5):
                    if c[T] != direct_row_combo(n, rows, mult, T):
                        raise ValueError("transform mismatch")
                spot[0] += 1
            leaves += 1
            for h in KERNEL.get((n, d), ()):
                if sum(ci * hv for ci, hv in zip(c, h)) != 0:
                    raise ValueError("c.h != 0 at a solution of A g = 0")
            m = sum(ci * (hi[T] if ci > 0 else lo[T]) for T, ci in enumerate(c))
            return m < 0
        raise ValueError("unknown line %r" % tok[0])

    ok = node(0)
    if pos != len(body):
        raise ValueError("lines after the tree")
    return ok, leaves, maxdepth


def moebius(table, n):
    a = list(table)
    for i in range(n):
        bit = 1 << i
        for T in range(1 << n):
            if T & bit:
                a[T] -= a[T ^ bit]
    return a


def nae3(x, y, z):
    return 0 if x == y == z else 1


class Run:
    def __init__(self):
        self.n = 0; self.fails = []

    def __call__(self, label, ok):
        self.n += 1
        if not ok:
            self.fails.append(label)
        print(f"  check {label}: {'PASS' if ok else 'FAIL'}")
        return ok


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    if len(args) != 1:
        print(__doc__); sys.exit(2)
    t0 = time.time()
    root = args[0]
    run = Run()
    print("37a: PR #283's statement, no Boolean function of degree 4 has sensitivity 10 (n = 10, d = 4)")
    raw_xz = open(os.path.join(root, "d10_4.cert.xz"), "rb").read()
    raw = lzma.decompress(raw_xz)
    run("certificate sha256: d10_4.cert.xz = 3be6fa8b...e7e9f4, decompressed = 7dc38b37...7e6c (as the PR states)",
        hashlib.sha256(raw_xz).hexdigest() == SHA["d10_4.cert.xz"] and hashlib.sha256(raw).hexdigest() == SHA["d10_4.cert"])
    texts = {"d10_4": raw.decode("ascii")}
    for name in ("d7_3", "d4_2"):
        b = open(os.path.join(root, name + ".cert"), "rb").read()
        run("certificate sha256: %s.cert = %s..." % (name, SHA[name + ".cert"][:8]),
            hashlib.sha256(b).hexdigest() == SHA[name + ".cert"])
        texts[name] = b.decode("ascii")
    # two solutions of A g = 0 for (10, 4), to test every leaf vector: Proskurin's symmetric point (times 12) and
    # NAE_3 o NAE_3 on the first 9 coordinates (it ignores the tenth, so g(e_10) = 0 and it is not a solution of (*))
    pp = [0, 12, 12, 7, 2, 0, 2, 7, 12, 12, 0]
    nn = []
    for x in range(1 << 10):
        b = [(x >> i) & 1 for i in range(9)]
        nn.append(nae3(nae3(*b[0:3]), nae3(*b[3:6]), nae3(*b[6:9])))
    KERNEL[(10, 4)] = ([pp[popcount(T)] for T in range(1 << 10)], nn)
    run("the two test vectors solve A g = 0 for (10, 4) (all 638 Moebius coefficients of degree > 4 vanish)",
        all(all(moebius(h, 10)[S] == 0 for S in rows_of(10, 4)) for h in KERNEL[(10, 4)]))
    for name, nd in (("d10_4", (10, 4)), ("d7_3", (7, 3)), ("d4_2", (4, 2))):
        n, d, body = parse(texts[name])
        nb = sum(1 for t in body if t[0] == "b")
        ok, leaves, depth = check_tree(n, d, body, run_spot=40)
        nz = max(int(t[1]) for t in body if t[0] == "l")
        run("(%d, %d): %d rows |S| > %d rebuilt from (n, d); %d branchings, %d leaves, depth %d; every leaf box "
            "infeasible (max c.g < 0, exact integers; up to %d multipliers per leaf%s); the leaves cover all 0/1 "
            "points: no g with deg g <= %d, g(0) = 0, g(e_i) = 1"
            % (n, d, len(rows_of(n, d)), d, nb, leaves, depth, nz,
               "; c.h = 0 at both test vectors at every leaf" if (n, d) in KERNEL else "", d),
            (n, d) == nd and ok and leaves == nb + 1)
    # LP relaxation is feasible: Proskurin's symmetric point
    p = [F(0), F(1), F(1), F(7, 12), F(1, 6), F(0), F(1, 6), F(7, 12), F(1), F(1), F(0)]
    g = [p[popcount(T)] for T in range(1 << 10)]
    a = moebius(g, 10)
    run("Proskurin's symmetric point solves the LP relaxation of (*) for (10, 4) exactly (all 638 a_S = 0, bounds and "
        "fixings hold), so the proof must use integrality",
        all(a[S] == 0 for S in rows_of(10, 4)) and g[0] == 0 and all(g[1 << i] == 1 for i in range(10))
        and all(0 <= v <= 1 for v in g))
    # sharpness: NAE3 o NAE3 on 9 variables
    tab = []
    for x in range(1 << 9):
        b = [(x >> i) & 1 for i in range(9)]
        tab.append(nae3(nae3(*b[0:3]), nae3(*b[3:6]), nae3(*b[6:9])))
    am = moebius(tab, 9)
    deg = max(popcount(S) for S in range(1 << 9) if am[S] != 0)
    sens = max(sum(tab[x] != tab[x ^ (1 << i)] for i in range(9)) for x in range(1 << 9))
    run("NAE_3 o NAE_3 on 9 variables: degree %d, sensitivity %d (the statement is sharp)" % (deg, sens),
        deg == 4 and sens == 9)
    if "--controls" in sys.argv:
        print("  -- red controls (each must be refused)")
        n, d, body = parse(texts["d10_4"])
        i = next(j for j, t in enumerate(body) if t[0] == "l")
        bad = [list(t) for t in body]
        bad[i][3] = str(-int(bad[i][3]))
        ok, _, _ = check_tree(n, d, bad)
        run("control: the first leaf's first multiplier negated", not ok)
        bad = [list(t) for t in body]
        j = max(jj for jj, t in enumerate(body) if t[0] == "l")
        del bad[j]
        try:
            check_tree(n, d, bad); refused = False
        except ValueError:
            refused = True
        run("control: the last leaf deleted (the tree is incomplete)", refused)
        bad = [list(t) for t in body]
        k = next(jj for jj, t in enumerate(body) if t[0] == "b")
        bad[k][1] = "0"
        try:
            check_tree(n, d, bad); refused = False
        except ValueError:
            refused = True
        run("control: the first branch moved to the fixed variable g_0", refused)
    print(f"  checks: {run.n}, failed: {len(run.fails)}")
    print(f"  run time: {time.time() - t0:.1f} s")
    print("RESULT: ALL CHECKS PASS" if not run.fails else "RESULT: FAIL " + "; ".join(run.fails))
    sys.exit(0 if not run.fails else 1)


if __name__ == "__main__":
    main()
