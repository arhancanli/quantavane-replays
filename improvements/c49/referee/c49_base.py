"""Base data for the 49a 54-point block, rebuilt from the text of Tilley's note (page 2) only.

Everything is integers / bitmasks. No file from the archive is read.
"""
import itertools

# ---- K5 edges, lexicographic order, bit 0 = edge 01 -------------------------------
EDGES = [(u, v) for u in range(5) for v in range(u + 1, 5)]
assert EDGES == [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
EIDX = {e: i for i, e in enumerate(EDGES)}


def ebit(u, v):
    if u > v:
        u, v = v, u
    return 1 << EIDX[(u, v)]


def emask(pairs):
    m = 0
    for (u, v) in pairs:
        m |= ebit(u, v)
    return m


def eset(m):
    return ["%d%d" % EDGES[i] for i in range(10) if m >> i & 1]


PORT_P = ebit(0, 1)   # p = 01
PORT_Q = ebit(0, 2)   # q = 02
E1 = ebit(0, 3)       # first pendant edge
E2 = ebit(0, 4)       # second pendant edge
ALL_E = (1 << 10) - 1


def permute_mask(m, perm):
    """Apply a vertex permutation (dict or list v->perm[v]) to an edge mask."""
    out = 0
    for i, (u, v) in enumerate(EDGES):
        if m >> i & 1:
            out |= ebit(perm[u], perm[v])
    return out


# ---- F: stars, K_{2,3}, complements of chairs --------------------------------------
def build_F():
    stars = []
    for v in range(5):
        stars.append(emask([(v, w) for w in range(5) if w != v]))
    k23 = []
    for a, b in itertools.combinations(range(5), 2):
        rest = [w for w in range(5) if w not in (a, b)]
        k23.append(emask([(x, y) for x in (a, b) for y in rest]))
    chairs = []
    for comb in itertools.combinations(range(10), 4):
        es = [EDGES[i] for i in comb]
        deg = [0] * 5
        for (u, v) in es:
            deg[u] += 1
            deg[v] += 1
        # connected?
        adj = {v: set() for v in range(5)}
        for (u, v) in es:
            adj[u].add(v)
            adj[v].add(u)
        seen = {0}
        stack = [0]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        if len(seen) == 5 and sorted(deg) == [1, 1, 1, 2, 3]:
            chairs.append(emask(es))
    assert emask([(0, 1), (0, 2), (0, 3), (1, 4)]) in chairs, "example chair missing"
    cchairs = [ALL_E ^ c for c in chairs]
    F = stars + k23 + cchairs
    assert len(stars) == 5 and len(k23) == 10 and len(chairs) == 60
    assert len(set(F)) == 75
    return F


# ---- H on Z_11 ---------------------------------------------------------------------
def build_H():
    QR = sorted({(x * x) % 11 for x in range(1, 11)})
    assert QR == [1, 3, 4, 5, 9]
    bases = [{0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 4, 6, 7}, {0, 1, 2, 6, 9}]
    H = set()
    for S in bases:
        for a in QR:
            for b in range(11):
                m = 0
                for s in S:
                    m |= 1 << ((a * s + b) % 11)
                H.add(m)
    H = sorted(H)
    assert len(H) == 121, len(H)
    return H


A1 = [251, 425, 499, 506, 570, 707, 761, 875, 883, 889, 923, 947, 985, 986, 1002]
A2 = [358, 471, 499, 501, 631, 707, 757, 758, 789, 883, 935, 947, 950, 982, 997]

SWAP_13 = [0, 3, 2, 1, 4]   # Ap: core vertices 1 and 3 interchanged
SWAP_23 = [0, 1, 3, 2, 4]   # Aq: core vertices 2 and 3 interchanged

# X*: core edges 12, 14, 23 + all ten points of the SECOND pendant; nothing in pendant 1; no port.
XSTAR_CORE = emask([(1, 2), (1, 4), (2, 3)])
PEND_FULL = (1 << 10) - 1
