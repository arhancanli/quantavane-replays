#!/usr/bin/env python3
"""Independent exact replay of the asterisked lower bounds for C_49 (teorth/optimizationproblems, constant 49a),
the Erdos-Szemeredi 3-sunflower-free capacity mu_3 = lim f(n)^(1/n).

Written from the definitions on the constant page and the claimants' data files, READMEs and format notes. The
claimants' verifiers were neither read nor run. Standard library only (Python 3.9+), run with python3 -I.

  sunflower_check.py ROW CERTS_ROOT [--controls] [--workers=N]   ROW in: irsi_v10 irsi_v11 irsi_v12 irsi_degkm cla_F
                                                                  cla_FH syl_unit syl_Bprime kos lemma selftest
  (--controls adds the red controls; --workers sets the number of forked processes, default: all cores)

What a row replay establishes (each step is exact: Python integers and fractions, no floating point decides anything)
1. the certificate file's SHA-256 (fetch.sh pins it; here it is compared with the fragments stated on the page / PRs);
2. the base family: members distinct, Sperner (no member strictly inside another) and 3-sunflower-free (no three
   distinct members A, B, C with A&B == A&C == B&C), decided exhaustively by one or two exact methods:
     pivot  for each member A (in index order) the later members S are grouped by S & A; a sunflower with A as its
            first member exists iff one group holds two members B, C with B & C == the group key. Every triple is
            examined from its first member, so this covers all triples. Containment is tested in the same pass.
     pair   for every pair A, B the admissible third members are exactly those C with A&B <= C and C & (A^B) == 0;
            they are found as an AND of per-point member bitsets. A second, independent generic method.
     box    the family is split exactly into boxes S_0 x S_1 x ... (products of block projections over a partition of
            the points; the split is verified to rebuild the family). Three members form a sunflower iff on every
            block their projections satisfy x&y == x&z == y&z (repetitions allowed) and every pair of members
            differs on some block; as the blocks of a box are independent, it suffices to know, per block and per
            triple of sides, which equality patterns occur. Exhaustive over all triples; fast for glued families.
     brute  every triple and pair, straight from the definition (small families only).
   Every reported witness is re-checked against the definition, and every method is run against planted defects.
   For the two 250k-set families (the [Syl2026] gadget G and the [Kos2026] block) the exhaustive method is the box
   method on two different point partitions; the generic methods would need about 10 CPU-hours there. As a second,
   independent path the [Kos2026] block is shown equal to its rebuild from the description, its 250659 product sets
   are valid by the gluing lemma, and every pair and triple through its 1633 added sets is scanned exhaustively.
3. the port counts N_ab = #{S : [u in S] = a, [v in S] = b} and the 2x2 transfer matrix of the chain;
4. the Perron root rho = (tr + sqrt((a-d)^2 + 4bc)) / 2 as an exact quadratic surd and the bound rho^(1/L), enclosed
   between consecutive decimals by exact rational comparisons q^L <= rho (no rounding involved);
5. the printed bound b is compared exactly with rho^(1/L): a lower bound must satisfy b^L <= rho.

The chain argument (gluing lemma of Deuber-Erdos-Gunderson-Kostochka-Meyer 1997, Lemma 4.4, and the transfer-matrix
count) is mathematics, not computation; `lemma` tests it on all small instances and on explicit glued chains.
"""
import hashlib
import itertools
import json
import lzma
import os
import random
import sys
import time
from collections import Counter
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import lru_cache
from math import comb, isqrt


def popcount(x):
    return bin(x).count("1")


def bits(x):
    """Indices of the set bits of x."""
    out = []
    while x:
        low = x & -x
        out.append(low.bit_length() - 1)
        x ^= low
    return out


def mask(points):
    m = 0
    for p in points:
        m |= 1 << p
    return m


# ----------------------------------------------------------------------------------------------------------------
# Reading the certificates (as data)
# ----------------------------------------------------------------------------------------------------------------

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_irsi(path):
    """c49-chain-v1 JSON: bases {name: {n, sets: [bitmasks]}}, word [{base, u, v}]."""
    d = json.load(open(path))
    assert d.get("format") == "c49-chain-v1", "unexpected format"
    bases = {}
    for name, b in d["bases"].items():
        n = int(b["n"])
        sets = [int(s) for s in b["sets"]]
        assert all(isinstance(s, int) and 0 <= s < (1 << n) for s in b["sets"]), "member outside the ground set"
        bases[name] = (n, sets)
    word = [(w["base"], int(w["u"]), int(w["v"])) for w in d["word"]]
    return d, bases, word


def read_mask_list(path):
    """First line: count; then whitespace-separated integer bitmasks."""
    toks = open(path).read().split()
    m = int(toks[0])
    sets = [int(t) for t in toks[1:]]
    assert len(sets) == m, "count line disagrees with the number of masks"
    return sets


def read_point_lines(path):
    """One member per line as point indices; lines starting with # are comments."""
    sets = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        sets.append(mask(int(t) for t in line.split()))
    return sets


def read_sunflower3(path):
    """The sunflower3 v1 format (FORMAT.md of [Kos2026]), parsed strictly. Returns [(n, p, q, [masks])]."""
    raw = lzma.open(path).read() if path.endswith(".xz") else open(path, "rb").read()
    text = raw.decode("ascii")
    assert "\r" not in text, "CR characters are not allowed"
    lines = [ln for ln in text.split("\n") if not ln.startswith("#")]
    if lines and lines[-1] == "":
        lines.pop()
    assert lines[0] == "sunflower3 v1", "bad header"

    def num(t):
        assert t.isdigit() and (t == "0" or t[0] != "0"), "bad token %r" % t
        return int(t)

    blocks, i = [], 1
    while i < len(lines):
        head = lines[i].split(" ")
        assert head[0] == "block" and len(head) == 5, "bad block line %r" % lines[i]
        n, p, q, m = (num(t) for t in head[1:])
        assert 3 <= n <= 64 and p < n and q < n and p != q
        members = []
        for ln in lines[i + 1:i + 1 + m]:
            if ln == "-":
                members.append(0)
                continue
            pts = [num(t) for t in ln.split(" ")]
            assert all(a < b for a, b in zip(pts, pts[1:])) and all(x < n for x in pts), "bad member line"
            members.append(mask(pts))
        assert len(members) == m, "block shorter than announced"
        blocks.append((n, p, q, members))
        i += 1 + m
    return blocks


# ----------------------------------------------------------------------------------------------------------------
# The families of [Cla2026], rebuilt from the words on the page, and the gluing of [DEGKM1997] Lemma 4.4
# ----------------------------------------------------------------------------------------------------------------

K5_EDGES = [(a, b) for a in range(5) for b in range(a + 1, 5)]  # lexicographic: (0,1),(0,2),...,(3,4)


def family_F():
    """75 edge sets of K5: the 5 stars K_{1,4}, the 10 graphs K_{2,3}, the 60 graphs whose complement is a chair
    (degree sequence 3,2,1,1,1: a centre c joined to three vertices, one of which, a, is joined to the fifth, e).
    Members are frozensets of edges (a, b) with a < b."""
    E = set(K5_EDGES)
    e = lambda x, y: (min(x, y), max(x, y))
    stars = {frozenset(e(v, w) for w in range(5) if w != v) for v in range(5)}
    k23 = {frozenset(e(x, y) for x in X for y in range(5) if y not in X) for X in itertools.combinations(range(5), 2)}
    chairs = set()
    for c in range(5):
        for ee in range(5):
            if ee == c:
                continue
            nbrs = [w for w in range(5) if w not in (c, ee)]
            for a in nbrs:
                chair = {e(c, w) for w in nbrs} | {e(a, ee)}
                chairs.add(frozenset(E - chair))
    assert len(stars) == 5 and len(k23) == 10 and len(chairs) == 60
    fam = stars | k23 | chairs
    assert len(fam) == 75
    return sorted(fam, key=lambda s: sorted(s))


def family_H():
    """121 subsets of Z_11: images of {0..6}, {0,1,2,3,4,6,7}, {0,1,2,6,9} under x -> a x + b, a a nonzero square."""
    squares = sorted({(x * x) % 11 for x in range(1, 11)})
    assert squares == [1, 3, 4, 5, 9]
    out = set()
    for base in ({0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 4, 6, 7}, {0, 1, 2, 6, 9}):
        for a in squares:
            for b in range(11):
                out.add(frozenset((a * x + b) % 11 for x in base))
    assert len(out) == 121
    return sorted(out, key=lambda s: sorted(s))


def relabel(fam, tag):
    """Members as frozensets of labels (tag, point)."""
    return [frozenset((tag, p) for p in S) for S in fam]


def glue(A, a, B, b):
    """[DEGKM1997] Lemma 4.4 as stated in the 49a references: {(C - a) | D : a in C in A, b not in D in B}
    union {C | (D - b) : a not in C in A, b in D in B}. A, B: lists of frozensets on disjoint label sets."""
    A1 = [C - {a} for C in A if a in C]
    A0 = [C for C in A if a not in C]
    B0 = [D for D in B if b not in D]
    B1 = [D - {b} for D in B if b in D]
    return [P | Q for P in A1 for Q in B0] + [P | Q for P in A0 for Q in B1]


def to_masks(fam, order):
    """Frozensets of labels -> bitmasks, with point i = order[i]."""
    index = {lab: i for i, lab in enumerate(order)}
    return [mask(index[x] for x in S) for S in fam]


def ground(fam):
    out = set()
    for S in fam:
        out |= S
    return out


# ----------------------------------------------------------------------------------------------------------------
# Exhaustive checks. Each returns (containment, sunflower): None, or a witness tuple of member masks.
# ----------------------------------------------------------------------------------------------------------------

def is_sunflower(a, b, c):
    return a != b and a != c and b != c and (a & b) == (a & c) == (b & c)


def brute_check(fam):
    """Every pair and triple, straight from the definition."""
    contain = sun = None
    for a, b in itertools.combinations(fam, 2):
        if (a & b) in (a, b):
            contain = contain or (a, b)
    for a, b, c in itertools.combinations(fam, 3):
        if (a & b) == (a & c) == (b & c):
            sun = (a, b, c)
            break
    return contain, sun


def group_pair(g, k):
    """Two members s != t of the group g (all with s & A == k, so s >= k) with s & t == k, or None.
    Small groups: every pair, scanned at C level. Large groups: per point p, the bitset of members whose outer part
    s ^ k avoids p; the partners of s are the AND of those bitsets over the points of its outer part."""
    m = len(g)
    if m < 2:
        return None
    if m == 2:
        return (g[0], g[1]) if g[0] & g[1] == k else None
    if m <= 40:
        if k in itertools.starmap(int.__and__, itertools.combinations(g, 2)):
            return next(p for p in itertools.combinations(g, 2) if p[0] & p[1] == k)
        return None
    full = (1 << m) - 1
    avoid = {}
    outs = [s ^ k for s in g]
    for j, o in enumerate(outs):
        if o == 0:  # s == k: with any other member t, s & t == k
            return g[j], g[1 - j if j < 2 else 0]
        for p in bits(o):
            avoid[p] = avoid.get(p, full) & ~(1 << j)
    for j, o in enumerate(outs):
        cand = full
        for p in bits(o):
            cand &= avoid[p]
            if not cand:
                break
        if cand:
            return g[j], g[(cand & -cand).bit_length() - 1]
    return None


def pivot_range(fam, lo, hi, step, everyone=False):
    """Pivots i in range(lo, hi, step). The later members (with everyone=True: all other members) are grouped by
    S & A; two members B, C of one group with B & C == key give the sunflower {A, B, C}.
    Returns (containment, sunflower, pivots scanned)."""
    contain = sun = None
    done = 0
    for i in range(lo, hi, step):
        a = fam[i]
        rest = fam[:i] + fam[i + 1:] if everyone else fam[i + 1:]
        keys = list(map(a.__and__, rest))
        done += 1
        if contain is None:
            if a in keys:  # a <= s
                contain = (a, rest[keys.index(a)])
            elif any(map(int.__eq__, rest, keys)):  # s <= a
                contain = (next(s for s, k in zip(rest, keys) if s == k), a)
        if sun:  # a sunflower is known: finish the containment scan only
            if contain:
                break
            continue
        groups = {}
        for k, s in zip(keys, rest):
            g = groups.get(k)
            if g is None:
                groups[k] = [s]
            else:
                g.append(s)
        for k, g in groups.items():
            if len(g) > 1:
                hit = group_pair(g, k)
                if hit:
                    sun = (a,) + hit
                    break
        if sun and contain:
            break
    return contain, sun, done


_SHARED = {}


def _worker(args):
    kind, lo, hi, step = args
    fam = _SHARED["fam"]
    if kind == "pivot":
        c, s, _ = pivot_range(fam, lo, hi, step)
    else:
        c, s = pair_range(fam, _SHARED["has"], _SHARED["hasnot"], lo, hi, step)
    return c, s


def _parallel(kind, fam, workers):
    """Deal the first index i round-robin to forked worker processes (copy-on-write: nothing is pickled)."""
    import multiprocessing as mp
    _SHARED["fam"] = fam
    if kind == "pair":
        _SHARED["has"], _SHARED["hasnot"] = point_bitsets(fam)
    ctx = mp.get_context("fork")
    chunks = workers * 8
    jobs = [(kind, r, len(fam), chunks) for r in range(chunks)]
    contain = sun = None
    with ctx.Pool(workers) as pool:
        for c, s in pool.imap_unordered(_worker, jobs):
            contain = contain or c
            sun = sun or s
    _SHARED.clear()
    return contain, sun


def pivot_check(fam, workers=1):
    """The pivot method over all members (every triple from its first member)."""
    if workers <= 1 or len(fam) < 3000:
        c, s, _ = pivot_range(fam, 0, len(fam), 1)
        return c, s
    return _parallel("pivot", fam, workers)


def point_bitsets(fam):
    n = len(fam)
    width = max(fam).bit_length() if fam else 0
    full = (1 << n) - 1
    has = [0] * width
    for i, s in enumerate(fam):
        for p in bits(s):
            has[p] |= 1 << i
    return has, [full ^ h for h in has]


def pair_range(fam, has, hasnot, lo, hi, step):
    """For every pair (A, B) with A first: the admissible third members are the C with A&B <= C and C & (A^B) == 0,
    found by ANDing per-point bitsets over the member indices; A and B themselves are removed from the candidates."""
    n = len(fam)
    full = (1 << n) - 1
    contain = sun = None
    for i in range(lo, hi, step):
        a = fam[i]
        for j in range(i + 1, n):
            b = fam[j]
            k = a & b
            if k == a or k == b:
                contain = contain or ((a, b) if k == a else (b, a))
                if sun:
                    return contain, sun
                continue
            if sun:  # a sunflower is known: finish the containment scan only
                continue
            cand = full ^ ((1 << i) | (1 << j))
            for p in bits(k):
                cand &= has[p]
                if not cand:
                    break
            else:
                for p in bits(a ^ b):
                    cand &= hasnot[p]
                    if not cand:
                        break
            if cand:
                sun = (a, b, fam[(cand & -cand).bit_length() - 1])
                if contain:
                    return contain, sun
    return contain, sun


def pair_check(fam, workers=1):
    """Independent of the pivot method; skipping comparable pairs loses nothing, since in a sunflower at most one
    member is comparable with the others (a member equal to the kernel), and the other two are incomparable."""
    if workers <= 1 or len(fam) < 1500:
        has, hasnot = point_bitsets(fam)
        return pair_range(fam, has, hasnot, 0, len(fam), 1)
    return _parallel("pair", fam, workers)


# ---- the box method ---------------------------------------------------------------------------------------------

def decompose(tuples, t):
    """Split a set of t-tuples exactly into boxes S_0 x ... x S_{t-1} (tuples of frozensets): group by the last
    coordinate, split each fibre recursively, and merge the last coordinates that share a sub-box."""
    if t == 1:
        return [(frozenset(x[0] for x in tuples),)]
    by_last = {}
    for x in tuples:
        by_last.setdefault(x[-1], set()).add(x[:-1])
    merged = {}
    for v, prefixes in by_last.items():
        for sub in decompose(prefixes, t - 1):
            merged.setdefault(sub, set()).add(v)
    return [sub + (frozenset(vs),) for sub, vs in merged.items()]


# separation pattern of a triple (x, y, z): bit 0 = x != y, bit 1 = x != z, bit 2 = y != z
def _sep(x, y, z):
    return (x != y) | ((x != z) << 1) | ((y != z) << 2)


def ws_patterns_naive(X, Y, Z):
    """8-bit mask: bit s is set iff some (x, y, z) in X x Y x Z has x&y == x&z == y&z and separation pattern s.
    Reference version, cubic in the side sizes (the self-test compares it with ws_patterns)."""
    Zl = list(Z)
    res = 0
    for x in X:
        for y in Y:
            k = x & y
            d = x ^ y
            for z in Zl:
                if not (z & d) and (z & k) == k:
                    res |= 1 << _sep(x, y, z)
    return res


def ws_patterns(X, Y, Z):
    """Same mask as ws_patterns_naive. For each pair (x, y), the admissible z (k = x&y inside z, z disjoint from x^y)
    are found as an AND of per-point bitsets over Z; then only z == x, z == y or 'another z' matter."""
    Zl = list(Z)
    pos = {z: j for j, z in enumerate(Zl)}
    full = (1 << len(Zl)) - 1
    has = {}
    for j, z in enumerate(Zl):
        for p in bits(z):
            has[p] = has.get(p, 0) | (1 << j)
    hasnot = {p: full & ~h for p, h in has.items()}
    res = 0
    for x in X:
        bit_x = 1 << pos[x] if x in pos else 0
        for y in Y:
            cand = full
            for p in bits(x & y):
                cand &= has.get(p, 0)
                if not cand:
                    break
            if cand:
                for p in bits(x ^ y):
                    cand &= hasnot.get(p, full)
                    if not cand:
                        break
            if not cand:
                continue
            bit_y = 1 << pos[y] if y in pos else 0
            neq = int(x != y)
            if cand & ~(bit_x | bit_y):
                res |= 1 << (neq | 6)  # z differs from x and from y
            if cand & bit_x:
                res |= 1 << (neq | (neq << 2))  # z == x
            if cand & bit_y:
                res |= 1 << (neq | (neq << 1))  # z == y
    return res


def ws_example(X, Y, Z, s):
    for x in X:
        for y in Y:
            k, d = x & y, x ^ y
            for z in Z:
                if not (z & d) and (z & k) == k and _sep(x, y, z) == s:
                    return x, y, z
    return None


def sub_patterns(X, Y):
    """bit 0: some x in X equals some y in Y; bit 1: some x in X is a proper subset of some y in Y."""
    res = 0
    for x in X:
        for y in Y:
            if x & y == x:
                res |= 1 if x == y else 2
                if res == 3:
                    return 3
    return res


def box_check(fam, blocks, boxes=None):
    """Exhaustive Sperner and 3-sunflower check by the box method. blocks: list of disjoint point masks covering
    every point used. Returns (containment, sunflower, info). `boxes` may be given (tests only); either way the
    split is verified to rebuild the family exactly before it is used."""
    t = len(blocks)
    allpts = 0
    for bm in blocks:
        assert allpts & bm == 0, "blocks overlap"
        allpts |= bm
    famset = set(fam)
    assert len(famset) == len(fam), "members are not distinct"
    assert all(s & ~allpts == 0 for s in fam), "a member uses a point outside the blocks"
    if boxes is None:
        boxes = decompose({tuple(s & bm for bm in blocks) for s in fam}, t)
    assert all(len(box) == t and all(x & ~bm == 0 for X, bm in zip(box, blocks) for x in X) for box in boxes), \
        "a box side leaves its block"
    # the split must rebuild the family exactly (no member twice, none missing, nothing extra)
    seen = set()
    total = 0
    for box in boxes:
        for parts in itertools.product(*box):
            s = 0
            for p in parts:
                s |= p
            assert s in famset and s not in seen, "box decomposition is not exact"
            seen.add(s)
            total += 1
    assert total == len(fam) == len(seen), "box decomposition is not exact"
    # side ids per block
    side_ids = [{} for _ in range(t)]
    for box in boxes:
        for i in range(t):
            side_ids[i].setdefault(box[i], len(side_ids[i]))
    sides = [[None] * len(side_ids[i]) for i in range(t)]
    for i in range(t):
        for S, k in side_ids[i].items():
            sides[i][k] = S
    bx = [tuple(side_ids[i][box[i]] for i in range(t)) for box in boxes]

    @lru_cache(maxsize=None)
    def ws(i, a, b, c):
        return ws_patterns(sides[i][a], sides[i][b], sides[i][c])

    @lru_cache(maxsize=None)
    def sp(i, a, b):
        return sub_patterns(sides[i][a], sides[i][b])

    @lru_cache(maxsize=None)
    def covers(pats):
        """Is there one separation pattern per block whose union separates all three pairs (value 7)?"""
        reach = 1  # bitmask over OR-values 0..7, starting from {0}
        for m in pats:
            new = 0
            for r in range(8):
                if reach >> r & 1:
                    for s in range(8):
                        if m >> s & 1:
                            new |= 1 << (r | s)
            reach = new
        return bool(reach >> 7 & 1)

    contain = sun = None
    nb = len(bx)
    # containment: ordered pairs of boxes
    for p in range(nb):
        for q in range(nb):
            rel = [sp(i, bx[p][i], bx[q][i]) for i in range(t)]
            if all(rel) and any(r & 2 for r in rel):
                contain = _contain_witness(boxes[p], boxes[q])
                break
        if contain:
            break
    # sunflowers: multisets of three boxes
    for p in range(nb):
        for q in range(p, nb):
            for r in range(q, nb):
                pats = []
                for i in range(t):
                    m = ws(i, bx[p][i], bx[q][i], bx[r][i])
                    if not m:
                        break
                    pats.append(m)
                else:
                    if covers(tuple(pats)):
                        sun = _sun_witness(boxes[p], boxes[q], boxes[r], pats)
                        break
            if sun:
                break
        if sun:
            break
    info = {"boxes": nb, "sides per block": [len(s) for s in sides],
            "box triples": comb(nb + 2, 3), "box pairs": nb * nb}
    return contain, sun, info


def _contain_witness(P, Q):
    """A proper containment between the boxes P and Q: a strict pair on the first block that has one, any
    included pair on the others."""
    small = big = 0
    used_strict = False
    for X, Y in zip(P, Q):
        strict = next(((x, y) for x in X for y in Y if x & y == x and x != y), None)
        x, y = strict if (strict and not used_strict) else next((x, y) for x in X for y in Y if x & y == x)
        used_strict = used_strict or x != y
        small, big = small | x, big | y
    return small, big


def _sun_witness(P, Q, R, pats):
    for combo in itertools.product(*[[s for s in range(8) if m >> s & 1] for m in pats]):
        acc = 0
        for s in combo:
            acc |= s
        if acc == 7:
            a = b = c = 0
            for X, Y, Z, s in zip(P, Q, R, combo):
                x, y, z = ws_example(X, Y, Z, s)
                a, b, c = a | x, b | y, c | z
            return a, b, c
    raise AssertionError("pattern bookkeeping error")


def verify_witnesses(fam, contain, sun):
    """Re-check reported witnesses against the definitions."""
    fs = set(fam)
    if contain is not None:
        x, y = contain
        assert x in fs and y in fs and x != y and (x & y) in (x, y), "bad containment witness"
    if sun is not None:
        assert all(s in fs for s in sun) and is_sunflower(*sun), "bad sunflower witness"


# ----------------------------------------------------------------------------------------------------------------
# Port counts, transfer matrices, the Perron root as an exact surd, and exact decimal enclosures
# ----------------------------------------------------------------------------------------------------------------

def port_counts(fam, u, v):
    """N[a][b] = #{S : [u in S] = a, [v in S] = b} (u = in-port, v = out-port, point indices)."""
    N = [[0, 0], [0, 0]]
    for s in fam:
        N[s >> u & 1][s >> v & 1] += 1
    return N


def T_of(N):
    """Transfer matrix of FORMAT.md / [Cla2026]: (avoid, contain) of the out-port, row vector times
    T = [[N10, N11], [N00, N01]] (members containing the old out-port pair with copy members avoiding the in-port)."""
    return ((N[1][0], N[1][1]), (N[0][0], N[0][1]))


def M_irsi(N):
    """The same map in the convention of the [IRSI2026] README, acting on column vectors (contain, avoid):
    M = [[#(v in S, u not in S), #(v in S, u in S)], [#(v not in S, u not in S), #(v not in S, u in S)]]."""
    return ((N[0][1], N[1][1]), (N[0][0], N[1][0]))


def matmul(X, Y):
    return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))


class Surd:
    """The larger root rho = (t + sqrt(D)) / 2 of x^2 - t x + det, D = t^2 - 4 det = (a - d)^2 + 4 b c."""

    def __init__(self, T):
        (a, b), (c, d) = T
        assert min(a, b, c, d) > 0, "the transfer matrix must be positive"
        self.T = T
        self.t = a + d
        self.det = a * d - b * c
        self.D = (a - d) ** 2 + 4 * b * c
        assert self.D == self.t ** 2 - 4 * self.det
        r = isqrt(self.D)
        self.square = r * r == self.D

    def le(self, q):
        """q <= rho, exactly, for a rational q."""
        w = 2 * Fraction(q) - self.t
        return w <= 0 or w * w <= self.D

    def lt(self, q):
        w = 2 * Fraction(q) - self.t
        return w < 0 or w * w < self.D

    def root_le(self, num, den, L):
        """(num/den)^L <= rho, in integers: 2 num^L - t den^L <= sqrt(D) den^L."""
        nl, dl = num ** L, den ** L
        w = 2 * nl - self.t * dl
        return w <= 0 or w * w <= self.D * dl * dl

    def approx(self, prec=60):
        getcontext().prec = prec
        return (Decimal(self.t) + Decimal(self.D).sqrt()) / 2

    def rate_floor(self, L, k):
        """The integer m with m/10^k <= rho^(1/L) < (m+1)/10^k, decided by exact integer comparisons."""
        getcontext().prec = k + 40
        guess = int((self.approx(k + 40).ln() / L).exp() * (Decimal(10) ** k))
        m = guess
        while not self.root_le(m, 10 ** k, L):
            m -= 1
        while self.root_le(m + 1, 10 ** k, L):
            m += 1
        assert self.root_le(m, 10 ** k, L) and not self.root_le(m + 1, 10 ** k, L)
        return m

    def rho_floor(self, k):
        """The integer m with m/10^k <= rho < (m+1)/10^k."""
        m = int(self.approx(k + 30) * (Decimal(10) ** k))
        while not self.le(Fraction(m, 10 ** k)):
            m -= 1
        while self.le(Fraction(m + 1, 10 ** k)):
            m += 1
        return m


def dec_str(m, k):
    s = str(m).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def check_printed(S, L, printed):
    """A printed lower bound b (decimal string) is valid iff b^L <= rho. Returns (valid, is_truncation, detail)."""
    whole, frac = printed.split(".")
    k = len(frac)
    num = int(whole + frac)
    valid = S.root_le(num, 10 ** k, L)
    trunc = valid and not S.root_le(num + 1, 10 ** k, L)
    return valid, trunc, k


# ----------------------------------------------------------------------------------------------------------------
# Reporting, hashes, red controls
# ----------------------------------------------------------------------------------------------------------------

class Report:
    def __init__(self, title):
        self.fails = []
        self.t0 = time.time()
        print(title, flush=True)

    def info(self, key, val):
        print("  %s: %s" % (key, val), flush=True)

    def check(self, desc, ok, detail=""):
        print("  check %s: %s%s" % (desc, "PASS" if ok else "FAIL", (" (" + detail + ")") if detail else ""),
              flush=True)
        if not ok:
            self.fails.append(desc)
        return ok

    def done(self):
        dt = time.time() - self.t0
        print("  run time: %.1f s" % dt)
        print("  RESULT: %s" % ("ALL CHECKS PASS" if not self.fails else "FAILED: " + "; ".join(self.fails)))
        return not self.fails


# hash fragments stated on the 49a page, in the pull requests and in the claimants' READMEs / SHA256SUMS
STATED = {
    "irsi_v10/cert_c49_chain20.json": [("PR #235", "c30a74fed646", "")],
    "irsi_v11/cert_c49_chain20_s17.json": [("49a page", "c12d0dd6", "fc7b"), ("PR #235", "c12d0dd62500", "")],
    "irsi_v12/cert_c49_chain20_v12.json": [("49a page", "7472112c", "4e82e"), ("PR #235", "7472112ced9d", "")],
    "syl_F/F.txt": [("Syl2026 SHA256SUMS", "2ce1fd40924e062e78fc3ee04a7edd21c1fe9a5dcee94408e9426f5c61b8699c", "")],
    "syl_H/H.txt": [("Syl2026 SHA256SUMS", "ead114dde62ca9472239aefddd69dbea00494d113aca9c93ad4b9e71ef036d83", "")],
    "syl_B/B.txt": [("Syl2026 SHA256SUMS", "340cef5e71124ed624425abf81e9621839406cbf44b3fd242594d8edcb7d37af", "")],
    "syl_Bprime/Bprime.txt": [
        ("Syl2026 README and SHA256SUMS", "2a72385704b8203269cda4f4b0dfdabc74d5adb8f4a2bc1dac1f570cb96b5cec", ""),
        ("PR #237 thread", "2a723857", "5cec")],
    "kos/two_pendant_extra_1.5678599.txt": [
        ("49a page", "10ee18f9", "4838"),
        ("Kos2026 README and SHA256SUMS", "10ee18f93204468f3611d2c909a9eee443f3b4769cf7fddd964209edfebd4838", "")],
    "kos/two_pendant_extra_1.5678599.txt.xz": [
        ("Kos2026 SHA256SUMS", "1304bc12d9c51be5c72233c56485746218b10dd3340a8dc9b1c1f351e2499af0", ""),
        ("PR #260 thread", "1304bc12", "9af0")],
}


def hash_checks(rep, root, rel):
    h = sha256(os.path.join(root, rel))
    rep.info("sha256 " + rel, h)
    stated = STATED.get(rel)
    if not stated:
        rep.info("stated hash", "none stated for this file (pinned by commit only)")
    for where, pre, suf in stated or []:
        rep.check("sha256 matches %s (%s%s)" % (where, pre, ("..." + suf) if suf else ("" if len(pre) == 64 else "...")),
                  h.startswith(pre) and h.endswith(suf))
    return h


def exhaustive(rep, label, fam, methods, blocks=None, workers=1):
    """Run the named exhaustive methods; all must find no containment and no sunflower."""
    fs = set(fam)
    rep.check("%s: members distinct" % label, len(fs) == len(fam), "%d members" % len(fam))
    for meth in methods:
        t = time.time()
        if meth == "brute":
            c, s = brute_check(fam)
            extra = "%d pairs, %d triples" % (comb(len(fam), 2), comb(len(fam), 3))
        elif meth == "pivot":
            c, s = pivot_check(fam, workers)
            extra = "%d pivots" % len(fam)
        elif meth == "pair":
            c, s = pair_check(fam, workers)
            extra = "%d pairs" % comb(len(fam), 2)
        else:
            c, s, info = box_check(fam, blocks[meth])
            extra = "blocks %s, %d boxes, %d box triples" % (
                "+".join(str(popcount(b)) for b in blocks[meth]), info["boxes"], info["box triples"])
        verify_witnesses(fam, c, s)
        rep.check("%s: Sperner and 3-sunflower-free, %s method" % (label, meth), c is None and s is None,
                  "%s; %.1f s%s" % (extra, time.time() - t,
                                    "" if c is None and s is None else "; witness %s %s" % (c, s)))


def plant_sunflower(fam, rng, tries=4000):
    """A new set C = (A & B) | R with R a nonempty subset outside A | B: {A, B, C} is a 3-sunflower. Prefer a C that
    is comparable with no member, so that the family stays Sperner and only the sunflower test can fire."""
    width = max(fam).bit_length()
    fs = set(fam)
    fallback = None
    for _ in range(tries):
        a, b = rng.sample(fam, 2)
        out = ((1 << width) - 1) & ~(a | b)
        if not out:
            continue
        r = 0
        while not r:
            r = out & rng.getrandbits(width)
        c = (a & b) | r
        if c in fs:
            continue
        fallback = fallback or (c, a, b)
        if not any((c & s) in (c, s) for s in fam):
            return c, a, b, True
    return fallback + (False,)


def plant_containment(fam, rng):
    width = max(fam).bit_length()
    fs = set(fam)
    while True:
        s = rng.choice(fam)
        outside = [p for p in range(width) if not s >> p & 1]
        if outside:
            t = s | (1 << rng.choice(outside))
            if t not in fs:
                return t, s


def controls(rep, label, fam, methods, blocks=None, workers=1, seed=49):
    """Red controls on the exhaustive methods: a planted 3-sunflower and a planted containment must be reported."""
    rng = random.Random(seed)
    c, a, b, clean = plant_sunflower(fam, rng)
    t, s = plant_containment(fam, rng)
    for meth in methods:
        for kind, extra in (("planted 3-sunflower", c), ("planted containment", t)):
            fam2 = fam + [extra]
            if meth == "brute":
                cc, ss = brute_check(fam2)
            elif meth == "pivot":
                cc, ss = pivot_check(fam2, workers)
            elif meth == "pair":
                cc, ss = pair_check(fam2, workers)
            else:
                cc, ss, _ = box_check(fam2, blocks[meth])
            verify_witnesses(fam2, cc, ss)
            if kind == "planted 3-sunflower":
                ok = ss is not None and extra in ss and (cc is None) == clean
                detail = "witness verified; the added set is comparable with no member, Sperner still reported" \
                    if clean else "witness verified"
            else:
                ok = cc is not None and extra in cc
                detail = "witness verified"
            rep.check("control %s: %s caught by the %s method" % (label, kind, meth), ok, detail)


def bound_checks(rep, S, L, printed, digits_claims=(), rho_claims=(), points=None):
    """Exact rate enclosure and comparison with the printed lower bound and the stated digits."""
    m10 = S.rho_floor(10)
    rep.info("rho", "(%d + sqrt(%d))/2, in [%s, %s)" % (S.t, S.D, dec_str(m10, 10), dec_str(m10 + 1, 10)))
    rep.check("rho is irrational (D not a perfect square)", not S.square)
    for approx in rho_claims:  # "approximately" values: must be the rounding or the truncation at their digits
        k = len(approx.split(".")[1])
        mk = S.rho_floor(k)
        num = int(approx.replace(".", ""))
        rep.check("stated rho ~ %s agrees" % approx, num in (mk, mk + 1),
                  "rho in [%s, %s)" % (dec_str(mk, k), dec_str(mk + 1, k)))
    m = S.rate_floor(L, 25)
    rep.info("bound rho^(1/%d)" % L, "in [%s, %s)" % (dec_str(m, 25), dec_str(m + 1, 25)))
    valid, trunc, k = check_printed(S, L, printed)
    rep.check("printed lower bound %s <= rho^(1/%d) (exact: %s^%d <= rho)" % (printed, L, printed, L), valid,
              "it is the truncation at %d decimals" % k if trunc else "NOT the truncation at %d decimals" % k)
    up = dec_str(int(printed.replace(".", "")) + 1, k)
    rep.check("control: %s (printed bound + 1e-%d) is refuted" % (up, k), not check_printed(S, L, up)[0])
    for claim, where in digits_claims:  # "1.5621220560..." style: the true value starts with these digits
        kk = len(claim.split(".")[1])
        mk = S.rate_floor(L, kk)
        if where == "page":
            rep.check("49a page's value %s... agrees (truncation)" % claim, dec_str(mk, kk) == claim,
                      "rho^(1/%d) in [%s, %s)" % (L, dec_str(mk, kk), dec_str(mk + 1, kk)))
            continue
        ok = dec_str(mk, kk) == claim or dec_str(mk + 1, kk) == claim
        how = "truncated" if dec_str(mk, kk) == claim else "rounded up, not truncated" if ok else "DISAGREES"
        rep.check("value %s... stated in %s agrees (%s)" % (claim, where, how), ok,
                  "rho^(1/%d) in [%s, %s)" % (L, dec_str(mk, kk), dec_str(mk + 1, kk)))
    m10r = S.rate_floor(L, 10)
    over = dec_str(m10r + 1, 10)
    rep.check("control: %s (10-decimal truncation + 1e-10) is refuted" % over, not check_printed(S, L, over)[0])
    return m


# ----------------------------------------------------------------------------------------------------------------
# Rows
# ----------------------------------------------------------------------------------------------------------------

def matrix_checks(rep, N, rows, label="", convention="T"):
    T = T_of(N) if convention == "T" else M_irsi(N)
    rep.info("port counts N00, N01, N10, N11%s" % label, "%d, %d, %d, %d" % (N[0][0], N[0][1], N[1][0], N[1][1]))
    if rows is not None:
        rep.check("transfer matrix rows %s, %s as on the 49a page%s" % (rows[0], rows[1], label), T == tuple(rows),
                  "computed %s, %s" % T)
    return T


def perturbation_controls(rep, N, rows, S, L, printed, convention="T"):
    """A changed count must change the matrix, the Perron root and the rate; a 0.1% smaller N11 must refute the
    printed bound."""
    conv = T_of if convention == "T" else M_irsi
    N2 = [r[:] for r in N]
    N2[1][1] += 1
    S2 = Surd(conv(N2))
    m1, m2 = S.rate_floor(L, 25), S2.rate_floor(L, 25)
    rep.check("control: N11 + 1 no longer gives the page's matrix and moves the exact rate",
              conv(N2) != tuple(rows) and m1 != m2 and (S2.t, S2.D) != (S.t, S.D),
              "rate %s -> %s" % (dec_str(m1, 25), dec_str(m2, 25)))
    N3 = [r[:] for r in N]
    N3[1][1] -= max(1, N3[1][1] // 1000)
    S3 = Surd(conv(N3))
    rep.check("control: with N11 lowered by 0.1%% (to %d) the printed bound %s is refuted" % (N3[1][1], printed),
              not check_printed(S3, L, printed)[0], "rate would be %s" % dec_str(S3.rate_floor(L, 10), 10))


IRSI_ROWS = {
    "irsi_v10": dict(rel="irsi_v10/cert_c49_chain20.json", tag="v1.0", size=6127, layers=None,
                     rows=((1884, 1360), (1081, 1802)), poly=None, closed=None, printed="1.5617847",
                     rho_approx=["3056.1945"], digits=[("1.5617847089", "PR #235"), ("1.56178470895053", "[IRSI2026] README")]),
    "irsi_v11": dict(rel="irsi_v11/cert_c49_chain20_s17.json", tag="v1.1", size=6150,
                     layers={7: 2, 10: 360, 11: 1108, 12: 3600, 13: 1080},
                     rows=((1802, 1360), (1080, 1908)), poly=(3710, 1969416),
                     closed=("1855 + sqrt(1471609)", 3710, 4 * 1471609), printed="1.5621220",
                     rho_approx=["3068.0989"], digits=[("1.5621220560", "page"), ("1.5621220560", "PR #235"), ("1.56212205600681", "[IRSI2026] README")]),
    "irsi_v12": dict(rel="irsi_v12/cert_c49_chain20_v12.json", tag="v1.2", size=6341,
                     layers={7: 2, 10: 551, 11: 1108, 12: 3600, 13: 1080},
                     rows=((1993, 1360), (1080, 1908)), poly=(3901, 2333844),
                     closed=("(3901 + sqrt(5882425))/2", 3901, 5882425), printed="1.5647730",
                     rho_approx=["3163.1856"], digits=[("1.5647730960", "page"), ("1.5647730960", "PR #235"), ("1.56477309608283", "[IRSI2026] README"), ("1.564773096082", "PR #200")]),
}


def row_irsi(root, key, do_controls, workers):
    spec = IRSI_ROWS[key]
    rep = Report("ROW %s: [IRSI2026] tag %s, mu_3 >= %s*" % (key, spec["tag"], spec["printed"]))
    hash_checks(rep, root, spec["rel"])
    d, bases, word = read_irsi(os.path.join(root, spec["rel"]))
    rep.check("one base family B on 20 points, chained with ports u = 18, v = 19",
              list(bases) == ["B"] and bases["B"][0] == 20 and word == [("B", 18, 19)], "word %s" % word)
    n, fam = bases["B"]
    rep.check("family size %d" % spec["size"], len(fam) == spec["size"], "%d sets" % len(fam))
    layers = dict(sorted(Counter(popcount(s) for s in fam).items()))
    rep.info("layer sizes (set size: count)", layers)
    if spec["layers"]:
        rep.check("layer sizes as on the 49a page", layers == spec["layers"])
    exhaustive(rep, "B", fam, ["pivot", "pair"], workers=workers)
    N = port_counts(fam, 18, 19)
    M = matrix_checks(rep, N, spec["rows"], " (u = 18, v = 19; [IRSI2026] README convention)", convention="M")
    T = T_of(N)
    rep.info("same map in the FORMAT.md / [Cla2026] convention", "T = %s, %s (same trace and determinant)" % T)
    if "M" in d:
        rep.check("matrix recorded inside the certificate equals the recomputed one", tuple(map(tuple, d["M"])) == M)
    S = Surd(M)
    rep.info("characteristic polynomial", "x^2 - %d x + %d" % (S.t, S.det))
    if spec["poly"]:
        rep.check("characteristic polynomial as on the 49a page", (S.t, S.det) == spec["poly"])
    if spec["closed"]:
        text, t, D = spec["closed"]
        rep.check("Perron root equals %s as on the 49a page" % text, (S.t, S.D) == (t, D))
    rep.check("same Perron root in both conventions", (Surd(T).t, Surd(T).D) == (S.t, S.D))
    bound_checks(rep, S, 18, spec["printed"], spec["digits"], spec["rho_approx"])
    rep.info("points per copy", "%d - 2 = 18, so the bound is rho^(1/18)" % n)
    if do_controls:
        controls(rep, "B", fam, ["pivot", "pair"], workers=workers)
        perturbation_controls(rep, N, spec["rows"], S, 18, spec["printed"], convention="M")
        dropped = port_counts(fam[1:], 18, 19)
        rep.check("control: dropping one member changes the port counts", dropped != N)
    return rep.done()


def row_degkm(root, do_controls, workers):
    rep = Report("ROW irsi_degkm: [IRSI2026] DEGKM control (not an asterisked row): DEGKM's families F, L chained")
    rel = "irsi_degkm/cert_degkm_control.json"
    hash_checks(rep, root, rel)
    d, bases, word = read_irsi(os.path.join(root, rel))
    P = ((1, 0), (0, 1))
    pts = 0
    for name, u, v in word:
        n, fam = bases[name]
        exhaustive(rep, name, fam, ["pivot", "pair"], workers=workers)
        N = port_counts(fam, u, v)
        rep.info("%s: %d sets on %d points, ports %d -> %d, N00, N01, N10, N11" % (name, len(fam), n, u, v),
                 "%d, %d, %d, %d" % (N[0][0], N[0][1], N[1][0], N[1][1]))
        P = matmul(P, T_of(N))
        pts += n - 2
    rep.info("period matrix (FORMAT.md convention)", "%s, %s" % P)
    S = Surd(P)
    bound_checks(rep, S, pts, "1.5529675")
    return rep.done()


def isomorphic(A, B, n):
    """Is there a permutation of range(n) mapping the family A (bitmasks) onto the family B? Backtracking on points,
    pruned by the multiset of traces of the members on the points mapped so far."""
    if len(A) != len(B) or sorted(map(popcount, A)) != sorted(map(popcount, B)):
        return None
    def traces(fam, pts):
        return Counter(tuple(s >> p & 1 for p in pts) for s in fam)
    def extend(srcs, dsts):
        if len(srcs) == n:
            return dict(zip(srcs, dsts))
        x = len(srcs)
        for y in range(n):
            if y in dsts:
                continue
            if traces(A, srcs + [x]) == traces(B, dsts + [y]):
                r = extend(srcs + [x], dsts + [y])
                if r:
                    return r
        return None
    phi = extend([], [])
    if phi:
        assert {mask(phi[p] for p in bits(s)) for s in A} == set(B)
    return phi


def F_masks():
    return to_masks(family_F(), K5_EDGES)


def H_masks():
    return to_masks(family_H(), list(range(11)))


def row_cla_F(root, do_controls, workers):
    rep = Report("ROW cla_F: [Cla2026] 75 sets on the edges of K5, mu_3 >= 1.5664286*")
    for rel in ("cla_k5/k5_family_masks.txt",):
        hash_checks(rep, root, rel)
    F = F_masks()
    filed = read_mask_list(os.path.join(root, "cla_k5/k5_family_masks.txt"))
    rep.check("family rebuilt from the page's words equals k5_family_masks.txt (edges in lexicographic order)",
              set(filed) == set(F) and len(filed) == 75)
    rep.info("set sizes", dict(sorted(Counter(map(popcount, F)).items())))
    perms = list(itertools.permutations(range(5)))
    idx = {e: i for i, e in enumerate(K5_EDGES)}
    def act(pi, s):
        return mask(idx[tuple(sorted((pi[a], pi[b])))] for a, b in (K5_EDGES[p] for p in bits(s)))
    rep.check("S5-invariant (all 120 vertex permutations)", all({act(pi, s) for s in F} == set(F) for pi in perms))
    exhaustive(rep, "F", F, ["brute", "pivot", "pair"])
    deg = [sum(s >> p & 1 for s in F) for p in range(10)]
    co = [sum(1 for s in F if s >> p & 1 and s >> q & 1) for p, q in itertools.combinations(range(10), 2)]
    rep.check("every edge lies in 44 members, every pair of edges in 24", set(deg) == {44} and set(co) == {24})
    Ts = {T_of(port_counts(F, u, v)) for u in range(10) for v in range(10) if u != v}
    rep.check("for every ordered pair of ports the matrix has rows (20, 24), (11, 20)", Ts == {((20, 24), (11, 20))})
    S = Surd(((20, 24), (11, 20)))
    rep.check("Perron root equals 20 + 2 sqrt(66)", S.t == 40 and S.D == 16 * 66)
    bound_checks(rep, S, 8, "1.5664286", [("1.5664286", "page"), ("1.5664286622", "PR #235 and [Cla2026] README")])
    if do_controls:
        controls(rep, "F", F, ["brute", "pivot", "pair"])
        perturbation_controls(rep, port_counts(F, 0, 1), ((20, 24), (11, 20)), S, 8, "1.5664286")
    return rep.done()


def row_cla_FH(root, do_controls, workers):
    rep = Report("ROW cla_FH: [Cla2026] F alternated with 121 sets on Z_11, mu_3 >= 1.5666969*")
    hash_checks(rep, root, "cla_z11/z11_family_masks.txt")
    F, H = F_masks(), H_masks()
    filed = read_mask_list(os.path.join(root, "cla_z11/z11_family_masks.txt"))
    rep.check("family rebuilt from the page's words equals z11_family_masks.txt", set(filed) == set(H) and
              len(filed) == 121)
    rep.info("set sizes", dict(sorted(Counter(map(popcount, H)).items())))
    exhaustive(rep, "H", H, ["brute", "pivot", "pair"])
    deg = [sum(s >> p & 1 for s in H) for p in range(11)]
    co = [sum(1 for s in H if s >> p & 1 and s >> q & 1) for p, q in itertools.combinations(range(11), 2)]
    rep.check("every point lies in 75 members, every pair of points in 44", set(deg) == {75} and set(co) == {44})
    Ts = {T_of(port_counts(H, u, v)) for u in range(11) for v in range(11) if u != v}
    rep.check("for every ordered pair of ports the matrix has rows (31, 44), (15, 31)", Ts == {((31, 44), (15, 31))})
    links_ok = True
    for p in range(11):
        link = [s & ~(1 << p) for s in H if s >> p & 1]
        others = [q for q in range(11) if q != p]
        compact = [mask(others.index(q) for q in bits(s)) for s in link]
        links_ok &= isomorphic(compact, F, 10) is not None
    rep.check("the link of every point is a copy of F (explicit isomorphism found for each of the 11 points)",
              links_ok)
    TF, TH = ((20, 24), (11, 20)), ((31, 44), (15, 31))
    P = matmul(TF, TH)
    rep.check("product of the two matrices has rows (980, 1624), (641, 1104)", P == ((980, 1624), (641, 1104)),
              "trace %d, det %d" % (P[0][0] + P[1][1], P[0][0] * P[1][1] - P[0][1] * P[1][0]))
    S = Surd(P)
    S2 = Surd(matmul(TH, TF))
    rep.check("Perron root equals 1042 + 2 sqrt(261207) (and is the same for the other order)",
              S.t == 2084 and S.D == 16 * 261207 and (S2.t, S2.D) == (S.t, S.D))
    rep.info("points per period", "(10 - 2) + (11 - 2) = 17")
    bound_checks(rep, S, 17, "1.5666969", [("1.5666969", "page"), ("1.5666969787", "PR #235 and [Cla2026] README")])
    if do_controls:
        controls(rep, "H", H, ["brute", "pivot", "pair"])
        N = port_counts(H, 0, 1)
        N2 = [r[:] for r in N]
        N2[1][1] += 1
        S3 = Surd(matmul(TF, T_of(N2)))
        rep.check("control: N11(H) + 1 moves the period's Perron root", (S3.t, S3.D) != (S.t, S.D))
    return rep.done()


# Sylla's labels: K5 edges in lexicographic order are F points 0..9; in B the F points other than (0,3) become 0..8
# and the points 1..10 of H become 9..18; the chain ports are 0 = (0,1) and 1 = (0,2).
SYL_B_ORDER = [("F", e) for e in K5_EDGES if e != (0, 3)] + [("H", x) for x in range(1, 11)]

SYL_ADDED = """02 04 13 14 23 | 01 02 04 12 13 23 24 | 01 02 12 13 14 23 24 | 01 04 12 13 14 23 24 | 01 02 04 12 13 14 34 |
01 04 12 23 34 | 01 02 12 13 14 23 34 | 01 02 13 24 34 | 01 04 12 13 14 24 34 | 02 04 12 13 14 24 34 |
02 04 12 13 23 24 34 | 01 02 04 14 23 24 34 | 01 02 12 14 23 24 34 | 02 04 12 14 23 24 34 | 01 04 13 14 23 24 34"""


def edge_sets(text):
    return [frozenset((int(w[0]), int(w[1])) for w in part.split()) for part in text.replace("\n", " ").split("|")]


def build_B():
    return glue(relabel(family_F(), "F"), ("F", (0, 3)), relabel(family_H(), "H"), ("H", 0))


def build_G():
    """[Syl2026] gadget: F with copies of H glued (at 0) at the edges (1,2) and (1,4); 28 points."""
    G1 = glue(relabel(family_F(), "F"), ("F", (1, 2)), relabel(family_H(), "H1"), ("H1", 0))
    G = glue(G1, ("F", (1, 4)), relabel(family_H(), "H2"), ("H2", 0))
    order = [("F", e) for e in K5_EDGES if e not in ((1, 2), (1, 4))] + \
            [("H1", x) for x in range(1, 11)] + [("H2", x) for x in range(1, 11)]
    return G, order


def row_syl_unit(root, do_controls, workers):
    rep = Report("ROW syl_unit: [Syl2026] @914de15, the unit chain, mu_3 >= 1.5675517*")
    for rel in ("syl_F/F.txt", "syl_H/H.txt", "syl_B/B.txt"):
        hash_checks(rep, root, rel)
    F, H = F_masks(), H_masks()
    rep.check("lemma hypotheses: F and H are Sperner and 3-sunflower-free (brute force over all pairs and triples)",
              brute_check(F) == (None, None) and brute_check(H) == (None, None))
    rep.check("data/F.txt and data/H.txt equal F and H rebuilt from the page's words",
              set(read_point_lines(os.path.join(root, "syl_F/F.txt"))) == set(F) and
              set(read_point_lines(os.path.join(root, "syl_H/H.txt"))) == set(H))
    B = to_masks(build_B(), SYL_B_ORDER)
    filed = read_point_lines(os.path.join(root, "syl_B/B.txt"))
    rep.check("data/B.txt equals the gluing of F at (0,3) and H at 0 (4349 sets on 19 points)",
              set(filed) == set(B) and len(filed) == len(B) == 4349)
    blocks = {"box": [mask(range(9)), mask(range(9, 19))]}
    exhaustive(rep, "B", filed, ["box", "pivot", "pair"], blocks, workers)
    N = port_counts(filed, 0, 1)
    T = matrix_checks(rep, N, ((1123, 1481), (622, 1123)), " of B, ports (0,1) -> (0,2)")
    SB = Surd(T)
    rep.check("first step: Perron root 1123 + sqrt(921182)", SB.t == 2246 and SB.D == 4 * 921182)
    mB = SB.rate_floor(17, 10)
    rep.check("first step: (1123 + sqrt(921182))^(1/17) = 1.5675245...", dec_str(mB, 7) == "1.5675245" or
              dec_str(mB, 10).startswith("1.5675245"), "in [%s, %s)" % (dec_str(mB, 10), dec_str(mB + 1, 10)))
    # the gadget
    Gl, order = build_G()
    G = to_masks(Gl, order)
    root_pt = order.index(("F", (0, 1)))
    avoid = sum(1 for s in G if not s >> root_pt & 1)
    rep.check("gadget G: 250659 sets on 28 points, 98308 avoiding and 152351 containing the root (0,1)",
              len(G) == len(set(G)) == 250659 and len(order) == 28 and (avoid, len(G) - avoid) == (98308, 152351),
              "%d sets, %d / %d" % (len(G), avoid, len(G) - avoid))
    gblocks = {"box": [mask(range(8)), mask(range(8, 18)), mask(range(18, 28))],
               "box2": [mask(range(4)), mask(range(4, 8)), mask(range(8, 18)), mask(range(18, 28))]}
    exhaustive(rep, "G", G, ["box", "box2"], gblocks)
    rep.info("G, second path", "G is two applications of the gluing lemma to F and H, which are valid")
    spot = sorted(random.Random(60).sample(range(len(G)), 400))
    t = time.time()
    c, s = additions_scan(G, spot, workers)
    verify_witnesses(G, c, s)
    rep.check("G, spot check by the generic pivot method: no containment and no 3-sunflower through 400 random "
              "members (each scanned against all others; not exhaustive)", c is None and s is None,
              "%.1f s" % (time.time() - t))
    # the unit: F core, ports (0,1) -> (0,2), copies of G glued (at their roots) at (0,3) and (0,4)
    w = {1: 98308, 0: 152351}  # glued point in C: the G-part avoids its root; otherwise it contains it
    Fl = family_F()
    U = [[0, 0], [0, 0]]
    for C in Fl:
        U[(0, 1) in C][(0, 2) in C] += w[(0, 3) in C] * w[(0, 4) in C]
    TU = matrix_checks(rep, U, None, " of the 62-point unit (gluing count)")
    rep.check("unit port counts as on the 49a page", (U[0][0], U[0][1], U[1][0], U[1][1]) ==
              (148811965456, 279351870092, 279351870092, 378843386627))
    # the same count formula against an explicit build, with F (root (0,1)) in place of G
    small = glue(glue(relabel(Fl, "c"), ("c", (0, 3)), relabel(Fl, "g1"), ("g1", (0, 1))),
                 ("c", (0, 4)), relabel(Fl, "g2"), ("g2", (0, 1)))
    sorder = sorted(ground(small))
    sm = to_masks(small, sorder)
    wF = {1: sum(1 for C in Fl if (0, 1) not in C), 0: sum(1 for C in Fl if (0, 1) in C)}
    pred = [[0, 0], [0, 0]]
    for C in Fl:
        pred[(0, 1) in C][(0, 2) in C] += wF[(0, 3) in C] * wF[(0, 4) in C]
    got = port_counts(sm, sorder.index(("c", (0, 1))), sorder.index(("c", (0, 2))))
    sblocks = {"box": [mask(i for i, x in enumerate(sorder) if x[0] == tag) for tag in ("c", "g1", "g2")]}
    rep.check("count formula against an explicit unit with F as the gadget (%d sets on %d points)" %
              (len(sm), len(sorder)), got == pred and len(set(sm)) == len(sm))
    exhaustive(rep, "explicit F-gadget unit", sm, ["box"], sblocks)
    S = Surd(TU)
    rep.check("Perron root equals 279351870092 + sqrt(378843386627 * 148811965456)",
              S.t == 2 * 279351870092 and S.D == 4 * 378843386627 * 148811965456)
    rep.info("points", "unit 10 + 2 (28 - 2) = 62, 60 added per unit")
    bound_checks(rep, S, 60, "1.5675517", [("1.5675517", "page"), ("1.5675517557", "PR #237"), ("1.5675517557582", "[Syl2026] README")])
    if do_controls:
        controls(rep, "B", filed, ["box", "pivot", "pair"], blocks, workers)
        controls(rep, "G", G, ["box"], gblocks)
        perturbation_controls(rep, U, TU, S, 60, "1.5675517")
    return rep.done()


def row_syl_Bprime(root, do_controls, workers):
    rep = Report("ROW syl_Bprime: [Syl2026] @8e21d02, mu_3 >= 1.5678015*")
    hash_checks(rep, root, "syl_Bprime/Bprime.txt")
    Bp = read_point_lines(os.path.join(root, "syl_Bprime/Bprime.txt"))
    B = to_masks(build_B(), SYL_B_ORDER)
    added = [s for s in Bp if s not in set(B)]
    listed = to_masks(edge_sets(SYL_ADDED), [e for _, e in SYL_B_ORDER[:9]])
    rep.check("B' = B (the 4349-set gluing) plus 15 sets", set(B) <= set(Bp) and len(Bp) == 4364 and len(added) == 15)
    rep.check("the 15 added sets lie in the F-part and are the 15 listed in the [Syl2026] README",
              all(s < (1 << 9) for s in added) and set(added) == set(listed) and Bp[-15:] == added)
    blocks = {"box": [mask(range(9)), mask(range(9, 19))]}
    exhaustive(rep, "B'", Bp, ["box", "pivot", "pair"], blocks, workers)
    N = port_counts(Bp, 0, 1)
    T = matrix_checks(rep, N, ((1127, 1488), (622, 1127)), " of B', ports (0,1) -> (0,2)")
    rep.check("start vector (u1, w1) = (1749, 2615)", (N[0][0] + N[1][0], N[0][1] + N[1][1]) == (1749, 2615))
    S = Surd(T)
    rep.check("Perron root equals 1127 + sqrt(925536)", S.t == 2254 and S.D == 4 * 925536)
    rep.info("points per copy", "19 - 2 = 17 (chain families on 17k + 2 points)")
    bound_checks(rep, S, 17, "1.5678015", [("1.5678015", "page"), ("1.5678015201", "PR #237"), ("1.5678015201687", "[Syl2026] README")])
    if do_controls:
        controls(rep, "B'", Bp, ["box", "pivot", "pair"], blocks, workers)
        perturbation_controls(rep, N, T, S, 17, "1.5678015")
        rep.check("control: dropping one added set changes the port counts", port_counts(Bp[:-1], 0, 1) != N)
    return rep.done()


KOS_EDGES = [(0, 1), (0, 2), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]  # points 0..7; 8..17 H1, 18..27 H2
KOS_A1 = """01 02 04 12 13 14 23 | 01 04 13 23 24 | 01 02 12 13 14 23 24 | 02 04 12 13 14 23 24 | 02 04 12 13 34 |
01 02 14 23 34 | 01 04 12 13 14 23 34 | 01 02 04 13 14 24 34 | 01 02 12 13 14 24 34 | 01 04 12 13 14 24 34 |
01 02 04 12 23 24 34 | 01 02 12 13 23 24 34 | 01 04 12 14 23 24 34 | 02 04 12 14 23 24 34 | 02 04 13 14 23 24 34"""
KOS_A2 = """02 03 13 14 24 | 01 02 03 12 14 23 24 | 01 02 12 13 14 23 24 | 01 03 12 13 14 23 24 | 01 02 03 12 13 14 34 |
01 02 14 23 34 | 01 03 12 13 14 23 34 | 02 03 12 13 14 23 34 | 01 03 12 24 34 | 01 02 12 13 14 24 34 |
01 02 03 13 23 24 34 | 01 02 12 13 23 24 34 | 02 03 12 13 23 24 34 | 02 03 12 14 23 24 34 | 01 03 13 14 23 24 34"""
KOS_W = """1 2 3 5 6 8 | 2 3 4 6 7 9 | 1 2 3 4 5 7 8 9 | 1 2 4 5 6 7 8 9 | 1 3 4 5 6 7 8 9 | 1 2 3 4 5 6 7 10 |
3 4 5 7 8 10 | 1 2 3 4 6 7 8 10 | 1 2 4 5 6 7 8 10 | 1 2 3 4 5 6 9 10 | 1 2 3 4 5 7 9 10 | 1 5 6 7 9 10 |
1 2 4 8 9 10 | 1 3 4 5 6 8 9 10 | 2 3 4 5 6 8 9 10 | 1 2 3 5 7 8 9 10 | 1 2 3 6 7 8 9 10 | 1 3 4 6 7 8 9 10 |
2 3 5 6 7 8 9 10 | 2 4 5 6 7 8 9 10"""


def build_kos():
    """[Kos2026] block rebuilt from the 49a page and the lists in its PROOF.md, as (part name, mask) pairs."""
    e1, e2 = (0, 3), (0, 4)
    eidx = {e: i for i, e in enumerate(KOS_EDGES)}
    H = family_H()
    h_avoid = [mask(x - 1 for x in D) for D in H if 0 not in D]  # 46 parts, if the glued edge is in C
    h_contain = [mask(x - 1 for x in D if x) for D in H if 0 in D]  # 75 parts, otherwise
    def parts(present):
        return h_avoid if present else h_contain
    def core(T):
        return mask(eidx[e] for e in T if e not in (e1, e2))
    out = []
    for C in family_F():
        for q1 in parts(e1 in C):
            for q2 in parts(e2 in C):
                out.append(("product", core(C) | q1 << 8 | q2 << 18))
    for T in edge_sets(KOS_A1):
        assert e1 not in T
        out += [("additions at H1", core(T) | q2 << 18) for q2 in parts(e2 in T)]
    for T in edge_sets(KOS_A2):
        assert e2 not in T
        out += [("additions at H2", core(T) | q1 << 8) for q1 in parts(e1 in T)]
    for part in KOS_W.replace("\n", " ").split("|"):
        out.append(("inside H1", mask(int(x) - 1 for x in part.split()) << 8))
    out.append(("extra set", core([(1, 2), (1, 4), (2, 3)]) | ((1 << 10) - 1) << 18))
    return out


def additions_scan(fam, pivots, workers):
    """Every pair and triple of members that involves one of the given members: the pivot method with those
    pivots, each scanning all other members."""
    import multiprocessing as mp
    _SHARED["fam"] = fam
    _SHARED["pivots"] = pivots
    ctx = mp.get_context("fork")
    chunks = [pivots[r::workers * 4] for r in range(workers * 4)]
    contain = sun = None
    with ctx.Pool(workers) as pool:
        for c, s in pool.imap_unordered(_scan_worker, chunks):
            contain = contain or c
            sun = sun or s
    _SHARED.clear()
    return contain, sun


def _scan_worker(idx):
    fam = _SHARED["fam"]
    contain = sun = None
    for i in idx:
        c, s, _ = pivot_range(fam, i, i + 1, 1, everyone=True)
        contain = contain or c
        sun = sun or s
        if sun:
            break
    return contain, sun


def row_kos(root, do_controls, workers):
    rep = Report("ROW kos: [Kos2026] @c4ee817, mu_3 >= 1.5678599*")
    hash_checks(rep, root, "kos/two_pendant_extra_1.5678599.txt.xz")
    hash_checks(rep, root, "kos/two_pendant_extra_1.5678599.txt")
    blocks_xz = read_sunflower3(os.path.join(root, "kos/two_pendant_extra_1.5678599.txt.xz"))
    blocks_txt = read_sunflower3(os.path.join(root, "kos/two_pendant_extra_1.5678599.txt"))
    rep.check("the .xz decompresses to the same certificate as the .txt", blocks_xz == blocks_txt)
    rep.check("one block: 28 points, in-port 0 = (0,1), out-port 1 = (0,2), 252292 members",
              len(blocks_txt) == 1 and blocks_txt[0][:3] == (28, 0, 1) and len(blocks_txt[0][3]) == 252292)
    n, p, q, fam = blocks_txt[0]
    rep.info("set sizes", dict(sorted(Counter(map(popcount, fam)).items())))
    kblocks = {"box": [mask(range(8)), mask(range(8, 18)), mask(range(18, 28))],
               "box2": [mask(range(4)), mask(range(4, 8)), mask(range(8, 18)), mask(range(18, 28))]}
    exhaustive(rep, "block", fam, ["box", "box2"], kblocks)
    N = port_counts(fam, p, q)
    T = matrix_checks(rep, N, ((64728, 88867), (33969, 64728)), " (ports 0 -> 1)")
    rep.check("port counts N00 = 33969, N01 = N10 = 64728, N11 = 88867 as on the 49a page",
              (N[0][0], N[0][1], N[1][0], N[1][1]) == (33969, 64728, 64728, 88867))
    # second path: rebuild from the description, gluing lemma for the product part, scan through the additions
    built = build_kos()
    rep.check("lemma hypotheses: F and H (rebuilt from the page's words) are Sperner and 3-sunflower-free (brute)",
              brute_check(F_masks()) == (None, None) and brute_check(H_masks()) == (None, None))
    rep.check("the certificate equals the block rebuilt from the page's description and the lists in PROOF.md",
              set(m for _, m in built) == set(fam) and len(built) == len(set(m for _, m in built)) == len(fam))
    comp = Counter(part for part, _ in built)
    rep.info("composition", dict(comp))
    rep.check("composition 250659 + 806 + 806 + 20 + 1 = 252292",
              (comp["product"], comp["additions at H1"], comp["additions at H2"], comp["inside H1"],
               comp["extra set"]) == (250659, 806, 806, 20, 1))
    by_part = {}
    for part, m in built:
        c = by_part.setdefault(part, [0, 0, 0])
        c[(m >> p & 1) + (m >> q & 1)] += 1
    rep.info("port counts by part (N00, N01 + N10, N11)", by_part)
    rep.check("port counts by part as in PR #260", by_part == {
        "product": [33948, 2 * 64360, 87991], "additions at H1": [0, 2 * 184, 438],
        "additions at H2": [0, 2 * 184, 438], "inside H1": [20, 0, 0], "extra set": [1, 0, 0]})
    index = {m: i for i, m in enumerate(fam)}
    adds = sorted(index[m] for part, m in built if part != "product")
    t = time.time()
    c, s = additions_scan(fam, adds, workers)
    verify_witnesses(fam, c, s)
    rep.check("second path: no containment and no 3-sunflower through any of the %d added sets "
              "(each scanned against all %d members; the 250659 product sets are valid by the gluing lemma)"
              % (len(adds), len(fam)), c is None and s is None, "%.1f s" % (time.time() - t))
    prods = [index[m] for part, m in built if part == "product"]
    spot = sorted(random.Random(26).sample(prods, 400))
    t = time.time()
    c, s = additions_scan(fam, spot, workers)
    verify_witnesses(fam, c, s)
    rep.check("spot check of the gluing-lemma part: no containment and no 3-sunflower through 400 random product "
              "sets (each scanned against all members; not exhaustive)", c is None and s is None,
              "%.1f s" % (time.time() - t))
    S = Surd(T)
    rep.check("Perron root equals 64728 + sqrt(3018723123), 3018723123 = 88867 * 33969",
              S.t == 2 * 64728 and S.D == 4 * 3018723123 and 3018723123 == 88867 * 33969)
    rep.info("points per copy", "28 - 2 = 26 (chain families on 26k + 2 points)")
    bound_checks(rep, S, 26, "1.5678599", [("1.5678599", "page"), ("1.5678599024", "PR #260"),
                                           ("1.5678599024245", "certificate header"),
                                           ("1.5678599024245898858529", "PR #260 (their independent check)")])
    if do_controls:
        controls(rep, "block", fam, ["box", "box2"], kblocks)
        perturbation_controls(rep, N, T, S, 26, "1.5678599")
        rep.check("control: dropping the extra set changes the port counts",
                  port_counts([m for m in fam if m != built[-1][1]], p, q) != N)
    return rep.done()


# ----------------------------------------------------------------------------------------------------------------
# The gluing lemma on instances, and the methods against brute force
# ----------------------------------------------------------------------------------------------------------------

def glue_masks(A, a, B, b):
    """Lemma 4.4 on bitmask families with disjoint supports; a, b are point indices."""
    A1 = [C & ~(1 << a) for C in A if C >> a & 1]
    A0 = [C for C in A if not C >> a & 1]
    B0 = [D for D in B if not D >> b & 1]
    B1 = [D & ~(1 << b) for D in B if D >> b & 1]
    return [P | Q for P in A1 for Q in B0] + [P | Q for P in A0 for Q in B1], len(A1) * len(B0) + len(A0) * len(B1)


def valid_families(n):
    """All Sperner, 3-sunflower-free families of subsets of an n-set (n <= 4), by brute force."""
    subsets = list(range(1 << n))
    out = []
    for code in range(1 << len(subsets)):
        fam = [s for s in subsets if code >> s & 1]
        if all((a & b) not in (a, b) for a, b in itertools.combinations(fam, 2)) and brute_check(fam) == (None, None):
            out.append(fam)
    return out


def chain(blocks_spec):
    """Explicit gluing chain: blocks_spec = [(tag, family as frozensets of labels, in-port, out-port)]."""
    tag, fam, _, out = blocks_spec[0]
    R = relabel(fam, tag)
    cur_out = (tag, out)
    for tag, fam, inp, out in blocks_spec[1:]:
        R = glue(R, cur_out, relabel(fam, tag), (tag, inp))
        cur_out = (tag, out)
    return R


def recursion_sizes(Ns):
    """|R_k| from the transfer-matrix recursion: (u1, w1) = (N00 + N10, N01 + N11), then (u, w) T for each block."""
    N = Ns[0]
    u, w = N[0][0] + N[1][0], N[0][1] + N[1][1]
    sizes = [u + w]
    for N in Ns[1:]:
        (a, b), (c, d) = T_of(N)
        u, w = u * a + w * c, u * b + w * d
        sizes.append(u + w)
    return sizes


def row_lemma(root, do_controls, workers):
    rep = Report("ROW lemma: the gluing lemma ([DEGKM1997] Lemma 4.4) and the chain count, tested on instances")
    t = time.time()
    V3, V4 = valid_families(3), valid_families(4)
    rep.info("valid (Sperner, 3-sunflower-free) families", "%d on 3 points, %d on 4 points" % (len(V3), len(V4)))
    bad = 0
    count = 0
    rng = random.Random(4)
    pairs = [(A, 3, B, 3) for A in V3 for B in V3] + [(A, 4, B, 3) for A in V4 for B in V3] + \
            [(rng.choice(V4), 4, rng.choice(V4), 4) for _ in range(3000)]
    for A, na, B, nb in pairs:
        Bs = [D << na for D in B]
        for a in range(na):
            for b in range(nb):
                G, predicted = glue_masks(A, a, Bs, na + b)
                count += 1
                if len(set(G)) != len(G) or len(G) != predicted or brute_check(G) != (None, None):
                    bad += 1
    rep.check("every gluing of two valid families is valid and has |A_a||B_not b| + |A_not a||B_b| members",
              bad == 0, "%d gluings: all valid families on 3+3 and 4+3 points, 3000 random 4+4 pairs; %.1f s"
              % (count, time.time() - t))
    # a control: the lemma needs Sperner; with a non-Sperner input the gluing can fail
    A = [0b01, 0b11]  # {0} inside {0,1}
    B = [0b1 << 2, 0b10 << 2, 0b100 << 2]
    fails = any(brute_check(glue_masks(A, a, B, 2 + b)[0]) != (None, None) for a in range(2) for b in range(3))
    rep.check("control: with a non-Sperner input the glued family can fail (so the hypothesis is used)", fails)
    # explicit chains of [Cla2026] families, checked exhaustively, against the recursion
    Fl, Hl = family_F(), family_H()
    NF = port_counts(F_masks(), K5_EDGES.index((0, 1)), K5_EDGES.index((0, 2)))
    NH = port_counts(H_masks(), 0, 1)
    specs = {
        "F+F": [("a", Fl, (0, 1), (0, 2)), ("b", Fl, (0, 1), (0, 2))],
        "F+H": [("a", Fl, (0, 1), (0, 2)), ("b", Hl, 0, 1)],
        "F+F+F": [("a", Fl, (0, 1), (0, 2)), ("b", Fl, (0, 1), (0, 2)), ("c", Fl, (0, 1), (0, 2))],
        "F+H+F": [("a", Fl, (0, 1), (0, 2)), ("b", Hl, 0, 1), ("c", Fl, (0, 1), (0, 2))],
    }
    stated = {"F+F": (2728, 18), "F+H": (4349, 19), "F+F+F": (98920, 26), "F+H+F": (155892, 27)}
    for name, spec in specs.items():
        R = chain(spec)
        order = sorted(ground(R))
        Rm = to_masks(R, order)
        Ns = [NF if fam is Fl else NH for _, fam, _, _ in spec]
        pred = recursion_sizes(Ns)[-1]
        rep.check("%s: %d sets on %d points, as the recursion and the [Cla2026] README say" % ((name,) + stated[name]),
                  (len(Rm), len(order)) == stated[name] and pred == len(Rm) == len(set(Rm)))
        blocks ={"box": [mask(i for i, x in enumerate(order) if x[0] == tag) for tag, _, _, _ in spec]}
        methods = ["box", "pivot"] if len(Rm) < 5000 else ["box"]
        exhaustive(rep, name, Rm, methods, blocks, workers)
    rep.check("two-copy chains by the recursion: B+B has 9,087,960 members, B'+B' 9,147,270 (as stated)",
              recursion_sizes([[[622, 1123], [1123, 1481]]] * 2)[-1] == 9087960 and
              recursion_sizes([[[622, 1127], [1127, 1488]]] * 2)[-1] == 9147270)
    return rep.done()


def row_selftest(root, do_controls, workers):
    rep = Report("ROW selftest: the pivot, pair and box methods against brute force on random families")
    rng = random.Random(2026)

    def verdicts(fam, blocks_list):
        out = {}
        out["brute"] = brute_check(fam)
        c, s, _ = pivot_range(fam, 0, len(fam), 1)
        out["pivot"] = (c, s)
        has, hasnot = point_bitsets(fam)
        out["pair"] = pair_range(fam, has, hasnot, 0, len(fam), 1)
        for j, blocks in enumerate(blocks_list):
            c, s, _ = box_check(fam, blocks)
            out["box%d" % j] = (c, s)
        for c, s in out.values():
            verify_witnesses(fam, c, s)
        return {k: (c is None, s is None) for k, (c, s) in out.items()}

    def random_blocks(n):
        assign = [rng.randrange(rng.randint(1, 4)) for _ in range(n)]
        return [mask(p for p in range(n) if assign[p] == j) for j in sorted(set(assign))]

    t = time.time()
    disagree, stats = 0, Counter()
    for trial in range(3000):
        n = rng.randint(4, 11)
        N = rng.randint(3, 45)
        if trial % 3 == 0:  # one layer: Sperner, so only the sunflower test can fire
            k = rng.randint(1, n - 1)
            pool = [s for s in range(1 << n) if popcount(s) == k]
            fam = rng.sample(pool, min(N, len(pool)))
        else:
            fam = list({rng.getrandbits(n) for _ in range(N)})
            rng.shuffle(fam)
        if len(fam) < 3:
            continue
        v = verdicts(fam, [random_blocks(n), random_blocks(n)])
        disagree += len(set(v.values())) > 1
        stats[v["brute"]] += 1
    rep.check("3000 random families on 4-11 points: brute, pivot, pair and box (two random point partitions) "
              "give the same verdicts", disagree == 0,
              "(Sperner, sunflower-free) outcomes %s; %.1f s" % (dict(stats), time.time() - t))
    # glued families (valid by construction), then one planted defect
    V4 = [f for f in valid_families(4) if len(f) >= 3]
    t = time.time()
    disagree, stats = 0, Counter()
    for trial in range(400):
        A, B, C = rng.choice(V4), rng.choice(V4), rng.choice(V4)
        G, _ = glue_masks(A, rng.randrange(4), [D << 4 for D in B], 4 + rng.randrange(4))
        support = sorted(set(p for s in G for p in bits(s)) | set(range(4, 8)))
        G, _ = glue_masks(G, rng.choice(support), [D << 8 for D in C], 8 + rng.randrange(4))
        if len(G) < 4:
            continue
        kind = trial % 4
        if kind == 1:
            G = G + [plant_sunflower(G, rng)[0]]
        elif kind == 2:
            G = G + [plant_containment(G, rng)[0]]
        elif kind == 3:
            x = rng.getrandbits(12)
            if x not in G:
                G = G + [x]
        aligned = [mask(range(0, 4)), mask(range(4, 8)), mask(range(8, 12))]
        v = verdicts(G, [aligned, random_blocks(12)])
        disagree += len(set(v.values())) > 1
        stats[(kind, v["brute"])] += 1
    rep.check("400 glued families on 12 points (as glued, and with a planted sunflower, a planted containment or a "
              "random extra set): all methods agree", disagree == 0,
              "(variant, verdict) outcomes %s; %.1f s" % (dict(sorted(stats.items())), time.time() - t))
    rep.check("the glued families themselves are all valid (the gluing lemma on these instances)",
              all(k[1] == (True, True) for k in stats if k[0] == 0) and any(k[0] == 0 for k in stats))
    # the fast pattern step against the cubic reference
    t = time.time()
    same = True
    for trial in range(1500):
        n = rng.randint(2, 9)
        pool = list(range(1 << n))
        X, Y, Z = (frozenset(rng.sample(pool, rng.randint(1, min(12, len(pool))))) for _ in range(3))
        if rng.random() < 0.3:
            Y = X
        if rng.random() < 0.3:
            Z = X
        same &= ws_patterns(X, Y, Z) == ws_patterns_naive(X, Y, Z)
    rep.check("1500 random side triples: the fast pattern step equals the cubic reference", same,
              "%.1f s" % (time.time() - t))
    # a decomposition that does not rebuild the family must be refused by box_check
    fam = F_masks()
    blocks = [0x1F, 0x3E0]
    boxes = decompose({(s & 0x1F, s & 0x3E0) for s in fam}, 2)
    big = max(range(len(boxes)), key=lambda i: len(boxes[i][0]))
    broken = list(boxes)
    broken[big] = (frozenset(sorted(boxes[big][0])[1:]),) + boxes[big][1:]
    refusals = 0
    for bad in (broken, boxes + [boxes[0]]):  # an element removed; a box listed twice
        try:
            box_check(fam, blocks, boxes=bad)
        except AssertionError:
            refusals += 1
    rep.check("control: a box split with an element removed, or with a box repeated, is refused", refusals == 2)
    return rep.done()


ROWS = {"irsi_v10": lambda r, c, w: row_irsi(r, "irsi_v10", c, w),
        "irsi_v11": lambda r, c, w: row_irsi(r, "irsi_v11", c, w),
        "irsi_v12": lambda r, c, w: row_irsi(r, "irsi_v12", c, w),
        "irsi_degkm": row_degkm, "cla_F": row_cla_F, "cla_FH": row_cla_FH, "syl_unit": row_syl_unit,
        "syl_Bprime": row_syl_Bprime, "kos": row_kos, "lemma": row_lemma, "selftest": row_selftest}


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in ROWS:
        print(__doc__)
        sys.exit(2)
    workers = os.cpu_count() or 1
    for a in sys.argv:
        if a.startswith("--workers="):
            workers = int(a.split("=")[1])
    ok = ROWS[sys.argv[1]](sys.argv[2], "--controls" in sys.argv, workers)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
