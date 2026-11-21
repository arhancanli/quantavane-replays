#!/usr/bin/env python3
"""Program A: exact check of the 54-point block of R. F. Tilley Jr., PR #379 of teorth/optimizationproblems (49a),
from the construction on page 2 of his note (Zenodo 10.5281/zenodo.24405229, c49-note.pdf). His verifier, Lean files
and scripts were neither read nor run; the note was read as mathematics.

Primitives from our own 4 November 49a checker (quantavane-replays c49/sunflower_check.py, imported unchanged):
family_F, family_H, ws_patterns, sub_patterns, decompose, Surd, check_printed, dec_str.

New here: a two-level box method. The 54 points are two halves (left copy minus q, right copy minus p); each half is
cut into three blocks (core, pendant 1, pendant 2). A family on the 54 points is a union of top-level boxes X x Y
(X a family on the left half, Y on the right); each half family is split exactly into sub-boxes over its three blocks.
Three members form a 3-sunflower iff on every block their projections x, y, z satisfy x&y == x&z == y&z, and the
members are pairwise distinct (some block separates each pair). For a triple of sub-boxes the separation patterns that
occur on the half are the ORs of one admissible pattern per block (the blocks are independent inside a box); the union
over all sub-box triples is the half's pattern set; a triple of top boxes holds a 3-sunflower iff some left pattern OR
some right pattern is 7 (all three pairs separated). Containment likewise with (equal, proper) relations. Exhaustive:
every multiset of three top boxes, every triple of sub-boxes, every member triple inside them through ws_patterns.

Usage: python3 -I tilley_a.py C49_CHECKER [--control plant-cross|plant-contain|drop-A1]
"""
import importlib.util
import itertools
import json
import sys
import time
from fractions import Fraction

args = sys.argv[1:]
spec = importlib.util.spec_from_file_location("sunflower_check", args[0])
sc = importlib.util.module_from_spec(spec)
sys.modules["sunflower_check"] = sc
spec.loader.exec_module(sc)


def opt(name, default=None):
    return args[args.index(name) + 1] if name in args else default


BLOCK = opt("--block", "tilley")
CONTROL = opt("--control", "")

# ---- the base families, in the note's terms ------------------------------------------------------------------
K5 = [(a, b) for a in range(5) for b in range(a + 1, 5)]          # bit i of an edge mask = K5[i] (01, 02, 03, ...)
CORE = [(0, 1), (0, 2), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]  # the 8 core points; p = 01 (0), q = 02 (1)
E1, E2 = (0, 3), (0, 4)
CI = {e: i for i, e in enumerate(CORE)}
A1_MASKS = [251, 425, 499, 506, 570, 707, 761, 875, 883, 889, 923, 947, 985, 986, 1002]
A2_MASKS = [358, 471, 499, 501, 631, 707, 757, 758, 789, 883, 935, 947, 950, 982, 997]


def edges_of(m):
    return frozenset(K5[i] for i in range(10) if m >> i & 1)


H = sc.family_H()
H_AVOID = [sc.mask(x - 1 for x in D) for D in H if 0 not in D]      # glued edge present: an H member avoiding 0
H_CONT = [sc.mask(x - 1 for x in D if x) for D in H if 0 in D]       # glued edge absent: Y - {0}, 0 in Y
assert len(H_AVOID) == 46 and len(H_CONT) == 75


def core_mask(T):
    return sc.mask(CI[e] for e in T if e not in (E1, E2))


def L_of(X, pendants=(1, 2)):
    """Glue an H at e1 (pendant 1, points 8..17) and/or e2 (pendant 2, points 18..27) of every core member."""
    out = []
    for T in X:
        s1 = (H_AVOID if E1 in T else H_CONT) if 1 in pendants else [0]
        s2 = (H_AVOID if E2 in T else H_CONT) if 2 in pendants else [0]
        c = core_mask(T)
        out += [c | a << 8 | b << 18 for a in s1 for b in s2]
    return out


def swap(T, u, v):
    f = lambda x: v if x == u else u if x == v else x
    return frozenset((min(f(a), f(b)), max(f(a), f(b))) for a, b in T)


F = [frozenset(S) for S in sc.family_F()]
A1 = [edges_of(m) for m in (A1_MASKS[1:] if CONTROL == "drop-A1" else A1_MASKS)]
A2 = [edges_of(m) for m in A2_MASKS]
assert all(E1 not in T for T in A1) and all(E2 not in T for T in A2)
B_tilley = L_of(F) + L_of(A1, (2,)) + L_of(A2, (1,))


def read_pinho(path):
    with open(path, "rb") as fh:
        d = json.loads(fh.read().decode("ascii"))
    assert d["n"] == 28 and d["ports"] == [0, 1]
    return d["masks"]


def read_koshta(path):
    return sc.read_sunflower3(path)[0][3]


if BLOCK == "tilley":
    B = B_tilley
elif BLOCK == "pinho":
    B = read_pinho(opt("--pinho-cert"))
elif BLOCK == "koshta":
    B = read_koshta(opt("--koshta-cert"))
else:
    raise SystemExit("unknown block")
assert len(set(B)) == len(B)

P = [h << 8 for h in H_CONT] + [h << 18 for h in H_CONT]
Ap = [swap(T, 1, 3) for T in A1]
Aq = [swap(T, 2, 3) for T in A1]
LAp, LAq = L_of(Ap), L_of(Aq)
PBIT, QBIT = 1 << CI[(0, 1)], 1 << CI[(0, 2)]

out = []
def line(s):
    print(s, flush=True)
    out.append(s)


line("PROGRAM A (two-level box method), block = %s, control = %s" % (BLOCK, CONTROL or "none"))
N = {(a, b): sum(1 for m in B if bool(m & PBIT) == a and bool(m & QBIT) == b) for a in (0, 1) for b in (0, 1)}
line("block B: %d members on 28 points; port counts (p,q) 00 %d, 01 %d, 10 %d, 11 %d"
     % (len(B), N[0, 0], N[0, 1], N[1, 0], N[1, 1]))
line("|P| = %d, |L(Ap)| = %d, |L(Aq)| = %d; L(Ap) avoids p: %s; L(Aq) avoids q: %s"
     % (len(P), len(LAp), len(LAq), all(not m & PBIT for m in LAp), all(not m & QBIT for m in LAq)))

BL, BR, PL, PR = B, B, P, P
# ---- half families (28-bit coordinates of their copy; the deleted inner port is never present) -------------------
L0 = [m for m in BL if not m & QBIT]
L1 = [m & ~QBIT for m in BL if m & QBIT]
R0 = [m for m in BR if not m & PBIT]
R1 = [m & ~PBIT for m in BR if m & PBIT]
TOP = [("G: L0 x R1", L0, R1), ("G: L1 x R0", L1, R0), ("P x L(Ap)", PL, LAp), ("L(Aq) x P", LAq, PR)]

if CONTROL == "plant-cross":
    # red control: one extra top box {x} x {y} that makes a 3-sunflower with two members of P x L(Ap):
    # take a, b in P with a & b == 0 and c in L(Ap); the set (0, c) is distinct and (a,c),(b,c),(0,c) pairwise meet in c
    TOP.append(("planted", [0], [LAp[0]]))
if CONTROL == "plant-contain":
    # red control: (P[0] minus its lowest point, L(Ap)[0]) lies strictly inside the member (P[0], L(Ap)[0])
    TOP.append(("planted", [P[0] & (P[0] - 1)], [LAp[0]]))
if CONTROL == "drop-A1":
    # red control: the note's counts must fail when one member of A1 is left out of the block
    pass

BLOCKS = [sc.mask(range(8)), sc.mask(range(8, 18)), sc.mask(range(18, 28))]


def subboxes(fam):
    boxes = sc.decompose({tuple(s & bm for bm in BLOCKS) for s in fam}, 3)
    total = 0
    seen = set()
    fs = set(fam)
    for box in boxes:
        for parts in itertools.product(*box):
            s = parts[0] | parts[1] | parts[2]
            assert s in fs and s not in seen
            seen.add(s)
            total += 1
    assert total == len(fam) == len(fs)
    return boxes


SIDE = {}
def sid(S):
    return SIDE.setdefault(S, len(SIDE))

SIDES = {}
WS = {}
SP = {}
def ws(a, b, c):
    k = (a, b, c)
    if k not in WS:
        WS[k] = sc.ws_patterns(SIDES[a], SIDES[b], SIDES[c])
    return WS[k]


def sp(a, b):
    k = (a, b)
    if k not in SP:
        SP[k] = sc.sub_patterns(SIDES[a], SIDES[b])
    return SP[k]


OR_SET = {}
def or_combine(masks):
    """All ORs of one pattern from each 8-bit pattern mask, as an 8-bit mask."""
    reach = 1
    for m in masks:
        new = 0
        for r in range(8):
            if reach >> r & 1:
                for s in range(8):
                    if m >> s & 1:
                        new |= 1 << (r | s)
        reach = new
    return reach


HALF_FAMS = {}
def half_id(fam):
    key = id(fam)
    if key not in HALF_FAMS:
        bx = subboxes(fam)
        ids = []
        for box in bx:
            t = []
            for S in box:
                i = sid(S)
                SIDES[i] = S
                t.append(i)
            ids.append(tuple(t))
        HALF_FAMS[key] = ids
    return HALF_FAMS[key]


HP = {}
def half_patterns(X, Y, Z):
    k = (id(X), id(Y), id(Z))
    if k in HP:
        return HP[k]
    bx, by, bz = half_id(X), half_id(Y), half_id(Z)
    reach = 0
    cache = {}
    for u in bx:
        for v in by:
            for w in bz:
                ms = (ws(u[0], v[0], w[0]), ws(u[1], v[1], w[1]), ws(u[2], v[2], w[2]))
                if 0 in ms:
                    continue
                if ms not in cache:
                    cache[ms] = or_combine(ms)
                reach |= cache[ms]
                if reach == 255:
                    break
            if reach == 255:
                break
        if reach == 255:
            break
    HP[k] = reach
    return reach


HR = {}
def half_relations(X, Y):
    """bit 0: some x in X equals some y in Y on the whole half; bit 1: some x is a proper subset of some y."""
    k = (id(X), id(Y))
    if k in HR:
        return HR[k]
    bx, by = half_id(X), half_id(Y)
    res = 0
    for u in bx:
        for v in by:
            rel = [sp(u[i], v[i]) for i in range(3)]
            if all(rel):
                if all(r & 1 for r in rel):
                    res |= 1
                if any(r & 2 for r in rel):
                    res |= 2
    HR[k] = res
    return res


t0 = time.time()
for name, X, Y in TOP:
    line("top box %-12s: %d x %d = %d members; sub-boxes %d (left) and %d (right)"
         % (name, len(X), len(Y), len(X) * len(Y), len(half_id(X)), len(half_id(Y))))
size = sum(len(X) * len(Y) for _, X, Y in TOP)
cols = {}
for name, X, Y in TOP:
    for a in (0, 1):
        for b in (0, 1):
            cols[name, a, b] = sum(1 for x in X if bool(x & PBIT) == a) * sum(1 for y in Y if bool(y & QBIT) == b)
for name, _, _ in TOP:
    line("  %-12s columns (p_left, q_right) 00 %d  01 %d  10 %d  11 %d"
         % (name, cols[name, 0, 0], cols[name, 0, 1], cols[name, 1, 0], cols[name, 1, 1]))
M = {(a, b): sum(cols[name, a, b] for name, _, _ in TOP) for a in (0, 1) for b in (0, 1)}
line("54-point block: %d members; totals 00 %d  01 %d  10 %d  11 %d" % (size, M[0, 0], M[0, 1], M[1, 0], M[1, 1]))

# containment: ordered pairs of top boxes
contain = []
for i, (n1, X1, Y1) in enumerate(TOP):
    for j, (n2, X2, Y2) in enumerate(TOP):
        rl, rr = half_relations(X1, X2), half_relations(Y1, Y2)
        if rl and rr and ((rl & 2) or (rr & 2)):
            contain.append((n1, n2))
line("containments between top boxes (first inside second): %s" % (contain or "none"))
# sunflowers: multisets of three top boxes
sun = []
for i, j, k in itertools.combinations_with_replacement(range(len(TOP)), 3):
    pl = half_patterns(TOP[i][1], TOP[j][1], TOP[k][1])
    pr = half_patterns(TOP[i][2], TOP[j][2], TOP[k][2])
    if any((pl >> a & 1) and (pr >> b & 1) and (a | b) == 7 for a in range(8) for b in range(8)):
        sun.append((TOP[i][0], TOP[j][0], TOP[k][0]))
line("3-sunflowers among top-box triples: %s" % (sun or "none"))
line("exhaustive check time %.1f s; half pattern triples computed %d" % (time.time() - t0, len(HP)))
valid = not contain and not sun
line("VALID (Sperner and 3-sunflower-free): %s" % valid)

# ---- rate --------------------------------------------------------------------------------------------------------
a, b, c = M[0, 1], M[0, 0], M[1, 1]
symmetric = M[0, 1] == M[1, 0]
line("a = N01 = N10: %s (%d); b = N00 = %d; c = N11 = %d" % (symmetric, a, b, c))
S = sc.Surd(((a, b), (c, a)))
line("Perron root rho = %d + sqrt(%d) (b*c = %d)" % (a, b * c, b * c))
assert S.t == 2 * a and S.D == 4 * b * c
r20 = S.rate_floor(52, 20)
line("rho^(1/52) in [%s, %s + 1e-20)" % (sc.dec_str(r20, 20), sc.dec_str(r20, 20)))
print("RATE20 %s VALID %s" % (sc.dec_str(r20, 20), valid))
for printed in ("1.567866388204975", "1.56786638820497592522"):
    v, tr, k = sc.check_printed(S, 52, printed)
    line("printed %s: valid lower bound %s; truncation at %d decimals %s" % (printed, v, k, tr))
if BLOCK == "tilley":
    ok_counts = (len(B) == 252271 and (N[0, 0], N[0, 1], N[1, 0], N[1, 1]) == (33948, 64728, 64728, 88867)
                 and len(P) == 150 and len(LAp) == len(LAq) == 42412
                 and (M[0, 0], M[0, 1], M[1, 0], M[1, 1]) == (4397311488, 7211663100, 7211663100, 11504366352)
                 and b * c == 50588282321810251776 and a == 7211663100)
    line("all counts as printed on page 2 and rho as in the Theorem: %s" % ok_counts)
    valid = valid and ok_counts and sc.check_printed(S, 52, "1.567866388204975")[0]
print("RESULT %s" % ("PASS" if valid else "FAIL"))
sys.exit(0 if valid else 1)
