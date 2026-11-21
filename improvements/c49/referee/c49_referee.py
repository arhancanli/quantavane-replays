#!/usr/bin/env python3
"""Independent exact referee for Problem 49a (Tilley PR #379 54-point block) and the X* variant.

Method (exact, complete):
  * The 54 points split into six blocks: left core (7 pts: K5 edges minus 03,04 and minus the
    deleted inner port q), left pendant 1 (10 pts), left pendant 2 (10), and the same three on the
    right (right core drops the deleted inner port p).
  * The family is a union of four products of half-families, Pi_i = U_i x V_i:
        Pi_1 = {x\\q : x in B, q in x} x {y : y in B, p notin y}     (glue, left carries q)
        Pi_2 = {x : x in B, q notin x} x {y\\p : y in B, p in y}     (glue, right carries p)
        Pi_3 = P x L(Ap)          Pi_4 = L(Aq) x P
    and each half-family is a union of 'boxes' (core singleton) x (pendant-1 family) x (pendant-2 family).
  * A triple (A,B,C) is a sunflower iff on EVERY block the projections are intersection-equal
    (a&b == a&c == b&c) and A,B,C are pairwise distinct.  Record, per block triple of families, the set
    of achievable 'difference patterns' (which of the pairs AB, AC, BC differ) over intersection-equal
    triples.  Because every box and every product is a Cartesian product over disjoint coordinates,
    the achievable patterns of a product triple are exactly the pairwise unions of the per-part
    patterns (operator OPLUS).  A sunflower exists iff some ordered triple of products reaches the
    pattern 'all three pairs differ'.  Sperner is decided the same way with the relation set
    {EQ, STRICT} for a subset-of-c pairs.  Nothing is sampled: every ordered box triple is visited.
  * Witnesses for refusals are rebuilt as explicit 54-point sets and re-checked against an explicit
    membership oracle built from the definitions by direct gluing (a separate code path).
"""
import sys
import itertools
import time
from c49_base import (build_F, build_H, A1, A2, SWAP_13, SWAP_23, XSTAR_CORE, PEND_FULL,
                      PORT_P, PORT_Q, E1, E2, permute_mask, eset)

# ------------------------------------------------------------------------------------------------
# pattern algebra
# ------------------------------------------------------------------------------------------------
FULLCODE = 7          # bit0: A!=B, bit1: A!=C, bit2: B!=C


def code3(a, b, c):
    return (a != b) | ((a != c) << 1) | ((b != c) << 2)


_OPLUS = {}


def oplus(s, t):
    key = (s, t)
    r = _OPLUS.get(key)
    if r is None:
        r = 0
        for i in range(8):
            if s >> i & 1:
                for j in range(8):
                    if t >> j & 1:
                        r |= 1 << (i | j)
        _OPLUS[key] = r
    return r


EQ, STRICT = 1, 2


def comb_rel(r, s):
    out = 0
    if (r & EQ) and (s & EQ):
        out |= EQ
    if ((r & STRICT) and s) or ((s & STRICT) and r):
        out |= STRICT
    return out


_FP_CACHE = {}
_FR_CACHE = {}


def fam_pattern(FA, FB, FC):
    key = (FA, FB, FC)
    r = _FP_CACHE.get(key)
    if r is None:
        r = _fam_pattern(FA, FB, FC)
        _FP_CACHE[key] = r
    return r


def fam_rel(FA, FC):
    key = (FA, FC)
    r = _FR_CACHE.get(key)
    if r is None:
        r = _fam_rel(FA, FC)
        _FR_CACHE[key] = r
    return r


def _fam_pattern(FA, FB, FC):
    """mask of codes over intersection-equal triples (a,b,c) in FA x FB x FC, plus one witness per code."""
    mask = 0
    wit = {}
    for a in FA:
        for b in FB:
            ab = a & b
            for c in FC:
                if (a & c) == ab and (b & c) == ab:
                    k = code3(a, b, c)
                    if not (mask >> k & 1):
                        mask |= 1 << k
                        wit[k] = (a, b, c)
    return mask, wit


def _fam_rel(FA, FC):
    r = 0
    wit = {}
    for a in FA:
        for c in FC:
            if a & c == a:
                k = EQ if a == c else STRICT
                if not (r & k):
                    r |= k
                    wit[k] = (a, c)
    return r, wit


class PendTables:
    def __init__(self, fams):
        self.fams = fams
        n = len(fams)
        self.PT, self.PW, self.RT, self.RW = {}, {}, {}, {}
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    m, w = fam_pattern(fams[i], fams[j], fams[k])
                    self.PT[(i, j, k)] = m
                    self.PW[(i, j, k)] = w
                r, w2 = fam_rel(fams[i], fams[j])
                self.RT[(i, j)] = r
                self.RW[(i, j)] = w2


# ------------------------------------------------------------------------------------------------
# half level: union of boxes (core singleton, pend1 family id, pend2 family id)
# ------------------------------------------------------------------------------------------------
def half_patterns(prod_boxes, T):
    flat = [(t, b) for t, lst in enumerate(prod_boxes) for b in lst]
    PT = T.PT
    acc = {}
    for (ta, (ca, a1, a2)) in flat:
        for (tb, (cb, b1, b2)) in flat:
            ab = ca & cb
            dab = 1 if ca != cb else 0
            for (tc, (cc, c1, c2)) in flat:
                if (ca & cc) != ab or (cb & cc) != ab:
                    continue
                m1 = PT[(a1, b1, c1)]
                if not m1:
                    continue
                m2 = PT[(a2, b2, c2)]
                if not m2:
                    continue
                code = dab | ((ca != cc) << 1) | ((cb != cc) << 2)
                m = oplus(oplus(1 << code, m1), m2)
                key = (ta, tb, tc)
                acc[key] = acc.get(key, 0) | m
    return acc


def half_rels(prod_boxes, T):
    flat = [(t, b) for t, lst in enumerate(prod_boxes) for b in lst]
    RT = T.RT
    acc = {}
    for (ta, (ca, a1, a2)) in flat:
        for (tc, (cc, c1, c2)) in flat:
            if ca & cc != ca:
                continue
            r0 = EQ if ca == cc else STRICT
            r1 = RT[(a1, c1)]
            if not r1:
                continue
            r2 = RT[(a2, c2)]
            if not r2:
                continue
            r = comb_rel(comb_rel(r0, r1), r2)
            key = (ta, tc)
            acc[key] = acc.get(key, 0) | r
    return acc


def box_size(box, T):
    return len(T.fams[box[1]]) * len(T.fams[box[2]])


def box_overlaps(boxes, T):
    """pairs of distinct boxes in one half-family sharing a member (would break counting)."""
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (ca, a1, a2), (cc, c1, c2) = boxes[i], boxes[j]
            if ca == cc and (T.RT[(a1, c1)] & EQ) and (T.RT[(a2, c2)] & EQ):
                bad.append((i, j))
    return bad


# ------------------------------------------------------------------------------------------------
# witness reconstruction (explicit sets, layout: core | pend1<<10 | pend2<<20 ; full = L | R<<30)
# ------------------------------------------------------------------------------------------------
def explicit_half(core, x1, x2):
    return core | (x1 << 10) | (x2 << 20)


def half_triple_witness(boxlists, T, target):
    """Find boxes (from three lists) and explicit half sets realising half-pattern code == target."""
    LA, LB, LC = boxlists
    for (ca, a1, a2) in LA:
        for (cb, b1, b2) in LB:
            for (cc, c1, c2) in LC:
                ab = ca & cb
                if (ca & cc) != ab or (cb & cc) != ab:
                    continue
                k0 = code3(ca, cb, cc)
                m1 = T.PT[(a1, b1, c1)]
                m2 = T.PT[(a2, b2, c2)]
                for k1 in range(8):
                    if not (m1 >> k1 & 1):
                        continue
                    for k2 in range(8):
                        if not (m2 >> k2 & 1):
                            continue
                        if k0 | k1 | k2 == target:
                            w1 = T.PW[(a1, b1, c1)][k1]
                            w2 = T.PW[(a2, b2, c2)][k2]
                            return (explicit_half(ca, w1[0], w2[0]), explicit_half(cb, w1[1], w2[1]),
                                    explicit_half(cc, w1[2], w2[2]))
    return None


def half_pair_witness(LA, LC, T, target_rel):
    for (ca, a1, a2) in LA:
        for (cc, c1, c2) in LC:
            if ca & cc != ca:
                continue
            r0 = EQ if ca == cc else STRICT
            r1, r2 = T.RT[(a1, c1)], T.RT[(a2, c2)]
            for k1 in (EQ, STRICT):
                if not (r1 & k1):
                    continue
                for k2 in (EQ, STRICT):
                    if not (r2 & k2):
                        continue
                    got = STRICT if STRICT in (r0, k1, k2) else EQ
                    if got == target_rel:
                        w1, w2 = T.RW[(a1, c1)][k1], T.RW[(a2, c2)][k2]
                        return explicit_half(ca, w1[0], w2[0]), explicit_half(cc, w1[1], w2[1])
    return None


# ------------------------------------------------------------------------------------------------
# full decision over products
# ------------------------------------------------------------------------------------------------
def decide(products, T, want_witness=True):
    """products: list of (left_boxes, right_boxes). Returns dict with verdicts and witnesses."""
    U = [p[0] for p in products]
    V = [p[1] for p in products]
    n = len(products)
    PL = half_patterns(U, T)
    PR = half_patterns(V, T)
    RL = half_rels(U, T)
    RR = half_rels(V, T)
    sf_hits = []
    for key in itertools.product(range(n), repeat=3):
        tot = oplus(PL.get(key, 0), PR.get(key, 0))
        if tot >> FULLCODE & 1:
            sf_hits.append(key)
    sp_hits, overlaps = [], []
    for key in itertools.product(range(n), repeat=2):
        tot = comb_rel(RL.get(key, 0), RR.get(key, 0))
        if tot & STRICT:
            sp_hits.append(key)
        if key[0] < key[1] and (tot & EQ):
            overlaps.append(key)
    box_ov = []
    for i in range(n):
        for side, L in (("L", U[i]), ("R", V[i])):
            for pr in box_overlaps(L, T):
                box_ov.append((i, side, pr))
    res = {"sunflower_product_triples": sf_hits, "sperner_product_pairs": sp_hits,
           "product_overlaps": overlaps, "box_overlaps": box_ov}
    if want_witness and sf_hits:
        key = sf_hits[0]
        found = None
        for s in range(8):
            if not (PL.get(key, 0) >> s & 1):
                continue
            for t in range(8):
                if (PR.get(key, 0) >> t & 1) and (s | t) == FULLCODE:
                    wl = half_triple_witness([U[key[0]], U[key[1]], U[key[2]]], T, s)
                    wr = half_triple_witness([V[key[0]], V[key[1]], V[key[2]]], T, t)
                    if wl and wr:
                        found = tuple(wl[i] | (wr[i] << 30) for i in range(3))
                        break
            if found:
                break
        res["sf_witness"] = (key, found)
    if want_witness and sp_hits:
        key = sp_hits[0]
        found = None
        for rl in (EQ, STRICT):
            if not (RL.get(key, 0) & rl):
                continue
            for rr in (EQ, STRICT):
                if (RR.get(key, 0) & rr) and STRICT in (rl, rr):
                    wl = half_pair_witness(U[key[0]], U[key[1]], T, rl)
                    wr = half_pair_witness(V[key[0]], V[key[1]], T, rr)
                    if wl and wr:
                        found = (wl[0] | (wr[0] << 30), wl[1] | (wr[1] << 30))
                        break
            if found:
                break
        res["sp_witness"] = (key, found)
    return res


def counts(products, T, lport, rport):
    """N[s][t] with s = left port bit, t = right port bit (ports are core bits)."""
    N = [[0, 0], [0, 0]]
    for (Ub, Vb) in products:
        cu = [0, 0]
        cv = [0, 0]
        for b in Ub:
            cu[1 if b[0] & lport else 0] += box_size(b, T)
        for b in Vb:
            cv[1 if b[0] & rport else 0] += box_size(b, T)
        for s in (0, 1):
            for t in (0, 1):
                N[s][t] += cu[s] * cv[t]
    return N


# ------------------------------------------------------------------------------------------------
# construction (from the note's page 2)
# ------------------------------------------------------------------------------------------------
F = build_F()
H = build_H()
H_IN = tuple(sorted(Y >> 1 for Y in H if Y & 1))       # Y\{0}, Y contains 0   (glued edge absent)
H_OUT = tuple(sorted(Y >> 1 for Y in H if not Y & 1))  # Y, 0 not in Y         (glued edge present)
assert len(H_IN) == 75 and len(H_OUT) == 46


class Fams:
    """pendant family registry"""

    def __init__(self):
        self.fams = []
        self.idx = {}

    def id(self, fam):
        fam = tuple(sorted(set(fam)))
        if fam not in self.idx:
            self.idx[fam] = len(self.fams)
            self.fams.append(fam)
        return self.idx[fam]


def L_boxes(X, R):
    """L(X): glue an H at each pendant edge (exactly-one rule with point 0)."""
    out = []
    for x in X:
        f1 = R.id(H_OUT) if x & E1 else R.id(H_IN)
        f2 = R.id(H_OUT) if x & E2 else R.id(H_IN)
        out.append((x & ~(E1 | E2), f1, f2))
    return out


def B_boxes(R, xstar=None, extra=()):
    """28-point block B. xstar in (None, 1, 2): pendant that X* fills. extra: list of boxes."""
    E = R.id([0])
    boxes = L_boxes(F, R)
    for x in A1:      # only the second pendant glued; 03 absent in every A1 mask (checked)
        assert not x & E1
        boxes.append((x & ~E2, E, R.id(H_OUT) if x & E2 else R.id(H_IN)))
    for x in A2:      # only the first pendant glued; 04 absent in every A2 mask (checked)
        assert not x & E2
        boxes.append((x & ~E1, R.id(H_OUT) if x & E1 else R.id(H_IN), E))
    # xstar: None, 1 (X' : pendant 1 full), 2 (X* : pendant 2 full), 12 (both X* and X')
    if xstar in (2, 12):
        boxes.append((XSTAR_CORE, E, R.id([PEND_FULL])))
    if xstar in (1, 12):
        boxes.append((XSTAR_CORE, R.id([PEND_FULL]), E))
    boxes.extend(extra)
    return boxes


def build54(R, xstar_left=None, xstar_right=None, perm_p=SWAP_13, perm_q=SWAP_23,
            extra_products=(), extraB_left=(), extraB_right=(), Ap_masks=None, Aq_masks=None):
    E = R.id([0])
    HI = R.id(H_IN)
    BL = B_boxes(R, xstar_left, extraB_left)
    BR = B_boxes(R, xstar_right, extraB_right)
    U1 = [(c & ~PORT_Q, a, b) for (c, a, b) in BL if c & PORT_Q]
    U2 = [(c, a, b) for (c, a, b) in BL if not c & PORT_Q]
    V1 = [(c, a, b) for (c, a, b) in BR if not c & PORT_P]
    V2 = [(c & ~PORT_P, a, b) for (c, a, b) in BR if c & PORT_P]
    Pb = [(0, HI, E), (0, E, HI)]
    Ap = Ap_masks if Ap_masks is not None else [permute_mask(x, perm_p) for x in A1]
    Aq = Aq_masks if Aq_masks is not None else [permute_mask(x, perm_q) for x in A1]
    LAp = L_boxes(Ap, R)
    LAq = L_boxes(Aq, R)
    struct_ok = all(not c & PORT_P for (c, _, _) in LAp) and all(not c & PORT_Q for (c, _, _) in LAq)
    products = [(U1, V1), (U2, V2), (Pb, LAp), (LAq, Pb)] + list(extra_products)
    return products, struct_ok, BL, BR


# ------------------------------------------------------------------------------------------------
# explicit oracle (separate code path: direct gluing with H, no H_IN/H_OUT shortcut)
# ------------------------------------------------------------------------------------------------
def glue_explicit(X, e1_glued=True, e2_glued=True):
    out = set()
    for x in X:
        core = x & ~(E1 | E2)
        Y1s = [Y for Y in H if (bool(x & E1) != bool(Y & 1))] if e1_glued else [None]
        Y2s = [Y for Y in H if (bool(x & E2) != bool(Y & 1))] if e2_glued else [None]
        for Y1 in Y1s:
            for Y2 in Y2s:
                s = core
                if Y1 is not None:
                    s |= (Y1 >> 1) << 10
                if Y2 is not None:
                    s |= (Y2 >> 1) << 20
                out.add(s)
    return out


def B_explicit(xstar=None, extra_sets=()):
    S = glue_explicit(F)
    S |= glue_explicit(A1, e1_glued=False)
    S |= glue_explicit(A2, e2_glued=False)
    if xstar in (2, 12):
        S.add(XSTAR_CORE | (PEND_FULL << 20))
    if xstar in (1, 12):
        S.add(XSTAR_CORE | (PEND_FULL << 10))
    S |= set(extra_sets)
    return S


class Oracle:
    def __init__(self, xstar_left, xstar_right, perm_p=SWAP_13, perm_q=SWAP_23, extra_members=(),
                 extraB_left=(), extraB_right=(), Ap_masks=None, Aq_masks=None):
        self.BL = B_explicit(xstar_left, extraB_left)
        self.BR = B_explicit(xstar_right, extraB_right)
        self.P = {(Y >> 1) << 10 for Y in H if Y & 1} | {(Y >> 1) << 20 for Y in H if Y & 1}
        Ap = Ap_masks if Ap_masks is not None else [permute_mask(x, perm_p) for x in A1]
        Aq = Aq_masks if Aq_masks is not None else [permute_mask(x, perm_q) for x in A1]
        self.LAp = glue_explicit(Ap)
        self.LAq = glue_explicit(Aq)
        self.extra = set(extra_members)

    def member(self, S):
        if S in self.extra:
            return True
        Lh = S & ((1 << 30) - 1)
        Rh = S >> 30
        if Lh & PORT_Q or Rh & PORT_P:
            return False
        if (Lh | PORT_Q) in self.BL and Rh in self.BR:
            return True
        if Lh in self.BL and (Rh | PORT_P) in self.BR:
            return True
        if Lh in self.P and Rh in self.LAp:
            return True
        if Lh in self.LAq and Rh in self.P:
            return True
        return False

    def counts(self):
        def split(S, port):
            c = [0, 0]
            for s in S:
                c[1 if s & port else 0] += 1
            return c
        N = [[0, 0], [0, 0]]
        # glue: left x (bits p,q), right y (bits p,q); exactly one of q in x, p in y
        cl = {}
        for x in self.BL:
            k = (1 if x & PORT_P else 0, 1 if x & PORT_Q else 0)
            cl[k] = cl.get(k, 0) + 1
        cr = {}
        for y in self.BR:
            k = (1 if y & PORT_P else 0, 1 if y & PORT_Q else 0)
            cr[k] = cr.get(k, 0) + 1
        for s in (0, 1):
            for t in (0, 1):
                N[s][t] += cl.get((s, 1), 0) * cr.get((0, t), 0) + cl.get((s, 0), 0) * cr.get((1, t), 0)
        lap = split(self.LAp, PORT_Q)
        laq = split(self.LAq, PORT_P)
        for t in (0, 1):
            N[0][t] += len(self.P) * lap[t]
        for s in (0, 1):
            N[s][0] += laq[s] * len(self.P)
        return N, cl, cr, lap, laq


def check_sf_witness(o, trip):
    A, B, C = trip
    ok_members = all(o.member(X) for X in trip)
    distinct = A != B and A != C and B != C
    ieq = (A & B) == (A & C) == (B & C)
    return ok_members and distinct and ieq


def check_sp_witness(o, pair):
    A, C = pair
    return o.member(A) and o.member(C) and A != C and (A & C) == A
