#!/usr/bin/env python3
"""Referee run: Tilley's 54-point family, the X* variant, red controls, counts and rate."""
import sys
import time
import itertools
from c49_base import (A1, A2, SWAP_13, SWAP_23, XSTAR_CORE, PEND_FULL, PORT_P, PORT_Q, E1, E2,
                      permute_mask, eset, ebit)
from c49_referee import (F, H, H_IN, H_OUT, Fams, PendTables, decide, counts, B_boxes, build54,
                         box_size, Oracle, B_explicit, check_sf_witness, check_sp_witness)

T0 = time.time()


def log(*a):
    print("[%6.1fs]" % (time.time() - T0), *a, flush=True)


def brute_valid(fam):
    fam = list(fam)
    sp = any(a != c and a & c == a for a in fam for c in fam)
    sf = any((a & b) == (a & c) == (b & c) for a, b, c in itertools.combinations(fam, 3))
    return (not sp) and (not sf)


# ---- 0. base families (brute force) -------------------------------------------------------------
log("H: 121 members, Sperner+SF-free on Z11 (brute):", brute_valid(H))
log("F: 75 members, Sperner+SF-free on K5 edges (brute):", brute_valid(F))

# ---- registry + tables ---------------------------------------------------------------------------
R = Fams()
E = R.id([0])
R.id(H_IN)
R.id(H_OUT)
R.id([PEND_FULL])


def tables():
    return PendTables(R.fams)


def block_report(name, xstar):
    T = tables()
    Bb = B_boxes(R, xstar)
    T = tables()
    triv = [(0, E, E)]
    res = decide([(Bb, triv)], T)
    c = {}
    for b in Bb:
        k = (1 if b[0] & PORT_P else 0, 1 if b[0] & PORT_Q else 0)
        c[k] = c.get(k, 0) + box_size(b, T)
    tot = sum(c.values())
    Bx = B_explicit(xstar)
    cx = {}
    for s in Bx:
        k = (1 if s & PORT_P else 0, 1 if s & PORT_Q else 0)
        cx[k] = cx.get(k, 0) + 1
    log("%s: boxes=%d |B|=%d (explicit %d) (p,q) 00 %d 01 %d 10 %d 11 %d | explicit %s" % (
        name, len(Bb), tot, len(Bx), c.get((0, 0), 0), c.get((0, 1), 0), c.get((1, 0), 0), c.get((1, 1), 0),
        [cx.get(k, 0) for k in ((0, 0), (0, 1), (1, 0), (1, 1))]))
    log("   sunflower=%s sperner_violation=%s box_overlaps=%s" % (
        bool(res["sunflower_product_triples"]), bool(res["sperner_product_pairs"]), res["box_overlaps"]))
    return res


def run54(name, expect_valid, **kw):
    products, struct_ok, BL, BR = build54(R, **{k: v for k, v in kw.items() if not k.startswith("_")})
    T = tables()
    t = time.time()
    res = decide(products, T)
    N = counts(products, T, PORT_P, PORT_Q)
    sf = bool(res["sunflower_product_triples"])
    sp = bool(res["sperner_product_pairs"])
    valid = struct_ok and not sf and not sp
    log("%s: struct_ok=%s sunflower=%s sperner_violation=%s product_overlaps=%s box_overlaps=%s (%.1fs)" % (
        name, struct_ok, sf, sp, res["product_overlaps"], res["box_overlaps"], time.time() - t))
    log("   counts N00 %d N01 %d N10 %d N11 %d total %d" % (N[0][0], N[0][1], N[1][0], N[1][1],
                                                          sum(map(sum, N))))
    if sf or sp:
        o = Oracle(kw.get("xstar_left"), kw.get("xstar_right"), kw.get("perm_p", SWAP_13),
                   kw.get("perm_q", SWAP_23), extraB_left=kw.get("_xB_left", ()),
                   extraB_right=kw.get("_xB_right", ()), Ap_masks=kw.get("Ap_masks"),
                   Aq_masks=kw.get("Aq_masks"))
        if sf:
            key, w = res["sf_witness"]
            log("   sunflower product triple", key, "(of", len(res["sunflower_product_triples"]), "hits)",
                "explicit witness verified by oracle:", w is not None and check_sf_witness(o, w))
            if w:
                for X in w:
                    log("      member %015x" % X)
        if sp:
            key, w = res["sp_witness"]
            log("   containment product pair", key, "explicit witness verified by oracle:",
                w is not None and check_sp_witness(o, w))
    verdict = "VALID" if valid else "REFUSED"
    if expect_valid is None:
        log("   verdict:", verdict, "| (informational)")
    else:
        log("   verdict:", verdict, "| expected:", "VALID" if expect_valid else "REFUSED",
            "| OK" if (valid == expect_valid) else "| *** UNEXPECTED ***")
    return res, N, valid


if __name__ == "__main__":
    log("=== 28-point blocks ===")
    block_report("Tilley B", None)
    block_report("B' = B + X* (pendant 2)", 2)

    log("=== Tilley 54-point family ===")
    res, N, ok = run54("Tilley54", True)
    o = Oracle(None, None)
    No, cl, cr, lap, laq = o.counts()
    log("   oracle recount", No, "| L(Ap) q-split (0,1) =", lap, "L(Aq) p-split =", laq, "|P| =", len(o.P))
    exp = [[4397311488, 7211663100], [7211663100, 11504366352]]
    log("   matches note's printed totals:", N == exp and No == exp)

    log("=== Modified: X* added to both copies ===")
    res, N2, ok2 = run54("X*-both", True, xstar_left=2, xstar_right=2)
    o2 = Oracle(2, 2)
    No2 = o2.counts()[0]
    log("   oracle recount", No2)
    exp2 = [[4397440944, 7211751967], [7211751967, 11504366352]]
    log("   matches claimed counts:", N2 == exp2 and No2 == exp2, "| total", sum(map(sum, N2)))
    # disjointness of the four parts, checked explicitly on the half sets
    BLq = {x & ~PORT_Q for x in o2.BL if x & PORT_Q}
    BLn = {x for x in o2.BL if not x & PORT_Q}
    BRp = {y & ~PORT_P for y in o2.BR if y & PORT_P}
    BRn = {y for y in o2.BR if not y & PORT_P}
    parts = [(BLq, BRn), (BLn, BRp), (o2.P, o2.LAp), (o2.LAq, o2.P)]
    inter = [(i, j) for i in range(4) for j in range(i + 1, 4) if (parts[i][0] & parts[j][0]) and (parts[i][1] & parts[j][1])]
    log("   explicit check: the four products pairwise disjoint:", not inter, inter)

    log("=== Second variant: X* AND X' added to both copies ===")
    block_report("B'' = B + X* + X'", 12)
    res, N3, ok3 = run54("X*+X'-both", True, xstar_left=12, xstar_right=12)
    No3 = Oracle(12, 12).counts()[0]
    log("   oracle recount", No3)
    log("   total", sum(map(sum, N3)))

    log("=== Variants (informational) ===")
    run54("X*-left-only", True, xstar_left=2)
    run54("X*-right-only", True, xstar_right=2)

    log("=== Red controls ===")
    # informational: X' alone (pendant 1 full) in both copies
    run54("X'-only (pendant 1) both copies", None, xstar_left=1, xstar_right=1)
    # C2: Sperner plant: add X* minus one pendant point to B (strictly inside X*)
    sub = (XSTAR_CORE, E, R.id([PEND_FULL & ~1]))
    subx = XSTAR_CORE | ((PEND_FULL & ~1) << 20)
    run54("C2 B' + (X* minus a point)", False, xstar_left=2, xstar_right=2, extraB_left=[sub], extraB_right=[sub],
          _xB_left=[subx], _xB_right=[subx])
    # C3: planted sunflower: take two explicit members A,B of the modified family and add
    #     C = (A & B) | (54 points outside A|B)  (a single extra member, as a one-box product)
    A_ = None
    # pick A,B from product 1 (glue) using explicit oracle sets
    xL = sorted(x for x in o2.BL if x & PORT_Q)[1000]
    yR = sorted(y for y in o2.BR if not y & PORT_P)[2000]
    xL2 = sorted(x for x in o2.BL if x & PORT_Q)[5000]
    yR2 = sorted(y for y in o2.BR if not y & PORT_P)[7000]
    A_ = (xL & ~PORT_Q) | (yR << 30)
    B_ = (xL2 & ~PORT_Q) | (yR2 << 30)
    assert o2.member(A_) and o2.member(B_)
    # ground set masks of each half: core bits (minus deleted port) + 20 pendant bits
    allL = (((1 << 10) - 1) & ~(E1 | E2 | PORT_Q)) | (((1 << 20) - 1) << 10)
    allR = (((1 << 10) - 1) & ~(E1 | E2 | PORT_P)) | (((1 << 20) - 1) << 10)
    ALL = allL | (allR << 30)
    C_ = (A_ & B_) | (ALL & ~(A_ | B_))
    CL, CR = C_ & ((1 << 30) - 1), C_ >> 30
    boxL = (CL & ((1 << 10) - 1), R.id([(CL >> 10) & 1023]), R.id([(CL >> 20) & 1023]))
    boxR = (CR & ((1 << 10) - 1), R.id([(CR >> 10) & 1023]), R.id([(CR >> 20) & 1023]))
    products, struct_ok, BL, BR = build54(R, xstar_left=2, xstar_right=2, extra_products=[([boxL], [boxR])])
    T = tables()
    r3 = decide(products, T)
    o3 = Oracle(2, 2, extra_members=[C_])
    sfw = r3.get("sf_witness")
    log("C3 planted sunflower: sunflower=%s witness verified=%s" % (
        bool(r3["sunflower_product_triples"]), sfw is not None and sfw[1] is not None and check_sf_witness(o3, sfw[1])),
        "| verdict:", "REFUSED" if (r3["sunflower_product_triples"] or r3["sperner_product_pairs"]) else "VALID",
        "| expected REFUSED")
    # C4: Ap/Aq from A1 or A2 under EVERY vertex permutation that sends the absent pendant edge onto
    #     the deleted inner port (so the structure is well formed); the note's swaps are two of them.
    for src_name, src, absent in (("A1", A1, (0, 3)), ("A2", A2, (0, 4))):
        for which, port in (("Ap", (0, 1)), ("Aq", (0, 2))):
            nval = 0
            tot = 0
            for perm in itertools.permutations(range(5)):
                sig = list(perm)
                img = tuple(sorted((sig[absent[0]], sig[absent[1]])))
                if img != port:
                    continue
                masks = [permute_mask(x, sig) for x in src]
                kw = dict(xstar_left=2, xstar_right=2)
                kw["Ap_masks" if which == "Ap" else "Aq_masks"] = masks
                note = (src_name == "A1" and ((which == "Ap" and sig == SWAP_13) or (which == "Aq" and sig == SWAP_23)))
                r_, n_, v_ = run54("C4 %s=%s o %s%s" % (which, src_name, sig, " (note's)" if note else ""),
                                   True if note else None, **kw)
                tot += 1
                nval += v_
            log("   C4 summary: %s from %s: %d of %d well-formed permutations VALID" % (which, src_name, nval, tot))
    # C5/C6: planted collisions in the added products: put a B core (no inner port) into Ap / Aq, so a member of
    #        P x L(Ap) shares its right half with a glued member (and likewise on the left).
    star2 = ebit(0, 2) | ebit(1, 2) | ebit(2, 3) | ebit(2, 4)    # star at vertex 2: in F, avoids p=01
    star1 = ebit(0, 1) | ebit(1, 2) | ebit(1, 3) | ebit(1, 4)    # star at vertex 1: in F, avoids q=02
    Ap0 = [permute_mask(x, SWAP_13) for x in A1]
    Aq0 = [permute_mask(x, SWAP_23) for x in A1]
    run54("C5 Ap + star(2) planted", False, xstar_left=2, xstar_right=2, Ap_masks=Ap0 + [star2])
    run54("C6 Aq + star(1) planted", False, xstar_left=2, xstar_right=2, Aq_masks=Aq0 + [star1])
    # C7: X*-type set with every 3-edge core from the six non-port, non-pendant core edges (pendant 2 full)
    inner = [ebit(1, 2), ebit(1, 3), ebit(1, 4), ebit(2, 3), ebit(2, 4), ebit(3, 4)]
    good = []
    for comb in itertools.combinations(inner, 3):
        core = comb[0] | comb[1] | comb[2]
        box = (core, E, R.id([PEND_FULL]))
        sx = core | (PEND_FULL << 20)
        r_, n_, v_ = run54("C7 core %s + full pendant 2" % eset(core), None, extraB_left=[box], extraB_right=[box],
                           _xB_left=[sx], _xB_right=[sx])
        if v_:
            good.append(eset(core))
    log("   C7 summary: valid 3-edge cores with full pendant 2:", good)
    log("done")
