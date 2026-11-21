#!/usr/bin/env python3
"""Validate the product/pattern reduction in c49_referee.decide against brute force on random small
instances (two halves, several products, boxes with core singletons and pendant families)."""
import random
import itertools
import sys
from c49_referee import decide, PendTables, explicit_half


def brute(products, fams):
    M = set()
    for (Ub, Vb) in products:
        Ls = set()
        for (c, a, b) in Ub:
            for x in fams[a]:
                for y in fams[b]:
                    Ls.add(explicit_half(c, x, y))
        Rs = set()
        for (c, a, b) in Vb:
            for x in fams[a]:
                for y in fams[b]:
                    Rs.add(explicit_half(c, x, y))
        for l in Ls:
            for r in Rs:
                M.add(l | (r << 30))
    M = sorted(M)
    sf = False
    for A, B, C in itertools.combinations(M, 3):
        ab = A & B
        if ab == (A & C) and ab == (B & C):
            sf = True
            break
    sp = False
    for A in M:
        for C in M:
            if A != C and (A & C) == A:
                sp = True
                break
        if sp:
            break
    return set(M), sf, sp


def main(trials, seed):
    rng = random.Random(seed)
    agree = 0
    stats = {}
    for t in range(trials):
        mode = t % 2
        w = rng.choice([2, 3]) if mode == 0 else 4
        nf = rng.randint(2, 5)
        fams = []
        if mode == 0:
            for _ in range(nf):
                k = rng.randint(1, 3)
                fams.append(tuple(sorted(set(rng.randrange(1 << w) for _ in range(k)))))
        else:
            # equal-size members everywhere -> Sperner holds; sunflower detection is the real test
            ksz = rng.randint(1, 3)
            pool = [m for m in range(1 << w) if bin(m).count("1") == ksz]
            for _ in range(nf):
                fams.append(tuple(sorted(set(rng.choice(pool) for _ in range(rng.randint(1, 3))))))
            cpool = [m for m in range(1 << w) if bin(m).count("1") == 2]
        T = PendTables(fams)
        nprod = rng.randint(1, 3)
        products = []
        pick_core = (lambda: rng.randrange(1 << w)) if mode == 0 else (lambda: rng.choice(cpool))
        for _ in range(nprod):
            Ub = [(pick_core(), rng.randrange(nf), rng.randrange(nf)) for _ in range(rng.randint(1, 3))]
            Vb = [(pick_core(), rng.randrange(nf), rng.randrange(nf)) for _ in range(rng.randint(1, 3))]
            products.append((Ub, Vb))
        M, sf, sp = brute(products, fams)
        res = decide(products, T)
        msf = bool(res["sunflower_product_triples"])
        msp = bool(res["sperner_product_pairs"])
        if msf != sf or msp != sp:
            print("MISMATCH trial", t, "brute", sf, sp, "method", msf, msp, products, fams)
            sys.exit(1)
        if msf:
            A, B, C = res["sf_witness"][1]
            assert A in M and B in M and C in M and len({A, B, C}) == 3 and (A & B) == (A & C) == (B & C)
        if msp:
            A, C = res["sp_witness"][1]
            assert A in M and C in M and A != C and (A & C) == A
        stats[(sf, sp)] = stats.get((sf, sp), 0) + 1
        agree += 1
    print("selftest: %d/%d agree; (sunflower, sperner-violation) outcome counts: %s" % (agree, trials, stats))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3000, int(sys.argv[2]) if len(sys.argv) > 2 else 1)
