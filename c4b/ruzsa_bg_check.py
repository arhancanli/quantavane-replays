#!/usr/bin/env python3
"""Exact replay of the two base-expansion rows of the 4b page (teorth/optimizationproblems, main 87b4ebb):

  [Ruz1984]  C_4b >= (1/2)(1 + log 7/log 65)  ~ 0.733077  (a 7-element residue set modulo 65)
  [BG2008]   C_4b >= (1/2)(1 + log 12/log 205) ~ 0.733412  (modulus 205, the set S printed in the row)

The finite hypothesis of the base-expansion lemma, as the page describes it (the [Jon2026a] row: a square-free
modulus m, and any nonzero z^2 mod m counts as a square): m is square-free and no difference s - s' of two distinct
elements of S is a square modulo m.  For 205 the set is read from the page; for 65 the page prints no set, so we
search: an exhaustive clique search over Z_65 (translations fixed by 0 in S) gives the largest such sets, of size 7,
and one witness.  The exponents are compared with the printed decimals (written with an approximately-equal sign, so
rounded to six decimals) in 60-digit decimal arithmetic (logarithms correctly rounded).  The lemma itself (Ruzsa 1984; Beigel and Gasarch, Theorem 3.7) is read, not checked.
usage: python3 -I ruzsa_bg_check.py CERTS_ROOT [--controls]   (CERTS_ROOT holds page/4b.md, from ./fetch.sh)
"""
import re
import sys
from decimal import Decimal, getcontext

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv
getcontext().prec = 60


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def squarefree(m):
    return all(m % (p * p) for p in range(2, int(m ** 0.5) + 2))


def squares(m):
    return {(x * x) % m for x in range(m)} - {0}


def conflicts(S, m):
    Q = squares(m)
    return [(a, b) for a in S for b in S if a != b and (a - b) % m in Q]


def largest(m):
    """All maximum sets containing 0 with no nonzero square difference mod m (exhaustive branch and bound)."""
    Q = squares(m)
    ok = [[a != b and (a - b) % m not in Q and (b - a) % m not in Q for b in range(m)] for a in range(m)]
    best = [[0]]

    def ext(cur, cand):
        if len(cur) + len(cand) < len(best[0]):
            return
        if len(cur) > len(best[0]):
            best[:] = [list(cur)]
        elif len(cur) == len(best[0]) and cur != best[0]:
            best.append(list(cur))
        for i, c in enumerate(cand):
            ext(cur + [c], [d for d in cand[i + 1:] if ok[c][d]])

    ext([0], [d for d in range(1, m) if ok[0][d]])
    return best


def exponent(k, m):
    return (1 + Decimal(k).ln() / Decimal(m).ln()) / 2


def trunc(x, digits):
    q = int(x * 10 ** digits)
    return f"{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


def rnd(x, digits):
    q = int((x * 10 ** digits).to_integral_value())
    return f"{q // 10 ** digits}.{q % 10 ** digits:0{digits}d}"


def page_set(path):
    text = open(path, encoding="utf-8").read()
    m = re.search(r"modulus \$(\d+)\$ and residue set \$S=\\\{([0-9,]+)\\\}\$", text)
    need(m is not None, "the [BG2008] set was not found on the page")
    return int(m.group(1)), [int(x) for x in m.group(2).split(",")]


def main():
    ok = True
    # [Ruz1984]
    sets = largest(65)
    w = sets[0]
    good = squarefree(65) and len(w) == 7 and not conflicts(w, 65) and rnd(exponent(7, 65), 6) == "0.733077"
    ok &= good
    print(f"[Ruz1984]: m = 65 = 5 * 13 is square-free; the largest residue sets with no nonzero square difference "
          f"(containing 0) have {len(w)} elements ({len(sets)} of them), e.g. {w}")
    print(f"RESULT [Ruz1984]: C_4b >= (1 + log 7/log 65)/2 = {trunc(exponent(7, 65), 15)}...; the printed ~0.733077 is "
          f"its value to six decimals: {'HOLDS' if good else 'FAILS'}")
    # [BG2008]
    m, S = page_set(f"{ROOT}/page/4b.md")
    c = conflicts(S, m)
    good = m == 205 and squarefree(m) and len(set(S)) == 12 and all(0 <= s < m for s in S) and not c and \
        rnd(exponent(12, 205), 6) == "0.733412"
    ok &= good
    print(f"[BG2008]: m = {m} = 5 * 41 is square-free: {squarefree(m)}; S = {S} ({len(set(S))} distinct residues); "
          f"pairs whose difference is a nonzero square mod {m}: {len(c)}")
    print(f"RESULT [BG2008]: C_4b >= (1 + log 12/log 205)/2 = {trunc(exponent(12, 205), 15)}...; the printed ~0.733412 "
          f"is its value to six decimals (rounded up; the row gives the exact expression): {'HOLDS' if good else 'FAILS'}")
    if CONTROLS:
        print("red controls (each must fail):")
        res = []
        S2 = [x for x in S if x != 181] + [180]
        res.append(("[BG2008] with 180 for 181", bool(conflicts(S2, m))))
        print(f"  S with 180 for 181: {len(conflicts(S2, m))} ordered pairs with a square difference")
        res.append(("an 8-element set modulo 65", len(w) < 8))
        res.append(("[BG2008] printed 0.733413", rnd(exponent(12, 205), 6) != "0.733413"))
        res.append(("modulus 65 * 5 = 325 is not square-free", not squarefree(325)))
        for label, r in res:
            print(f"  {label}: {'fails, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
