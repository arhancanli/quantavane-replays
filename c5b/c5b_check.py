#!/usr/bin/env python3
"""Independent replay of the exact small values of f(n) on the 5b page (Erdos problem #757), [SK2026] (PR #253),
from the thirteen (4,5)-sets of its data files and the page's definitions alone. The repository's verify.py and
check_all.sh were neither read nor run. Standard library only (Python 3.9+); run with python3 -I.

Definitions (5b page). A finite set A of reals is a (4,5)-set if every four of its elements determine at least five
distinct differences (of the six positive differences of a 4-set, at most one value repeats once). h(A) is the size
of a largest Sidon subset of A (all pairwise sums distinct, equivalently all positive differences distinct), and
f(n) = min h(A) over (4,5)-sets of size n. The page states f(n) = 4, 5, 5, 6, 6, 7, 7, 8, 8, 9 for 6 <= n <= 15 and
f(17) = 10, f(19) = 11, f(21) = 12: each data file is a (4,5)-set of n integers with h(A) = n - floor((8n - 6)/17),
and f(n) >= n - floor((8n - 6)/17) is the argument of Ma and Tang (arXiv:2602.23282, Lemma 2.4 and Theorem 5.3,
Henning-Yeo) with the -6 kept. That lower bound is mathematics and is not computed here.

What is computed, for each file:
  1. the file: n canonical integers (n from the file name), pairwise distinct;
  2. the (4,5) property by two routes that share no code: (a) every 4-subset, its six differences, at least five
     distinct; (b) every pair of pairs with equal difference, as an index set of size 3 or 4, attributed to the
     4-sets containing it; a 4-set fails exactly when it receives two or more (one value three times gives three);
  3. h(A) by two routes: (a) branch and bound over elements in order with the set of used differences (difference
     test); (b) brute force over all subsets of size h and h + 1 with the sum test (all pairwise sums distinct);
     route (b) must find a Sidon h-subset and no Sidon (h + 1)-subset, and route (a) must return h;
  4. h(A) against the formula n - floor((8n - 6)/17) and the page's table.
Red controls (--controls): the largest element of n21 replaced by a1 + a2 - a0 (a parallelogram: one 4-set with two
coincidences; both routes of 2 refuse); a Sidon set of 21 powers of two (h = 21, so the claim is refused); a file
with a repeated element (refused); every file with h claimed one larger (refused: no Sidon subset of that size); a
self-test on 200 random small sets (the two (4,5) routes agree, and so do the two h routes).

Usage: python3 -I c5b_check.py CERTS_ROOT [--controls]     (CERTS_ROOT holds sk2026/n6.txt ... n21.txt)
"""
import hashlib
import itertools
import os
import re
import sys
import time

TABLE = {6: 4, 7: 5, 8: 5, 9: 6, 10: 6, 11: 7, 12: 7, 13: 8, 14: 8, 15: 9, 17: 10, 19: 11, 21: 12}
INT = re.compile(r"-?(0|[1-9][0-9]*)\Z")


def read_set(path, n):
    toks = open(path, "rb").read().decode("ascii").split()
    if len(toks) != n:
        raise ValueError(f"{path}: {len(toks)} tokens, expected {n}")
    for t in toks:
        if not INT.match(t) or t == "-0":
            raise ValueError(f"{path}: token {t!r} is not a canonical integer")
    a = [int(t) for t in toks]
    if len(set(a)) != n:
        raise ValueError(f"{path}: repeated element")
    return sorted(a)


def is45_route_a(a):
    """Every 4-subset has at least five distinct positive differences."""
    for q in itertools.combinations(a, 4):
        d = {abs(x - y) for x, y in itertools.combinations(q, 2)}
        if len(d) < 5:
            return False, q
    return True, None


def is45_route_b(a):
    """Equal-difference pairs of pairs, attributed to the 4-sets containing them; a 4-set fails iff it gets >= 2."""
    n = len(a)
    by_diff = {}
    for i in range(n):
        for j in range(i + 1, n):
            by_diff.setdefault(a[j] - a[i], []).append((i, j))
    count = {}
    for pairs in by_diff.values():
        for p, q in itertools.combinations(pairs, 2):
            u = frozenset(p) | frozenset(q)
            if len(u) == 4:
                count[u] = count.get(u, 0) + 1
            else:  # a 3-term progression x, y, z with y - x = z - y: lies in every 4-set that adds one element
                for x in range(n):
                    if x not in u:
                        w = u | {x}
                        count[w] = count.get(w, 0) + 1
    bad = [k for k, v in count.items() if v >= 2]
    return (not bad), (tuple(sorted(a[i] for i in min(bad, key=sorted))) if bad else None)


def h_route_a(a):
    """Largest Sidon subset by branch and bound over the elements in order, with the set of used differences."""
    n = len(a)
    best = [0, ()]

    def rec(k, chosen, diffs):
        if len(chosen) + (n - k) <= best[0]:
            return
        if k == n:
            best[0], best[1] = len(chosen), tuple(chosen)
            return
        x = a[k]
        new = [x - y for y in chosen]
        if len(set(new)) == len(new) and not (diffs & set(new)):
            chosen.append(x)
            rec(k + 1, chosen, diffs | set(new))
            chosen.pop()
        rec(k + 1, chosen, diffs)

    rec(0, [], frozenset())
    return best[0], best[1]


def sidon_by_sums(s):
    sums = [x + y for x, y in itertools.combinations_with_replacement(s, 2)]
    return len(set(sums)) == len(sums)


def h_route_b(a, h):
    """A Sidon h-subset exists and no (h + 1)-subset is Sidon, by brute force with the sum test."""
    witness = next((s for s in itertools.combinations(a, h) if sidon_by_sums(s)), None)
    larger = next((s for s in itertools.combinations(a, h + 1) if sidon_by_sums(s)), None)
    return witness, larger


def check_one(a, n, stated=None, label=""):
    ok45a, wa = is45_route_a(a)
    ok45b, wb = is45_route_b(a)
    t = time.time()
    ha, wit_a = h_route_a(a)
    formula = n - (8 * n - 6) // 17
    want = formula if stated is None else stated
    wit_b, larger = h_route_b(a, want)
    hb_ok = wit_b is not None and larger is None
    lines = [f"{label} n = {n}: (4,5) route a {'holds' if ok45a else 'FAILS at ' + str(wa)}, "
             f"route b {'holds' if ok45b else 'FAILS at ' + str(wb)}",
             f"   h(A) route a (branch and bound, difference test) = {ha}, a largest Sidon subset {list(wit_a)}",
             f"   route b (all {want}- and {want + 1}-subsets, sum test): Sidon {want}-subset "
             f"{'found ' + str(list(wit_b)) if wit_b else 'NONE'}; Sidon {want + 1}-subset "
             f"{'FOUND ' + str(list(larger)) if larger else 'none'}",
             f"   n - floor((8n - 6)/17) = {formula}; page's f({n}) = {TABLE.get(n)}; {time.time() - t:.2f} s"]
    good = ok45a and ok45b and ha == want and hb_ok and want == formula == TABLE.get(n)
    return good, ok45a, ok45b, ha, hb_ok, lines


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    controls = "--controls" in sys.argv[2:]
    t0 = time.time()
    status = 0
    sets = {}
    for n in sorted(TABLE):
        path = os.path.join(root, "sk2026", f"n{n}.txt")
        print(f"file sk2026/n{n}.txt sha256 {hashlib.sha256(open(path, 'rb').read()).hexdigest()}")
        a = read_set(path, n)
        sets[n] = a
        good, *_, lines = check_one(a, n, label="[SK2026]")
        print("\n".join(lines))
        print(f"   {'CONFIRMED' if good else 'NOT CONFIRMED'}: f({n}) <= {TABLE[n]} by this set, equal to the lower bound")
        if not good:
            status = 1
    # the table against the formula, for every listed n
    ftab = all(TABLE[n] == n - (8 * n - 6) // 17 for n in TABLE)
    print(f"formula n - floor((8n-6)/17) equals the page's table at all 13 sizes: {ftab}")
    if not ftab:
        status = 1
    if controls:
        print("RED CONTROLS (each must be refused)")
        ok_all = True
        # R1: the largest element of the n = 21 set replaced by v = a1 + a2 - a0, so {a0, a1, a2, v} is a parallelogram
        # (v - a2 = a1 - a0 and v - a1 = a2 - a0): two coincidences in one 4-set
        a = list(sets[21])
        v = a[1] + a[2] - a[0]
        b = sorted(a[:-1] + [v])
        r = (v not in a) and (not is45_route_a(b)[0]) and (not is45_route_b(b)[0])
        print(f"  R1 the largest element of n21 replaced by {v} = a1 + a2 - a0: both (4,5) routes refuse: {r}")
        ok_all &= r
        # R2: a Sidon set of size 21 (powers of two): h = 21 > 12, so the claim f(21) = 12 by this set is refused
        s = [2 ** i for i in range(21)]
        good, ok45a, ok45b, ha, hb_ok, _ = check_one(s, 21, stated=12, label="control")
        r = (not good) and ha == 21 and not hb_ok
        print(f"  R2 a Sidon set of 21 powers of two: h = {ha}, claim refused: {r}")
        ok_all &= r
        # R3: a repeated element is refused by the reader
        path = os.path.join(root, "sk2026", "n21.txt")
        toks = open(path).read().split()
        try:
            tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "control_repeat.txt")
            os.makedirs(os.path.dirname(tmp), exist_ok=True)
            open(tmp, "w").write(" ".join(toks[:-1] + [toks[0]]) + "\n")
            read_set(tmp, 21)
            r = False
        except ValueError:
            r = True
        finally:
            try:
                os.remove(tmp)
            except OSError:
                pass
        print(f"  R3 n21 with its last element replaced by its first (a repeat): refused: {r}")
        ok_all &= r
        # R4: every file with its stated value raised by one: no Sidon subset of that size, so the claim is refused
        r = True
        for n, a in sets.items():
            good = check_one(a, n, stated=TABLE[n] + 1, label="control")[0]
            r &= (not good) and h_route_b(a, TABLE[n] + 1)[0] is None
        print(f"  R4 each of the 13 files with h claimed one larger: refused, no Sidon subset of that size: {r}")
        ok_all &= r
        # R5: the two Sidon tests and the two h routes agree on 200 random small sets (self-test)
        import random
        rng = random.Random(20261110)
        agree = True
        for _ in range(200):
            m = rng.randint(4, 11)
            b = sorted(rng.sample(range(0, 60), m))
            ha, wa = h_route_a(b)
            wb, lb = h_route_b(b, ha)
            agree &= (wb is not None and lb is None and sidon_by_sums(wa))
            agree &= (is45_route_a(b)[0] == is45_route_b(b)[0])
        print(f"  R5 self-test, 200 random sets: the two (4,5) routes and the two h routes agree: {agree}")
        ok_all &= agree
        print(f"  controls {'all as expected' if ok_all else 'FAILED'}")
        if not ok_all:
            status = 1
    print(f"RESULT {'every value CONFIRMED' if status == 0 else 'NOT CONFIRMED'}; run time {time.time() - t0:.1f} s")
    return status


if __name__ == "__main__":
    sys.exit(main())
