#!/usr/bin/env python3
"""Program B for the 49a block of PR #370 (252,303 sets on 28 points): written for this folder, sharing no code with
../c49/sunflower_check.py or with c49p370_a.py.

Method. The family is stored as a quasi-reduced decision diagram over the points 0, 1, ..., n-1 in this order: a
node at level L is the set of tails (points L..n-1) of all members with one given head (points 0..L-1); equal sets of
tails get one node, so the diagram of a glued family is small. Then two joins over the diagram, each a memoised
search over tuples of nodes, one level (point) at a time:

  containment  pairs (A, B) with, at every point, (a, b) in {00, 11, 01} (so A <= B), and A != B somewhere;
  sunflower    triples (A, B, C) with, at every point, (a, b, c) in {000, 111, 100, 010, 001}, and the three
               members pairwise different (three bits remember which pairs already differ).

A triple of distinct sets is a 3-sunflower exactly when no point lies in exactly two of them (then A&B, A&C and B&C
are all equal to A&B&C), so the join answers the question for every triple of members at once. Every reported
witness is checked again against the definition. Then: port counts, the transfer matrix and its larger
eigenvalue as the root of x^2 - tr x + det, and the rate rho^(1/(n-2)), all by integer comparisons.

Standard library only; run with python3 -I. Usage: c49p370_b.py FILE.json [--control NAME]
"""
import json
import sys
import time
from math import isqrt

sys.setrecursionlimit(10000)


def load(path):
    with open(path, "rb") as fh:
        doc = json.loads(fh.read().decode("ascii"))
    n = doc["n"]
    ports = doc["ports"]
    masks = doc["masks"]
    assert isinstance(n, int) and 3 <= n <= 62
    assert all(isinstance(m, int) and 0 <= m < (1 << n) for m in masks)
    return n, ports, masks


class Diagram:
    """Quasi-reduced decision diagram: level L in 0..n, node ids per level; id 0 = no member (at every level),
    id 1 at level n = the empty tail."""

    def __init__(self, n, masks):
        self.n = n
        # bottom-up: at level L the key of a member is its head (bits 0..L-1); node of a head = (child0, child1)
        heads = {m: 1 for m in masks}  # level n: each full member is the terminal node 1
        self.kids = [None] * (n + 1)  # kids[L][id] = (id of child with point L absent, id with point L present)
        self.kids[n] = [None, None]
        for L in range(n - 1, -1, -1):
            pairs = {}
            for h, node in heads.items():
                low = h & ~(1 << L)
                bit = (h >> L) & 1
                pr = pairs.setdefault(low, [0, 0])
                assert pr[bit] == 0, "duplicate member"
                pr[bit] = node
            table = {}
            ids = [None]
            new = {}
            for low, pr in pairs.items():
                key = (pr[0], pr[1])
                if key not in table:
                    table[key] = len(ids)
                    ids.append(key)
                new[low] = table[key]
            self.kids[L] = ids
            heads = new
        assert list(heads) == [0] and len(self.kids[0]) == 2
        self.root = heads[0]

    def width(self):
        return [len(k) - 1 for k in self.kids]

    def count(self):
        memo = [dict() for _ in range(self.n + 1)]

        def c(L, v):
            if v == 0:
                return 0
            if L == self.n:
                return 1
            if v not in memo[L]:
                a, b = self.kids[L][v]
                memo[L][v] = c(L + 1, a) + c(L + 1, b)
            return memo[L][v]
        return c(0, self.root)


SUN_PATTERNS = [(0, 0, 0), (1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
SUB_PATTERNS = [(0, 0), (1, 1), (0, 1)]


def sunflower_join(D):
    """A 3-sunflower (A, B, C) of the family, as masks, or None. diff bits: 1 = A != B, 2 = A != C, 4 = B != C."""
    n, kids = D.n, D.kids
    dead = [set() for _ in range(n + 1)]
    states = [0]

    def walk(L, a, b, c, diff):
        if L == n:
            return 0 if diff == 7 else None
        key = (a, b, c, diff)
        if key in dead[L]:
            return None
        states[0] += 1
        ka, kb, kc = kids[L][a], kids[L][b], kids[L][c]
        for x, y, z in SUN_PATTERNS:
            na, nb, nc = ka[x], kb[y], kc[z]
            if na and nb and nc:
                nd = diff | (x != y) | ((x != z) << 1) | ((y != z) << 2)
                r = walk(L + 1, na, nb, nc, nd)
                if r is not None:
                    return r | (x << L) | (y << (n + L)) | (z << (2 * n + L))
        dead[L].add(key)
        return None

    r = walk(0, D.root, D.root, D.root, 0)
    full = (1 << n) - 1
    found = None if r is None else (r & full, (r >> n) & full, (r >> (2 * n)) & full)
    return found, states[0]


def containment_join(D):
    """(A, B) with A a proper subset of B, both members, or None."""
    n, kids = D.n, D.kids
    dead = [set() for _ in range(n + 1)]
    states = [0]

    def walk(L, a, b, diff):
        if L == n:
            return 0 if diff else None
        key = (a, b, diff)
        if key in dead[L]:
            return None
        states[0] += 1
        ka, kb = kids[L][a], kids[L][b]
        for x, y in SUB_PATTERNS:
            na, nb = ka[x], kb[y]
            if na and nb:
                r = walk(L + 1, na, nb, diff | (x != y))
                if r is not None:
                    return r | (x << L) | (y << (n + L))
        dead[L].add(key)
        return None

    r = walk(0, D.root, D.root, 0)
    full = (1 << n) - 1
    return (None if r is None else (r & full, (r >> n) & full)), states[0]


def larger_root_at_least(q_num, q_den, tr, det):
    """Is q = q_num/q_den <= rho, the larger root of x^2 - tr x + det (two real roots)? For q <= tr/2 yes; above,
    q <= rho iff the polynomial is <= 0 at q. Integers only."""
    if 2 * q_num <= tr * q_den:
        return True
    return q_num * q_num - tr * q_num * q_den + det * q_den * q_den <= 0


def rate_ok(num, k, L, tr, det):
    """(num / 10^k)^L <= rho, exactly."""
    return larger_root_at_least(num ** L, 10 ** (k * L), tr, det)


def selftest(trials=3000):
    """Both joins against the definitions (every pair and triple) on random small families, including glued ones
    (a product of two random families on disjoint point sets, plus a few extra sets), so that both answers occur."""
    import random
    rng = random.Random(370)
    seen = {"sun": [0, 0], "sub": [0, 0]}
    for t in range(trials):
        n = rng.randint(3, 9)
        if t % 3 == 0 and n >= 6:
            h = n // 2
            f1 = rng.sample(range(1 << h), rng.randint(1, 5))
            f2 = rng.sample(range(1 << (n - h)), rng.randint(1, 5))
            fam = sorted({a | (b << h) for a in f1 for b in f2} | set(rng.sample(range(1 << n), rng.randint(0, 2))))
        else:
            fam = rng.sample(range(1 << n), rng.randint(1, min(14, 1 << n)))
        D = Diagram(n, fam)
        assert D.count() == len(fam)
        sun, _ = sunflower_join(D)
        sub, _ = containment_join(D)
        brute_sun = any(a & b == a & c == b & c for i, a in enumerate(fam) for j, b in enumerate(fam[i + 1:], i + 1)
                        for c in fam[j + 1:])
        brute_sub = any(a != b and a & b == a for a in fam for b in fam)
        assert (sun is not None) == brute_sun, (n, fam)
        assert (sub is not None) == brute_sub, (n, fam)
        seen["sun"][brute_sun] += 1
        seen["sub"][brute_sub] += 1
    print("B: check self-test: both joins agree with the definitions on %d random small families "
          "(sunflower-free %d, with a sunflower %d; Sperner %d, with a containment %d): PASS"
          % (trials, seen["sun"][0], seen["sun"][1], seen["sub"][0], seen["sub"][1]))


def main():
    args = sys.argv[1:]
    if args[0] == "--selftest":
        selftest()
        print("B: VERDICT PASS")
        return
    path = args[0]
    control = args[args.index("--control") + 1] if "--control" in args else None
    t0 = time.time()
    n, ports, masks = load(path)
    printed = "1.5678643"
    expect = {"members": 252303, "N": (33980, 64728, 64728, 88867), "D": 3019700660}
    if control:
        print("B: CONTROL", control)
    if control == "sunflower":
        # a 3-sunflower planted: the two-point set {2, 3} with two points of the second pendant copy, and two more
        # sets on the same kernel with disjoint petals; the family is no longer sunflower-free whatever else holds
        masks = masks + [(1 << 2) | (1 << 5), (1 << 2) | (1 << 9), (1 << 2) | (1 << 20)]
    elif control == "sunflower2":
        # one new set C = (A & B) | P with P outside A | B, so that {A, B, C} is a 3-sunflower, chosen comparable
        # with no member (the family stays Sperner): only the sunflower test can refuse it
        import random
        rng = random.Random(49)
        S = set(masks)
        while True:
            a, b = rng.sample(masks, 2)
            free = [i for i in range(n) if not (a | b) >> i & 1]
            if len(free) < 2:
                continue
            c = (a & b)
            for i in rng.sample(free, rng.randint(1, len(free))):
                c |= 1 << i
            if c not in S and all(c & m not in (c, m) for m in masks):
                break
        print("B: planted %d with %d and %d" % (c, a, b))
        masks = masks + [c]
    elif control == "subset":
        # a member with one point removed (proper subset of a member)
        m = masks[1000]
        masks = masks + [m & (m - 1)]
    elif control == "drop":
        # one of the 32 added seven-point sets (64768, inside the first pendant) removed: still a valid family,
        # but N00 and the bound change
        assert 64768 in masks
        masks = [m for m in masks if m != 64768]
    elif control == "printed":
        printed = "1.5678644"
    ok = True

    def check(desc, good, detail=""):
        nonlocal ok
        ok = ok and good
        print("B: check %s: %s%s" % (desc, "PASS" if good else "FAIL", (" (%s)" % detail) if detail else ""))

    check("n = 28 points, ports (0, 1), %d members, all distinct" % len(masks),
          n == 28 and ports == [0, 1] and len(set(masks)) == len(masks), "%d distinct" % len(set(masks)))
    D = Diagram(n, masks)
    check("the diagram holds exactly the members", D.count() == len(set(masks)),
          "%d nodes, widest level %d" % (sum(D.width()), max(D.width())))
    t = time.time()
    sub, st2 = containment_join(D)
    if sub is not None:
        a, b = sub
        assert a in set(masks) and b in set(masks) and a != b and a & b == a
    check("Sperner (no member a proper subset of another), pair join over the diagram", sub is None,
          "%d states, %.1f s" % (st2, time.time() - t) if sub is None else "witness %s < %s" % sub)
    t = time.time()
    sun, st3 = sunflower_join(D)
    if sun is not None:
        a, b, c = sun
        S = set(masks)
        assert a in S and b in S and c in S and len({a, b, c}) == 3 and a & b == a & c == b & c
    check("3-sunflower-free, triple join over the diagram", sun is None,
          "%d states, %.1f s" % (st3, time.time() - t) if sun is None else "witness %s" % (sun,))
    p, q = ports
    N = [[0, 0], [0, 0]]
    for m in masks:
        N[(m >> p) & 1][(m >> q) & 1] += 1
    n00, n01, n10, n11 = N[0][0], N[0][1], N[1][0], N[1][1]
    print("B: port counts N00 N01 N10 N11 = %d %d %d %d" % (n00, n01, n10, n11))
    check("members and port counts as PR #370 states (252303; 33980, 64728, 64728, 88867)",
          (len(masks), (n00, n01, n10, n11)) == (expect["members"], expect["N"]))
    # transfer matrix of the chain (rows (N10, N11), (N00, N01)): trace and determinant
    tr = n10 + n01
    det = n10 * n01 - n11 * n00
    disc = tr * tr - 4 * det
    r = isqrt(disc)
    check("rho = 64728 + sqrt(3019700660): trace 2 * 64728 and discriminant 4 * 3019700660, not a square",
          tr == 2 * 64728 and disc == 4 * expect["D"] and r * r != disc and expect["D"] == n00 * n11)
    L = n - len(ports)
    # the rate to 22 decimals by bisection on integers
    k = 22
    lo, hi = 10 ** k, 2 * 10 ** k
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rate_ok(mid, k, L, tr, det):
            lo = mid
        else:
            hi = mid
    s = str(lo)
    print("B: rho^(1/%d) in [%s.%s, %s.%s + 1e-22)" % (L, s[0], s[1:], s[0], s[1:]))
    check("the stated 1.567864384572806979 are the first 18 decimals", s.startswith("1567864384572806979"))
    whole, frac = printed.split(".")
    good = rate_ok(int(whole + frac), len(frac), L, tr, det)
    check("printed lower bound %s <= rho^(1/26) (exact: %s^26 <= rho)" % (printed, printed), good)
    print("B: run time %.1f s" % (time.time() - t0))
    print("B: VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
