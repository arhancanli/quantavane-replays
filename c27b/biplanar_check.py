#!/usr/bin/env python3
"""Exact replay of the two 27b lower-bound rows (teorth/optimizationproblems, main 87b4ebb; both added in PR #23):

  [R1959]  C_27b >= 8: K_8 is biplanar (Ringel), and chi(K_8) = 8.
  [G1980]  C_27b >= 9: "the join of a 6-vertex complete graph and a 5-vertex cycle graph" (Sulanke, in Gardner's
           column) is biplanar, and its chromatic number is 6 + 3 = 9.

C_27b is the largest chromatic number of a biplanar graph (thickness at most 2: the edge set is the union of the edge
sets of two planar graphs).  For each graph G we
  1. search for a split E(G) = E_1 + E_2 with both parts planar (randomised greedy, fixed seeds, with a planarity
     routine of Demoucron-Malgrange-Pertuiset type that also returns the faces of an embedding);
  2. verify the result with code that shares nothing with the search: each part is connected and spanning, every edge
     of the part appears exactly once in each direction in the face list and no other dart does, the faces define a
     rotation system (at every vertex v the successor map on its neighbours, read off the faces, is one cycle of
     length deg(v)), and V - E + F = 2, so the embedding is a sphere and the part is planar; E_1 and E_2 are disjoint
     and their union is E(G);
  3. compute chi(G): no proper colouring with chi - 1 colours (exhaustive backtracking with the first clique's
     colours fixed, which loses nothing by symmetry) and a proper colouring with chi colours.
The faces are printed, so the embeddings can be checked by hand.  No download.
Red controls (--controls): a part with one more edge than 3n - 6 is refused by the Euler bound and by the planarity
routine; a face list with two vertices swapped in one face is refused by the verifier; K_6 + C_5 has no 8-colouring
(a claimed 8 is refused); K_5 as a part is refused (10 > 3*5 - 6 edges).
usage: python3 -I biplanar_check.py [--controls]
"""
import itertools
import random
import sys

CONTROLS = "--controls" in sys.argv


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


# ------------------------------------------------------------------ graphs
def k_n(n):
    return list(range(n)), [frozenset(e) for e in itertools.combinations(range(n), 2)]


def k6_join_c5():
    V = list(range(11))
    E = [frozenset(e) for e in itertools.combinations(range(6), 2)]
    E += [frozenset((6 + i, 6 + (i + 1) % 5)) for i in range(5)]
    E += [frozenset((a, b)) for a in range(6) for b in range(6, 11)]
    return V, E


# ------------------------------------------------------------------ search side: DMP planarity with faces
def adjacency(V, E):
    adj = {v: set() for v in V}
    for e in E:
        a, b = tuple(e)
        adj[a].add(b)
        adj[b].add(a)
    return adj


def find_cycle(adj):
    start = next(v for v in adj if adj[v])
    parent = {start: None}
    stack = [(start, iter(sorted(adj[start])))]
    while stack:
        v, it = stack[-1]
        for w in it:
            if w == parent[v]:
                continue
            if w in parent:
                cyc = [v]
                u = v
                while u != w:
                    u = parent[u]
                    cyc.append(u)
                    if u is None:
                        break
                if cyc[-1] == w:
                    return cyc
                continue
            parent[w] = v
            stack.append((w, iter(sorted(adj[w]))))
            break
        else:
            stack.pop()
    return None


def dmp_faces(V, E):
    """Faces of a planar embedding of the 2-connected graph (V, E), or None if it is not planar."""
    adj = adjacency(V, E)
    cyc = find_cycle(adj)
    if cyc is None:
        return None
    Hv = set(cyc)
    He = {frozenset((cyc[i], cyc[(i + 1) % len(cyc)])) for i in range(len(cyc))}
    faces = [list(cyc), list(reversed(cyc))]
    Eset = set(E)
    while len(He) < len(Eset):
        frags = []
        for e in Eset - He:
            a, b = tuple(e)
            if a in Hv and b in Hv:
                frags.append(("edge", {a, b}, [a, b]))
        seen = set()
        for v in V:
            if v in Hv or v in seen:
                continue
            comp, stack = {v}, [v]
            while stack:
                u = stack.pop()
                for w in adj[u]:
                    if w not in Hv and w not in comp:
                        comp.add(w)
                        stack.append(w)
            seen |= comp
            att = {w for u in comp for w in adj[u] if w in Hv}
            frags.append(("comp", att, comp))
        choice = None
        for f in frags:
            adm = [i for i, fc in enumerate(faces) if f[1] <= set(fc)]
            if not adm:
                return None
            if choice is None or len(adm) == 1:
                choice = (f, adm[0])
                if len(adm) == 1:
                    break
        (kind, att, body), fi = choice
        if kind == "edge":
            path = body
        else:
            a = sorted(att)[0]
            # a path from a into the component and out to another attachment vertex
            comp = body
            start = [w for w in adj[a] if w in comp][0]
            prev = {start: a}
            queue = [start]
            end = None
            while queue and end is None:
                u = queue.pop(0)
                for w in sorted(adj[u]):
                    if w in Hv and w != a:
                        end = (u, w)
                        break
                    if w in comp and w not in prev:
                        prev[w] = u
                        queue.append(w)
            u, b = end
            path = [b, u]
            while path[-1] != a:
                path.append(prev[path[-1]])
            path.reverse()
        a, b = path[0], path[-1]
        f = faces[fi]
        i, j = f.index(a), f.index(b)
        rot = f[i:] + f[:i]
        j = rot.index(b)
        inner = path[1:-1]
        f1 = rot[:j + 1] + list(reversed(inner))
        f2 = rot[j:] + [rot[0]] + inner
        faces[fi] = f1
        faces.append(f2)
        Hv |= set(path)
        He |= {frozenset((path[k], path[k + 1])) for k in range(len(path) - 1)}
    return faces


def two_connected(V, E):
    adj = adjacency(V, E)
    for x in V:
        rest = [v for v in V if v != x]
        seen, stack = {rest[0]}, [rest[0]]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w != x and w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) != len(rest):
            return False
    return True


def search_split(V, E, seeds=range(1, 2001)):
    for seed in seeds:
        rng = random.Random(seed)
        order = list(E)
        rng.shuffle(order)
        E1, E2 = [], []
        for e in order:
            trial = E1 + [e]
            sub_v = sorted({v for f in trial for v in f})
            if len(trial) <= 3 * len(V) - 6 and (len(trial) < 3 or not two_connected(sub_v, trial)
                                                  or dmp_faces(sub_v, trial) is not None):
                E1 = trial
            else:
                E2.append(e)
        if len(E1) + len(E2) != len(E):
            continue
        if not (two_connected(V, E1) and two_connected(V, E2)):
            continue
        f1, f2 = dmp_faces(V, E1), dmp_faces(V, E2)
        if f1 is not None and f2 is not None:
            return seed, (E1, f1), (E2, f2)
    return None


# ------------------------------------------------------------------ verification side (no code shared with the search)
def verify_planar(V, E, faces):
    """True iff the faces form an orientable 2-cell embedding of the connected spanning graph (V, E) with
    V - E + F = 2."""
    Eset = {tuple(sorted(e)) for e in E}
    darts = {}
    for fc in faces:
        if len(fc) < 3:
            return False, "a face with fewer than 3 vertices"
        for k in range(len(fc)):
            u, v = fc[k], fc[(k + 1) % len(fc)]
            if tuple(sorted((u, v))) not in Eset:
                return False, f"dart {u}->{v} is not an edge"
            if (u, v) in darts:
                return False, f"dart {u}->{v} twice"
            darts[(u, v)] = (fc, k)
    if len(darts) != 2 * len(Eset):
        return False, "not every dart appears"
    # rotation: in the face, dart (u -> v) is followed by (v -> w); successor at v maps u to w
    succ = {}
    for (u, v), (fc, k) in darts.items():
        w = fc[(k + 2) % len(fc)]
        succ[(v, u)] = w
    nbrs = {v: set() for v in V}
    for a, b in Eset:
        nbrs[a].add(b)
        nbrs[b].add(a)
    for v in V:
        if not nbrs[v]:
            return False, f"vertex {v} isolated"
        start = next(iter(nbrs[v]))
        cyc, u = [start], succ[(v, start)]
        while u != start:
            cyc.append(u)
            u = succ[(v, u)]
            if len(cyc) > len(nbrs[v]):
                break
        if len(cyc) != len(nbrs[v]) or set(cyc) != nbrs[v]:
            return False, f"the rotation at {v} is not one cycle"
    seen, stack = {V[0]}, [V[0]]
    while stack:
        x = stack.pop()
        for y in nbrs[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    if len(seen) != len(V):
        return False, "not connected"
    chi = len(V) - len(Eset) + len(faces)
    return chi == 2, f"V - E + F = {len(V)} - {len(Eset)} + {len(faces)} = {chi}"


def colourable(V, E, k):
    adj = adjacency(V, E)
    order = sorted(V, key=lambda v: -len(adj[v]))
    col = {}
    # fix the colours of a largest clique among the first vertices (symmetry): greedy clique on the order
    clique = []
    for v in order:
        if all(v in adj[u] for u in clique):
            clique.append(v)
    if len(clique) > k:
        return None
    for i, v in enumerate(clique):
        col[v] = i
    rest = [v for v in order if v not in col]

    def rec(i):
        if i == len(rest):
            return True
        v = rest[i]
        used = {col[w] for w in adj[v] if w in col}
        for c in range(k):
            if c not in used:
                col[v] = c
                if rec(i + 1):
                    return True
                del col[v]
        return False

    return dict(col) if rec(0) else None


def report(name, V, E, chi, tag, rows):
    found = search_split(V, E)
    need(found is not None, f"{name}: no split found")
    seed, (E1, f1), (E2, f2) = found
    ok1, why1 = verify_planar(V, E1, f1)
    ok2, why2 = verify_planar(V, E2, f2)
    union = {tuple(sorted(e)) for e in E1} | {tuple(sorted(e)) for e in E2}
    disjoint = not ({tuple(sorted(e)) for e in E1} & {tuple(sorted(e)) for e in E2})
    cover = union == {tuple(sorted(e)) for e in E}
    lower = colourable(V, E, chi - 1)
    upper = colourable(V, E, chi)
    good = ok1 and ok2 and disjoint and cover and lower is None and upper is not None
    print(f"{name}: {len(V)} vertices, {len(E)} edges; split found with seed {seed}: |E_1| = {len(E1)}, "
          f"|E_2| = {len(E2)} (each at most 3n - 6 = {3 * len(V) - 6})")
    for lab, Ei, fi, oki, whyi in (("E_1", E1, f1, ok1, why1), ("E_2", E2, f2, ok2, why2)):
        print(f"  {lab} = {sorted(tuple(sorted(e)) for e in Ei)}")
        print(f"  faces of {lab}: {fi}")
        print(f"  {lab} verified planar: {oki} ({whyi})")
    print(f"  E_1 and E_2 disjoint: {disjoint}; union = E(G): {cover}")
    print(f"  proper {chi - 1}-colouring: {'NONE' if lower is None else lower}; proper {chi}-colouring: "
          f"{upper if upper else 'NONE'}")
    print(f"RESULT [{tag}]: {name} is biplanar with chromatic number {chi}, so C_27b >= {chi}: "
          f"{'HOLDS' if good else 'FAILS'}")
    rows.append((E1, f1, E2, f2))
    return good


def main():
    ok = True
    rows = []
    V8, E8 = k_n(8)
    ok &= report("K_8", V8, E8, 8, "R1959", rows)
    V11, E11 = k6_join_c5()
    ok &= report("K_6 + C_5 (join)", V11, E11, 9, "G1980", rows)
    if CONTROLS:
        print("red controls (each must be refused):")
        res = []
        E1, f1, E2, f2 = rows[1]
        moved = E1 + [E2[0]]
        sub_ok = len(moved) <= 3 * len(V11) - 6
        res.append(("one more edge in E_1", not sub_ok or dmp_faces(V11, moved) is None))
        print(f"  E_1 plus one edge of E_2: {len(moved)} edges > 3n - 6 = {3 * len(V11) - 6}: {not sub_ok}")
        bad = [list(fc) for fc in f1]
        bad[0][0], bad[0][1] = bad[0][1], bad[0][0]
        okb, whyb = verify_planar(V11, E1, bad)
        print(f"  two vertices swapped in one face of E_1: {whyb}")
        res.append(("swapped face", not okb))
        res.append(("a claimed 8-colouring of K_6 + C_5", colourable(V11, E11, 8) is None))
        V5, E5 = k_n(5)
        res.append(("K_5 as a part", len(E5) > 3 * 5 - 6 and dmp_faces(V5, E5) is None))
        for label, r in res:
            print(f"  {label}: {'refused, as expected' if r else 'NOT AS EXPECTED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
