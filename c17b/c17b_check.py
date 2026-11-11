#!/usr/bin/env python3
"""Independent replay of R(3,33) >= 241 ([SK2026], cited on the 17b page by PR #252; Zenodo 10.5281/zenodo.23552014),
from the certificate data alone: the circulant graph G = Cay(Z_240, S) of R3_33_n240_Z240.json. The bundle's
checkers (check_ramsey.py, mcq.c, bbmis.c, drive.py, basic_checks.py) and its search code were neither read nor run.

The claim: G is triangle-free and has no independent set of 33 vertices, so colouring the edges of K_240 by G and its
complement has no blue K_3 and no red K_33, and R(3,33) > 240.

Checks (standard library; the exact independence search is our C program c17b_indep.c):
  1. the file: n = 240; the connection set S (both of its copies agree) is closed under negation, misses 0, |S| = 32;
     the edge list is exactly {x, x + s} (3,840 edges, no repeats); the group table is addition mod 240;
  2. triangle-free by two routes: (a) no s, t in S with s + t in S; (b) for every edge uv, N(u) and N(v) are disjoint;
     so N(0) is an independent set of 32 vertices and alpha(G) >= 32;
  3. the multipliers x -> a x, a in {13, 37, 71, 169, 203, 227, 239}, are units with a S = S (so they are
     automorphisms fixing 0), they and 1 form a group of order 8, and the file's stabilizer maps are exactly these;
     their orbits on R(0) (the 207 non-neighbours of 0) are computed;
  4. alpha(G) <= 32: c17b_indep.c with K = 33 must report NONE for every orbit (orbital branching, see its header);
  5. the second file (k = 34) holds the same graph, so it has no independent 34-set either.
Red controls (--controls): K = 32 must find an independent 32-set, which is then checked; S without the pair +-14
(no longer multiplier-invariant, so every vertex of R(0) is its own orbit) must give an independent 33-set, checked;
the C program against brute force on 60 random small circulants (both answers at K = alpha and alpha + 1).

Usage: python3 -I c17b_check.py CERTS_ROOT BINDIR [--controls] [--main] [--polite]
       (--main runs the K = 33 search, the heavy part, one orbit per process, resumable through results/orbits_K33.log;
        --polite waits while the machine is busy, see run_orbits; without --main only checks 1-3, 5 and the controls run)
"""
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from itertools import combinations
from math import gcd

MULT = [13, 37, 71, 169, 203, 227, 239]


def load(root, name):
    path = os.path.join(root, "sk2026", name)
    raw = open(path, "rb").read()
    print(f"file sk2026/{name} sha256 {hashlib.sha256(raw).hexdigest()}")
    return json.loads(raw)


def circulant_adj(n, S):
    Sset = set(S)
    return [set((u + s) % n for s in Sset) for u in range(n)]


def structural(d):
    n = d["n"]
    S = sorted(d["connection_set_Z240"])
    ok = True
    ok &= (n == 240)
    ok &= S == sorted(d["connection_set_indices"]) == sorted(x[0] for x in d["connection_set_labels"])
    ok &= len(S) == len(set(S)) == 32 and 0 not in S and all(0 < s < n for s in S)
    ok &= all((n - s) % n in S for s in S)
    ok &= all(d["element_labels"][i] == [i] for i in range(n))
    ok &= all(d["group_table"][i][j] == (i + j) % n for i in range(n) for j in range(n))
    edges = set()
    dup = 0
    for u, v in d["blue_edges"]:
        e = frozenset((u, v))
        dup += e in edges or u == v
        edges.add(e)
    want = {frozenset((x, (x + s) % n)) for x in range(n) for s in S}
    ok &= dup == 0 and edges == want and len(edges) == 3840
    return ok, n, S


def triangle_free(n, S):
    Sset = set(S)
    a = not any((s + t) % n in Sset for s in S for t in S)
    adj = circulant_adj(n, S)
    b = all(not (adj[u] & adj[v]) for u in range(n) for v in adj[u])
    return a, b, adj


def multipliers(n, S, d):
    Sset = set(S)
    units = all(gcd(a, n) == 1 for a in MULT)
    inv = all({(a * s) % n for s in S} == Sset for a in MULT)
    grp = set([1] + MULT)
    closed = all((a * b) % n in grp for a in grp for b in grp)
    maps = d.get("stabilizer_generators", [])
    file_maps = sorted(tuple(m) for m in maps)
    ours = sorted(tuple((a * x) % n for x in range(n)) for a in MULT)
    return units and inv, closed and len(grp) == 8, file_maps == ours, sorted(grp)


def orbits(n, R0, grp):
    seen, orbs = set(), []
    for v in sorted(R0):
        if v in seen:
            continue
        o = sorted({(a * v) % n for a in grp})
        assert set(o) <= R0
        seen |= set(o)
        orbs.append(o)
    return orbs


def run_indep(binary, n, S, K, orbs, label, cap=None):
    R0 = {v for o in orbs for v in o}
    orb = [-1] * n
    for i, o in enumerate(orbs):
        for v in o:
            orb[v] = i
    lines = [f"{n} {K}", f"{len(S)} " + " ".join(map(str, sorted(S))), str(len(orbs)),
             " ".join(map(str, orb)), " ".join(str(o[0]) for o in orbs)]
    t = time.time()
    p = subprocess.run([binary], input="\n".join(lines) + "\n", capture_output=True, text=True, timeout=cap)
    out = p.stdout
    found = None
    for line in out.splitlines():
        if "FOUND independent set" in line:
            found = [int(x) for x in line.split("{")[1].rstrip("}").split(",")]
    return p.returncode, out, found, time.time() - t


def load1():
    """the 1-minute load average (macOS sysctl), or None where it cannot be read"""
    try:
        out = subprocess.run(["sysctl", "-n", "vm.loadavg"], capture_output=True, text=True, timeout=10).stdout
        return float(out.replace("{", " ").split()[0])
    except Exception:
        return None


def run_orbits(binary, n, S, K, orbs, log_path, polite):
    """The same search as run_indep, one orbit per process (c17b_indep i 100 runs orbit i alone; each subproblem still
    excludes every orbit of smaller index), from the cheapest orbit (the last) to the first. Each finished orbit is
    appended to log_path with the SHA-256 of the input, so an interrupted run resumes where it stopped. With polite,
    an orbit starts only while the 1-minute load is below 36, and a running orbit is suspended (SIGSTOP) while the load
    is above 48 or a file PAUSE sits next to this script, and resumed (SIGCONT) when the load is below 36 again."""
    import hashlib as _h
    import signal
    orb = [-1] * n
    for i, o in enumerate(orbs):
        for v in o:
            orb[v] = i
    text = "\n".join([f"{n} {K}", f"{len(S)} " + " ".join(map(str, sorted(S))), str(len(orbs)),
                      " ".join(map(str, orb)), " ".join(str(o[0]) for o in orbs)]) + "\n"
    tag = _h.sha256(text.encode()).hexdigest()
    pause = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PAUSE")
    done = {}
    if os.path.exists(log_path):
        for line in open(log_path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 3 and parts[0] == tag:
                done[int(parts[1])] = parts[2]

    def busy():
        L = load1() if polite else None
        return (L is not None and L > 48) or os.path.exists(pause)

    def calm():
        L = load1() if polite else None
        return (L is None or L < 36) and not os.path.exists(pause)

    for i in reversed(range(len(orbs))):
        if i in done:
            continue
        while polite and not calm():
            time.sleep(60)
        p = subprocess.Popen([binary, str(i), "100"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        p.stdin.write(text)
        p.stdin.close()
        stopped, wait = False, 0.25
        while True:
            try:
                p.wait(timeout=wait)
                break
            except subprocess.TimeoutExpired:
                wait = min(20.0, 2 * wait)
            if polite and not stopped and busy():
                p.send_signal(signal.SIGSTOP)
                stopped = True
            elif stopped and calm():
                p.send_signal(signal.SIGCONT)
                stopped = False
        out = p.stdout.read()
        line = next((x for x in out.splitlines() if x.startswith("orbit")), None)
        if line is None or p.returncode not in (0, 1):
            raise RuntimeError(f"orbit {i}: the search did not report (exit {p.returncode})")
        done[i] = line
        with open(log_path, "a") as f:
            f.write(f"{tag}\t{i}\t{line}\n")
    lines = [done[i] for i in range(len(orbs))]
    nodes = sum(int(x.split("nodes ")[1].split()[0]) for x in lines)
    secs = sum(float(x.split("nodes ")[1].split()[1]) for x in lines)
    found = None
    for x in lines:
        if "FOUND independent set" in x:
            found = [int(y) for y in x.split("{")[1].rstrip("}").split(",")]
    res = (f"RESULT {'FOUND' if found else 'NONE'}: an independent set of size {K} "
           f"{'exists' if found else 'does not exist'} (nodes {nodes}, {secs:.1f} s)")
    return (1 if found else 0), "\n".join(lines + [res]), found, secs


def independent(adj, X):
    return all(v not in adj[u] for u, v in combinations(X, 2))


def alpha_brute(n, adj):
    best = 0
    order = list(range(n))

    def rec(i, chosen):
        nonlocal best
        if len(chosen) + (n - i) <= best:
            return
        if i == n:
            best = max(best, len(chosen))
            return
        v = order[i]
        if all(v not in adj[u] for u in chosen):
            chosen.append(v)
            rec(i + 1, chosen)
            chosen.pop()
        rec(i + 1, chosen)

    rec(0, [])
    return best


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    root, bindir = sys.argv[1], sys.argv[2]
    controls, do_main = "--controls" in sys.argv[3:], "--main" in sys.argv[3:]
    binary = os.path.join(bindir, "c17b_indep")
    t0 = time.time()
    status = 0
    d33, d34 = load(root, "R3_33_n240_Z240.json"), load(root, "R3_34_n240_Z240.json")
    ok1, n, S = structural(d33)
    print(f"1. structure: n = {n}, |S| = {len(S)}, S = -S, edge list = Cay(Z_240, S) (3,840 edges), group table = "
          f"addition mod 240: {ok1}")
    print(f"   S = {S}")
    ta, tb, adj = triangle_free(n, S)
    deg = {len(a) for a in adj}
    N0 = sorted(adj[0])
    print(f"2. triangle-free, route a (S + S misses S): {ta}; route b (every edge's endpoints share no neighbour): {tb}; "
          f"degrees {sorted(deg)}; N(0) independent: {independent(adj, N0)} (so alpha >= 32)")
    m_ok, g_ok, f_ok, grp = multipliers(n, S, d33)
    R0 = set(range(1, n)) - set(N0)
    orbs = orbits(n, R0, grp)
    sizes = sorted(len(o) for o in orbs)
    print(f"3. multipliers {MULT}: units with aS = S: {m_ok}; with 1 a group of order 8: {g_ok}; the file's 7 "
          f"stabilizer maps are exactly x -> ax: {f_ok}; |R(0)| = {len(R0)}, {len(orbs)} orbits, sizes {sizes}")
    ok5 = sorted(d34["connection_set_Z240"]) == S and d34["n"] == n and d34["k"] == 34
    print(f"5. the k = 34 file holds the same graph (n = 240, same S): {ok5}")
    if not (ok1 and ta and tb and deg == {32} and m_ok and g_ok and f_ok and ok5):
        status = 1
    if do_main:
        log_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "orbits_K33.log")
        rc, out, found, sec = run_orbits(binary, n, S, 33, orbs, log_path, "--polite" in sys.argv[3:])
        print("4. alpha(G) <= 32: independent 33-sets through 0, by orbit (c17b_indep, K = 33, one orbit per process)")
        print("   " + out.strip().replace("\n", "\n   "))
        main_ok = rc == 0 and found is None and "RESULT NONE" in out
        print(f"   alpha(G) = 32: {main_ok}  ({sec:.0f} s of search CPU)")
        if not main_ok:
            status = 1
    else:
        print("4. the K = 33 search was not run in this call (use --main)")
    if controls:
        print("RED CONTROLS")
        ok_all = True
        rc, out, found, sec = run_indep(binary, n, S, 32, orbs, "c1")
        r = rc == 1 and found is not None and len(set(found)) == 32 and independent(adj, found)
        print(f"  C1 K = 32 (must find): found {sorted(found) if found else None}, independent in G: {r} ({sec:.1f} s)")
        ok_all &= r
        S2 = [s for s in S if s not in (14, n - 14)]
        adj2 = circulant_adj(n, S2)
        R02 = set(range(1, n)) - adj2[0]
        single = [[v] for v in sorted(R02)]
        rc, out, found, sec = run_indep(binary, n, S2, 33, single, "c2", cap=3600)
        r = rc == 1 and found is not None and len(set(found)) == 33 and independent(adj2, found)
        print(f"  C2 S without +-14, K = 33 (must find): found {sorted(found) if found else None}, independent: {r} "
              f"({sec:.1f} s)")
        ok_all &= r
        rng = random.Random(20261110)
        agree, tested = True, 0
        for _ in range(60):
            m = rng.randint(9, 34)
            half = rng.sample(range(1, m // 2 + 1), rng.randint(1, max(1, m // 4)))
            S3 = sorted({s % m for h in half for s in (h, m - h)})
            adj3 = circulant_adj(m, S3)
            a = alpha_brute(m, adj3)
            R03 = set(range(1, m)) - adj3[0]
            single3 = [[v] for v in sorted(R03)]
            for K in (a, a + 1):
                if K < 3 or not single3:
                    continue
                rc, out, found, _ = run_indep(binary, m, S3, K, single3, "c3")
                want_found = K <= a
                ok = (found is not None) == want_found and (found is None or (len(set(found)) == K and independent(adj3, found)))
                agree &= ok
                tested += 1
        print(f"  C3 c17b_indep against brute force on 60 random circulants ({tested} decisions at K = alpha, alpha + 1): {agree}")
        ok_all &= agree
        print(f"  controls {'all as expected' if ok_all else 'FAILED'}")
        if not ok_all:
            status = 1
    print(f"RESULT {'R(3,33) >= 241 CONFIRMED' if (status == 0 and do_main) else ('checks 1-3, 5 hold (search not run)' if status == 0 else 'NOT CONFIRMED')}; "
          f"run time {time.time() - t0:.0f} s")
    return status


if __name__ == "__main__":
    sys.exit(main())
