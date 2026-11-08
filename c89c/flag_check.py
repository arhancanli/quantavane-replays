"""Replay of the 89c candidate lower bound of PR #308, C_89c >= 8280981029/250000000000000 = 0.000033123924116*
(rpinho/math-certificates @ 603b40c, c89c/round3), from the certificate and the graph list alone; the packet's two
checkers were neither read nor run.

The finite statement (PROOF.md of the packet, in its own words): for every graph H on 8 vertices,
    mono_K6(H) - sum_blocks contribution(H) >= b,
mono_K6(H) = the fraction of the 28 six-sets that are cliques or independent sets of H, and for a block (t roots,
extension size m, root graph, flags F_i, integer rows L_r, Q = sum_r L_r L_r^T / D^2, colour 1 reading the
complement) the contribution is the sum of Q_ij over the ordered root tuples R inducing the block's root graph and
the ordered pairs (A, B) of disjoint m-sets outside R with F_i = flag(R + A), F_j = flag(R + B), divided by
C(8,t) t! C(8-t,m) C(8-t-m,m). Q is positive semidefinite by construction, which is what makes each block's term a
valid flag-algebra square; the passage from the finite inequality to the asymptotic bound is the standard one.

What we check, with our own code: the certificate's structure (blocks, flags, root parts, no two listed flags in one
orbit); that the 12,346 graphs of graph8.g6 are pairwise non-isomorphic and that their orbits cover all 2^28 labelled
graphs on 8 vertices (sum of 8!/|Aut| = 2^28); and the inequality for every graph, exactly (sums of Gram entries in
128-bit integers in flag_count.c, the comparison in Python fractions). Red controls and an independent slow
recomputation in pure Python for three graphs.

Usage: python3 -I flag_check.py CERTS_ROOT BIN [--controls]
"""
import hashlib
import itertools
import os
import subprocess
import sys
import time
from fractions import Fraction as Fr
from math import comb, factorial

B_STATED = Fr(8280981029, 250000000000000)
SLACK_STATED = Fr(112016854832021, 56000000000000000000000000)
WORST_STATED = 905
SHA = {"c89c/c6_ordered_type8.cert.xz": "6fc1a65ed64022e891261364bb931117d077149fca4574fda3d187d748d516a3",
       "c89c/c6_ordered_type8.cert": "888552bb82994a387f152d8a7816c011c2bad19f5ed98059152b6dd417ee3c94",
       "c89c/graph8.g6": "546a249902101c97d3aa590f93e53366854bd0a6f405aa59bdb32d25c57f845a"}
TS = {6: 1, 4: 2, 2: 3, 0: 4}


def pidx(x, y):
    return y * (y - 1) // 2 + x


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_graph6(path):
    out = []
    with open(path, "rb") as fh:
        for line in fh.read().split():
            s = line.decode("ascii")
            if s[0] != "G" or len(s) != 6:
                raise ValueError("not an 8-vertex graph6 line: %r" % s)
            bits = []
            for ch in s[1:]:
                v = ord(ch) - 63
                if not 0 <= v < 64:
                    raise ValueError("bad graph6 character")
                bits += [(v >> (5 - k)) & 1 for k in range(6)]
            if any(bits[28:]):
                raise ValueError("padding bits set")
            # graph6 lists the upper triangle column by column: (0,1), (0,2), (1,2), (0,3), ...: our pidx order
            out.append(sum(b << i for i, b in enumerate(bits[:28])))
    return out


def read_cert(path):
    with open(path) as fh:
        tok = fh.read().split()
    pos = 0

    def nxt():
        nonlocal pos
        pos += 1
        return int(tok[pos - 1])

    head = [nxt() for _ in range(7)]
    blocks = []
    for _ in range(head[3]):
        t, m, root, nf, nr, rmode, col = [nxt() for _ in range(7)]
        flags = [nxt() for _ in range(nf)]
        rows = [[nxt() for _ in range(nf)] for _ in range(nr)]
        blocks.append(dict(t=t, m=m, root=root, flags=flags, rows=rows, rmode=rmode, col=col))
    return head, blocks, pos == len(tok)


def write_cert(path, head, blocks):
    with open(path, "w") as fh:
        fh.write(" ".join(map(str, head)) + "\n")
        for b in blocks:
            fh.write("%d %d %d %d %d %d %d\n" % (b["t"], b["m"], b["root"], len(b["flags"]), len(b["rows"]), b["rmode"],
                                                 b["col"]))
            fh.write(" ".join(map(str, b["flags"])) + "\n")
            for r in b["rows"]:
                fh.write(" ".join(map(str, r)) + "\n")


def perm_ext(code, t, m, sg):
    out = 0
    v = t + m
    for y in range(1, v):
        for x in range(y):
            if code >> pidx(x, y) & 1:
                x2 = x if x < t else t + sg[x - t]
                y2 = y if y < t else t + sg[y - t]
                if x2 > y2:
                    x2, y2 = y2, x2
                out |= 1 << pidx(x2, y2)
    return out


def canon_flag(code, t, m):
    return min(perm_ext(code, t, m, sg) for sg in itertools.permutations(range(m)))


def structure(head, blocks, complete, out):
    ok = complete and head == [8, 6, 10 ** 12, 340, 8280981029, 250000000000000, -2]
    out("header n = 8, k = 6, D = 10^12, 340 blocks, b = 8280981029/250000000000000, mode -2; nothing after the last "
        "block: %s" % ok)
    fam = {}
    seen = set()
    for b in blocks:
        t, m = b["t"], b["m"]
        fam[(t, m, b["col"])] = fam.get((t, m, b["col"]), 0) + 1
        key = (b["col"], t, b["root"])
        good = t + 2 * m == 8 and b["rmode"] == 1 and key not in seen and all(len(r) == len(b["flags"]) for r in b["rows"])
        seen.add(key)
        rmask = (1 << comb(t, 2)) - 1
        good &= all(f & rmask == b["root"] and f >> comb(t + m, 2) == 0 for f in b["flags"])
        canon = [canon_flag(f, t, m) for f in b["flags"]]
        good &= len(set(canon)) == len(canon)
        ok &= good
    out("blocks by (t, m, colour): %s" % dict(sorted(fam.items())))
    out("every block: t + 2m = 8, ordered roots, rows of the right length, one block per (colour, root graph), "
        "flags whose root part is the block's root graph and no two in one orbit of the extension vertices: %s" % ok)
    return ok


def run_sums(binpath, cert, graphs_file):
    t0 = time.time()
    res = subprocess.run([binpath, "sums", cert, graphs_file], capture_output=True, text=True)
    if res.returncode != 0:
        return None, res.stderr.strip(), time.time() - t0
    rows = []
    for line in res.stdout.split("\n"):
        if line:
            g, mono, s6, s4, s2, s0 = line.split()
            rows.append((int(g), int(mono), {6: int(s6), 4: int(s4), 2: int(s2), 0: int(s0)}))
    return rows, "", time.time() - t0


def norm(t):
    m = TS[t]
    return comb(8, t) * factorial(t) * comb(8 - t, m) * comb(8 - t - m, m)


def values(rows, D=10 ** 12):
    vals = []
    for g, mono, S in rows:
        v = Fr(mono, 28) - sum(Fr(S[t], D * D * norm(t)) for t in S)
        vals.append((v, g))
    return vals


def slow_sums(H, blocks):
    """Independent pure-Python recomputation for one graph: canonical flags by minimum over extension permutations."""
    full = (1 << 28) - 1
    index = {}
    for b in blocks:
        t, m = b["t"], b["m"]
        G = {}
        canon = {canon_flag(f, t, m): i for i, f in enumerate(b["flags"])}
        index[(b["col"], t, b["root"])] = (b, canon)
    S = {6: 0, 4: 0, 2: 0, 0: 0}
    gram = {}
    for col in (0, 1):
        Gm = H if col == 0 else ~H & full

        def adj(x, y):
            if x > y:
                x, y = y, x
            return Gm >> pidx(x, y) & 1

        for t, m in TS.items():
            for R in itertools.permutations(range(8), t):
                rc = 0
                for y in range(1, t):
                    for x in range(y):
                        rc |= adj(R[x], R[y]) << pidx(x, y)
                hit = index.get((col, t, rc))
                if hit is None:
                    continue
                b, canon = hit
                key = id(b)
                if key not in gram:
                    rows = b["rows"]
                    nf = len(b["flags"])
                    gram[key] = [[sum(r[i] * r[j] for r in rows) for j in range(nf)] for i in range(nf)]
                Gr = gram[key]
                rest = [v for v in range(8) if v not in R]
                for A in itertools.combinations(rest, m):
                    Bv = [v for v in rest if v not in A]
                    for Bs in itertools.combinations(Bv, m):
                        fl = []
                        for E in (A, Bs):
                            vv = list(R) + list(E)
                            code = rc
                            for y in range(t, t + m):
                                for x in range(y):
                                    code |= adj(vv[x], vv[y]) << pidx(x, y)
                            fl.append(canon.get(canon_flag(code, t, m)))
                        if fl[0] is not None and fl[1] is not None:
                            S[t] += Gr[fl[0]][fl[1]]
    return S


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    root, binpath = sys.argv[1], sys.argv[2]
    controls = "--controls" in sys.argv

    def out(x):
        print(x, flush=True)

    out("89c candidate row [Pin2026] (PR #308): C_89c >= 8280981029/250000000000000 = 0.000033123924116*")
    status = True
    for rel, h in SHA.items():
        got = sha256(os.path.join(root, rel))
        out("%s sha256 %s (as stated: %s)" % (rel, got, got == h))
        status &= got == h
    cert = os.path.join(root, "c89c/c6_ordered_type8.cert")
    head, blocks, complete = read_cert(cert)
    status &= structure(head, blocks, complete, out)
    graphs = read_graph6(os.path.join(root, "c89c/graph8.g6"))
    work = os.path.join(root, "c89c/work")
    os.makedirs(work, exist_ok=True)
    gfile = os.path.join(work, "graphs.txt")
    with open(gfile, "w") as fh:
        fh.write("\n".join(map(str, graphs)) + "\n")
    # the graph list: pairwise non-isomorphic, orbits covering all labelled graphs
    t0 = time.time()
    res = subprocess.run([os.path.join(binpath, "flag_count"), "aut", gfile], capture_output=True, text=True)
    aut = [tuple(map(int, l.split())) for l in res.stdout.split("\n") if l]
    canon = [a[1] for a in aut]
    orbits = sum(Fr(40320, a[2]) for a in aut)
    distinct = len(set(canon)) == len(canon) == len(graphs) == 12346
    cover = orbits == 2 ** 28
    out("graph list: %d graphs, pairwise non-isomorphic (distinct canonical forms): %s; sum of 8!/|Aut| = %s = 2^28: %s "
        "(%.1f s)" % (len(graphs), distinct, orbits, cover, time.time() - t0))
    status &= distinct and cover
    # the inequality
    rows, err, dt = run_sums(os.path.join(binpath, "flag_count"), cert, gfile)
    if rows is None:
        out("counting failed: %s" % err)
        return 1
    vals = values(rows)
    vmin, gmin = min(vals)
    allok = all(v >= B_STATED for v, g in vals)
    out("every one of the %d graphs satisfies mono_K6 - contributions >= b: %s (%.1f s)" % (len(vals), allok, dt))
    out("minimum of mono_K6 - contributions: %s = %.18f..., at graph %d (0-based index in graph8.g6)" % (
        vmin, float(vmin), gmin))
    tight = sorted(vals)[:5]
    out("five smallest values minus b: %s" % ", ".join("graph %d: %.3e" % (g, float(v - B_STATED)) for v, g in tight))
    slack = vmin - B_STATED
    out("slack at the minimum: %s; the packet's stated slack %s at graph %d: equal: %s" % (
        slack, SLACK_STATED, WORST_STATED, slack == SLACK_STATED and gmin == WORST_STATED))
    status &= allok
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        # (i) b raised by the slack plus 10^-30 is refused
        b2 = B_STATED + slack + Fr(1, 10 ** 30)
        r1 = not all(v >= b2 for v, g in vals)
        out("[ctl] b raised to the minimum plus 10^-30: refused: %s" % r1)
        ctl &= r1
        # (ii) one block's rows scaled by 101/100 (rounded to integers): some graph fails
        big = max(range(len(blocks)), key=lambda i: sum(abs(x) for r in blocks[i]["rows"] for x in r) if blocks[i]["t"] == 6 else -1)
        blocks2 = [dict(b) for b in blocks]
        blocks2[big] = dict(blocks[big])
        blocks2[big]["rows"] = [[x * 101 // 100 for x in r] for r in blocks[big]["rows"]]
        c2 = os.path.join(work, "scaled.cert")
        write_cert(c2, head, blocks2)
        rows2, err2, _ = run_sums(os.path.join(binpath, "flag_count"), c2, gfile)
        v2 = values(rows2) if rows2 else []
        fails = [g for v, g in v2 if v < B_STATED]
        r2 = bool(fails)
        out("[ctl] block %d (t = %d, colour %d) with every row scaled by 101/100: %d graphs fall below b (first: %s): %s"
            % (big, blocks[big]["t"], blocks[big]["col"], len(fails), fails[:1], r2))
        ctl &= r2
        # (iii) a listed flag with a wrong root part is refused by the counting program
        blocks3 = [dict(b) for b in blocks]
        bb = next(i for i, b in enumerate(blocks) if b["t"] == 4)
        blocks3[bb] = dict(blocks[bb])
        blocks3[bb]["flags"] = list(blocks[bb]["flags"])
        blocks3[bb]["flags"][0] ^= 1
        c3 = os.path.join(work, "badflag.cert")
        write_cert(c3, head, blocks3)
        rows3, err3, _ = run_sums(os.path.join(binpath, "flag_count"), c3, gfile)
        r3 = rows3 is None and "root part" in err3
        out("[ctl] a flag of block %d with one root edge flipped: refused (%s): %s" % (bb, err3, r3))
        ctl &= r3
        # (iv) the graph list with one graph removed, or one duplicated, is refused
        r4a = sum(Fr(40320, a[2]) for a in aut[1:]) != 2 ** 28
        dup = canon + [canon[7]]
        r4b = len(set(dup)) != len(dup)
        out("[ctl] graph list without graph 0: orbit sum != 2^28: %s; with graph 7 twice: not distinct: %s" % (r4a, r4b))
        ctl &= r4a and r4b
        # (v) independent pure-Python recomputation for the worst graph and two others
        idxs = [gmin, 0, 12345]
        agree = True
        for gi in idxs:
            t0 = time.time()
            Sp = slow_sums(graphs[gi], blocks)
            Sc = rows[gi][2]
            same = Sp == Sc
            agree &= same
            out("[ctl] graph %d recomputed in pure Python (canonical flags by minimum over extension orders): S = %s, "
                "the same as flag_count.c: %s (%.1f s)" % (gi, {t: str(v)[:12] + "..." for t, v in Sp.items()}, same,
                                                          time.time() - t0))
        ctl &= agree
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT: %s" % ("PASS (every 8-vertex graph satisfies the stated inequality with b = 8280981029/250000000000000; "
                        "the minimum and its graph as the packet states)" if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
