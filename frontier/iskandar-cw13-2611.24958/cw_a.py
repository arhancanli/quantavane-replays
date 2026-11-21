#!/usr/bin/env python3
"""Program A (Python, exact): the thirteen binary constant-weight codes of P. R. A. Iskandar,
arXiv 2611.24958v1 (ancillary files read as data only; the paper's verifiers were not read or run).

For each code A(n,d,w) >= M listed in CLAIMS.csv:
  * every line of codes/A{n}_{d}_{w}_{M}.txt is a 0/1 string of length n and weight w;
  * no word repeats and there are exactly M words;
  * minimum distance >= d, by the t-subset method: with t = w - d/2 + 1, every t-subset of the support
    of every word is stored in one hash set; a t-subset met twice means two words share t points, i.e.
    distance <= 2(w - t) = d - 2.  (Two distinct words at distance < d share >= t points.)
  * orbit rebuild: the coordinates are split into consecutive blocks of the manifest's lengths; g shifts
    every block of length > 1 by one place; the union of the <g>-orbits of the representatives in
    reps/... must equal the full word list as a set, the orbits must be disjoint (sizes sum to M) and the
    orbit-size distribution must equal the manifest's.
  * SHA-256 of each code file equals CLAIMS.csv.
Usage: cw_a.py ANC_DIR [--control distance|repeat|weight]   (the control edits the first code in memory)
Prints 'A: VERDICT PASS' (exit 0) or 'A: VERDICT REFUSED' (exit 1).
"""
import sys
import os
import csv
import hashlib
from itertools import combinations


def main():
    args = sys.argv[1:]
    control = None
    if "--control" in args:
        i = args.index("--control"); control = args[i + 1]; del args[i:i + 2]
    anc = args[0]
    claims = list(csv.DictReader(open(os.path.join(anc, "CLAIMS.csv"))))
    manifest = {(int(r["n"]), int(r["d"]), int(r["w"]), int(r["M"])): r
                for r in csv.DictReader(open(os.path.join(anc, "ORBIT_MANIFEST.csv")))}
    fails = []
    totw = totr = 0
    for ci, c in enumerate(claims):
        n, d, w = [int(x) for x in c["parameter"].strip('"A()').split(",")]
        M = int(c["new_lower_bound"])
        path = os.path.join(anc, c["code_file"])
        raw = open(path, "rb").read()
        if hashlib.sha256(raw).hexdigest() != c["code_sha256"]:
            fails.append("A(%d,%d,%d): SHA-256 differs from CLAIMS.csv" % (n, d, w))
        lines = [l.strip() for l in raw.decode("ascii").splitlines() if l.strip()]
        if ci == 0 and control:
            w0 = list(lines[0])
            if control == "distance":   # move one 1 of word 0 to a 0 position: distance 2 from word 0
                i1 = w0.index("1"); i0 = w0.index("0"); w0[i1], w0[i0] = "0", "1"
                lines[1] = "".join(w0)
            elif control == "repeat":
                lines[1] = lines[0]
            elif control == "weight":
                i0 = w0.index("0"); w0[i0] = "1"; lines[1] = "".join(w0)
            print("A: CONTROL %s: A(%d,%d,%d) word 2 replaced by %s" % (control, n, d, w, lines[1]))
        supports = []
        bad_len = bad_w = 0
        for s in lines:
            if len(s) != n or set(s) - {"0", "1"}:
                bad_len += 1; continue
            sup = tuple(i for i, ch in enumerate(s) if ch == "1")
            if len(sup) != w:
                bad_w += 1
            supports.append(sup)
        nd = len(set(supports))
        if bad_len: fails.append("A(%d,%d,%d): %d lines not 0/1 strings of length %d" % (n, d, w, bad_len, n))
        if bad_w: fails.append("A(%d,%d,%d): %d words of weight != %d" % (n, d, w, bad_w, w))
        if nd != len(supports): fails.append("A(%d,%d,%d): %d repeated words" % (n, d, w, len(supports) - nd))
        if len(lines) != M: fails.append("A(%d,%d,%d): %d words, claimed %d" % (n, d, w, len(lines), M))
        t = w - d // 2 + 1
        seen = {}
        clashes = 0
        for idx, sup in enumerate(set(supports)):
            for sub in combinations(sup, t):
                if sub in seen:
                    clashes += 1
                else:
                    seen[sub] = idx
        if clashes:
            fails.append("A(%d,%d,%d): %d t-subsets (t=%d) lie in two words: distance < %d" % (n, d, w, clashes, t, d))
        # orbit rebuild
        m = manifest[(n, d, w, M)]
        blocks = [int(x) for x in m["blocks"].split()]
        if sum(blocks) != n:
            fails.append("A(%d,%d,%d): manifest blocks do not sum to n" % (n, d, w))
        perm = list(range(n))
        start = 0
        for b in blocks:
            for j in range(b):
                perm[start + j] = start + (j + 1) % b
            start += b
        reps = [l.strip() for l in open(os.path.join(anc, m["rep_file"])) if l.strip()]
        union = set()
        sizes = {}
        total = 0
        for r in reps:
            sup = frozenset(i for i, ch in enumerate(r) if ch == "1")
            orb = set()
            cur = sup
            while cur not in orb:
                orb.add(cur)
                cur = frozenset(perm[i] for i in cur)
            sizes[len(orb)] = sizes.get(len(orb), 0) + 1
            total += len(orb)
            union |= orb
        want_sizes = dict(tuple(int(y) for y in x.split(":")) for x in m["orbit_sizes"].split(","))
        full = set(frozenset(s) for s in supports)
        same = union == full
        if not same: fails.append("A(%d,%d,%d): orbit union differs from the word list" % (n, d, w))
        if total != len(union) or total != M: fails.append("A(%d,%d,%d): orbits not disjoint or wrong total" % (n, d, w))
        if sizes != want_sizes: fails.append("A(%d,%d,%d): orbit sizes %s != manifest %s" % (n, d, w, sizes, want_sizes))
        if len(reps) != int(m["rep_count"]): fails.append("A(%d,%d,%d): rep count" % (n, d, w))
        totw += len(lines); totr += len(reps)
        print("A: A(%d,%d,%d) >= %d: words %d, distinct %d, length/weight ok %s, t=%d subsets clashing %d, "
              "orbit rebuild (%d reps, |<g>| from blocks %s, sizes %s) equals list: %s"
              % (n, d, w, M, len(lines), nd, "yes" if not (bad_len or bad_w) else "NO", t, clashes,
                 len(reps), m["blocks"], ",".join("%d:%d" % kv for kv in sorted(sizes.items())), "yes" if same else "NO"))
    print("A: total words %d, total representatives %d" % (totw, totr))
    for x in fails[:20]:
        print("A: FAIL " + x)
    if fails:
        print("A: VERDICT REFUSED (%d failures)" % len(fails))
        return 1
    print("A: VERDICT PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
