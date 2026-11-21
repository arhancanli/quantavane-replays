#!/usr/bin/env python3
"""Program A (Python, exact integers): upper-bound side of Rohan Arun's Erdos problem 838 minima
(github.com/rohanarun/erdos-838-convex-polygon-minima @ c301008).

For each published integer-coordinate witness: the points are distinct and no three are collinear
(integer orientation determinant != 0), and the number of subsets of size >= 3 in convex position is
counted.  Criterion (method A): a finite point set in general position is in convex position iff no
point lies inside a triangle formed by three others, i.e. iff it contains no 'bad' 4-subset (one point
strictly inside the triangle of the other three).  Subsets are processed in increasing bitmask order:
S is convex iff S minus its highest point is convex and no bad 4-subset with that highest point lies in S.

Witness specs (any number): PATH[:FORMAT[:INDEX]]=TOTAL[,b3,b4,...][@d1/d2/...]
  FORMAT: xy (two integers per line, default), jsonl (record INDEX of a JSON Lines file
  [sign_string, [[x, y], ...]]; the sign string is also checked: bit 1 = positive determinant, triples
  in ascending lexicographic order), lean (lines '| i => (x, y)').
  TOTAL: claimed number of convex subsets of size >= 3; b3,b4,...: claimed counts by size from 3;
  @d1/d2/...: claimed sorted totals after deleting each single point.
Options: --control collinear|move (applied to the first witness: point 2 replaced by 2*P1 - P0, which is
collinear with P0 and P1; or point 0 moved by (-1000, 37)).
Prints 'A: VERDICT PASS' (exit 0) or 'A: VERDICT REFUSED' (exit 1).
"""
import sys
import json
import re
from itertools import combinations


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def read_points(path, fmt, index):
    sign = None
    if fmt == "xy":
        pts = []
        for line in open(path):
            t = line.split()
            if t:
                pts.append((int(t[0]), int(t[1])))
    elif fmt == "jsonl":
        recs = [json.loads(l) for l in open(path) if l.strip()]
        sign, raw = recs[index][0], recs[index][1]
        pts = [(int(x), int(y)) for x, y in raw]
    elif fmt == "lean":
        pts = []
        for line in open(path):
            m = re.match(r"\s*\|\s*(\d+)\s*=>\s*\((-?\d+),\s*(-?\d+)\)", line)
            if m:
                assert int(m.group(1)) == len(pts)
                pts.append((int(m.group(2)), int(m.group(3))))
    else:
        raise ValueError(fmt)
    return pts, sign


def analyse(pts):
    n = len(pts)
    collinear = [t for t in combinations(range(n), 3) if orient(*(pts[i] for i in t)) == 0]
    bad_by_top = [[] for _ in range(n)]
    for quad in combinations(range(n), 4):
        for k in range(4):
            p = quad[k]
            a, b, c = [quad[j] for j in range(4) if j != k]
            o = orient(pts[a], pts[b], pts[c])
            if o == 0:
                continue
            s1 = orient(pts[a], pts[b], pts[p]); s2 = orient(pts[b], pts[c], pts[p]); s3 = orient(pts[c], pts[a], pts[p])
            if (o > 0 and s1 > 0 and s2 > 0 and s3 > 0) or (o < 0 and s1 < 0 and s2 < 0 and s3 < 0):
                m = (1 << quad[0]) | (1 << quad[1]) | (1 << quad[2]) | (1 << quad[3])
                bad_by_top[quad[3]].append(m)
                break
    N = 1 << n
    convex = bytearray(N)
    convex[0] = 1
    for S in range(1, N):
        t = S.bit_length() - 1
        if convex[S ^ (1 << t)]:
            ok = 1
            for m in bad_by_top[t]:
                if S & m == m:
                    ok = 0
                    break
            convex[S] = ok
    bysize = [0] * (n + 1)
    dele = [0] * n
    for S in range(N):
        if convex[S]:
            k = bin(S).count("1")
            if k >= 3:
                bysize[k] += 1
                for i in range(n):
                    if not S >> i & 1:
                        dele[i] += 1
    return collinear, bysize, dele


def main():
    args = sys.argv[1:]
    control = None
    if "--control" in args:
        i = args.index("--control"); control = args[i + 1]; del args[i:i + 2]
    fails = []
    for wi, spec in enumerate(args):
        lhs, rhs = spec.rsplit("=", 1)
        parts = lhs.split(":")
        path = parts[0]; fmt = parts[1] if len(parts) > 1 else "xy"; index = int(parts[2]) if len(parts) > 2 else 0
        dels = None
        if "@" in rhs:
            rhs, d = rhs.split("@"); dels = [int(x) for x in d.split("/")]
        nums = [int(x) for x in rhs.split(",")]
        total, breakdown = nums[0], nums[1:]
        pts, sign = read_points(path, fmt, index)
        tag = "%s%s" % (path.split("/")[-1], ("[%d]" % index) if fmt == "jsonl" else "")
        if wi == 0 and control == "collinear":
            pts[2] = (2 * pts[1][0] - pts[0][0], 2 * pts[1][1] - pts[0][1])
            print("A: CONTROL collinear: %s point 2 replaced by 2*P1 - P0 = %s" % (tag, pts[2]))
        if wi == 0 and control == "move":
            pts[0] = (pts[0][0] - 1000, pts[0][1] + 37)
            print("A: CONTROL move: %s point 0 moved to %s" % (tag, pts[0]))
        n = len(pts)
        if len(set(pts)) != n:
            fails.append("%s: repeated point" % tag)
        collinear, bysize, dele = analyse(pts)
        if collinear:
            fails.append("%s: %d collinear triples, e.g. %s" % (tag, len(collinear), collinear[0]))
        tot = sum(bysize[3:])
        if sign is not None:
            mine = "".join("1" if orient(pts[a], pts[b], pts[c]) > 0 else "0" for a, b, c in combinations(range(n), 3))
            if mine != sign:
                fails.append("%s: determinant signs differ from the record's sign string" % tag)
        prof = ",".join(str(x) for x in bysize[3:])
        print("A: %s: n=%d, collinear triples %d, convex subsets by size 3..%d = %s, total %d (claimed %d)%s%s"
              % (tag, n, len(collinear), n, prof, tot, total,
                 "" if sign is None else ", sign string %s" % ("matches" if sign == mine else "DIFFERS"),
                 "" if dels is None else ", deletion totals %s" % "/".join(map(str, sorted(dele)))))
        if tot != total:
            fails.append("%s: total %d != claimed %d" % (tag, tot, total))
        if breakdown and bysize[3:3 + len(breakdown)] != breakdown or any(bysize[3 + len(breakdown):]) and breakdown:
            fails.append("%s: breakdown %s != claimed %s" % (tag, bysize[3:], breakdown))
        if dels is not None and sorted(dele) != dels:
            fails.append("%s: deletion totals %s != claimed %s" % (tag, sorted(dele), dels))
    for x in fails[:20]:
        print("A: FAIL " + x)
    if fails:
        print("A: VERDICT REFUSED (%d failures)" % len(fails))
        return 1
    print("A: VERDICT PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
