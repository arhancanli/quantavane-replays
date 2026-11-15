#!/usr/bin/env python3
"""compare.py OUT_A OUT_B: every line 'NAME in [lo, hi]' printed by both programs of a pair must give overlapping
intervals (two rigorous enclosures of the same number). Prints one line per quantity; exit 1 on any disjoint pair."""
import re
import sys
from fractions import Fraction as Fr

PAT = re.compile(r"^(.*?) in \[([-+0-9.Ee]+), ([-+0-9.Ee]+)\]")


def grab(path):
    out = {}
    for line in open(path):
        m = PAT.match(line.strip())
        if m:
            out[m.group(1)] = (Fr(m.group(2)), Fr(m.group(3)))
    return out


a, b = grab(sys.argv[1]), grab(sys.argv[2])
bad = 0
common = [k for k in a if k in b]
for k in common:
    (l1, h1), (l2, h2) = a[k], b[k]
    ok = max(l1, l2) <= min(h1, h2)
    bad += not ok
    print("%s: %s" % (k, "the two enclosures overlap" if ok else "DISJOINT"))
print("compared %d quantities, %d disjoint" % (len(common), bad))
sys.exit(1 if bad or not common else 0)
