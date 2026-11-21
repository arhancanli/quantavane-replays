"""Red-control inputs for program A: copies of the colour strings with one forbidden set planted (data edits only)."""
import itertools, sys, os
src, dst = sys.argv[1], sys.argv[2]
os.makedirs(dst, exist_ok=True)
def plant(fname, n, t, S, colour, out):
    s = list(open(os.path.join(src, fname)).read().strip())
    idx = {c: i for i, c in enumerate(itertools.combinations(range(n), t))} if t == 3 else None
    if t == 3:
        for c in itertools.combinations(sorted(S), 3):
            s[idx[c]] = str(colour)
    open(os.path.join(dst, out), "w").write("".join(s))
plant("c37.txt", 37, 3, [0, 5, 9, 13], 0, "c37_plant.txt")          # a colour-0 K4
plant("c70.txt", 70, 3, [1, 7, 20, 33, 50, 64], 1, "c70_plant.txt")  # a colour-1 K6
plant("c98.txt", 98, 3, [2, 11, 40, 63, 97], 0, "c98_plant.txt")    # a colour-0 K5 through infinity
# 5-subsets: in the short form, every 5-set through 0,1 of the 6-set {0,1,10,20,30,40} set to colour 1; the
# colouring then fails the 20-normalisation consistency (program A must report inconsistencies)
s = list(open(os.path.join(src, "s83.txt")).read().strip())
rank = {c: i for i, c in enumerate(itertools.combinations(range(2, 83), 3))}
for c in itertools.combinations([10, 20, 30, 40], 3):
    s[rank[c]] = "1"
open(os.path.join(dst, "s83_plant.txt"), "w").write("".join(s))
t = list(open(os.path.join(src, "t599.txt")).read().strip()); t[2] = "1" if t[2] == "0" else "0"
open(os.path.join(dst, "t599_flip.txt"), "w").write("".join(t))
