# Erdős problem 838: minimum numbers of convex polygons on 12, 13 and 14 points, the upper-bound side (Rohan Arun)

Checked by two independent programs; not peer reviewed. Partial check: the witnesses (upper bounds) only.

## The claim (credit first)

Rohan Arun (Super Powers AI), "Minimum numbers of convex polygons on twelve through fourteen points",
[rohanarun/erdos-838-convex-polygon-minima](https://github.com/rohanarun/erdos-838-convex-polygon-minima). The first
commit is from 20 November 2026 05:52 UTC (12 and 13 points). The 14-point result followed at 06:58 UTC, and the check
uses commit `c301008` (13:52 UTC). Among n points in the plane with no three collinear, the minimum number of subsets of
size at least 3 in convex position is 408 for n = 12, 570 for n = 13 (five order types) and 780 for n = 14 (one order
type, with 364 triangles, 350 quadrilaterals, 64 pentagons and 2 hexagons; removing any point leaves at least 574). The
repository also re-derives Nikol Savova's values 189 (n = 10) and 283 (n = 11). The author writes: "OpenAI Codex
assisted with mathematical derivations, code, certificate generation, verification, and manuscript preparation."

## One command

```sh
./fetch.sh && ./run_all.sh
```

or, from the repository root, `./check.sh frontier/erdos838-convex-polygons` (from this folder:
`../../check.sh frontier/erdos838-convex-polygons`). `fetch.sh` downloads seven files from the pinned commit into
`certs/` and checks their SHA-256 against `SHA256SUMS`: the integer-coordinate witnesses for 10, 11, 12 and 14 points,
the five 13-point realizations in `real13.jsonl`, and `Witness13.lean`, which is read only for its list of coordinates.
Python 3.9+ and a C compiler; a few seconds. Output: `results/erdos838.txt`. No script, C++ or Lean file of the
repository was run.

## Verdict: partial. The upper-bound side is confirmed

Every published witness has no three collinear points and exactly the claimed number of convex polygons. Both programs
give the same counts by size:

| witness | points | triangles, quadrilaterals, pentagons, hexagons | total | claimed |
|---|---|---|---|---|
| `witness10.txt` | 10 | 120, 66, 3, 0 | 189 | 189 (Savova) |
| `witness11.txt` | 11 | 165, 108, 10, 0 | 283 | 283 (Savova) |
| `witness12.txt` | 12 | 220, 165, 22, 1 | 408 | 408 |
| `real13.jsonl`, record 0 | 13 | 286, 241, 41, 2 | 570 | 570 |
| `real13.jsonl`, records 1 to 4 | 13 | 286, 247, 35, 2 (each) | 570 (each) | 570 |
| `Witness13.lean` | 13 | 286, 241, 41, 2 | 570 | 570 |
| `witness780-small.txt` | 14 | 364, 350, 64, 2 | 780 | 780 = 364 + 350 + 64 + 2 |
| `witness780.txt` (coordinates up to 4 * 10^12) | 14 | 364, 350, 64, 2 | 780 | 780 |

No set has a convex heptagon or larger. Deleting each single point of the 14-point witness (both files) leaves
574, 574, 574, 575, 575, 575, 576, 578, 578, 580, 581, 582, 586 and 588 convex polygons. This is the repository's list,
so the smallest is 574. Program A also confirms that the sign strings in `real13.jsonl` are the determinant signs of
their coordinates. So the minima are at most 189, 283, 408, 570 and 780.

## The two programs (written separately, no shared code)

- **A** = `erdos838_a.py` (Python, exact integers). It tests collinearity with the integer orientation determinant. A set
  with no three collinear is in convex position exactly when no point lies inside a triangle of three others, and A
  evaluates this for every subset by a bitmask recursion.
- **B** = `erdos838_b.c` (C, 128-bit orientation). It builds the convex hull of every subset by Andrew's monotone chain.
  A subset is in convex position exactly when its hull has as many vertices as the subset has points.

## Red controls (both programs refuse each)

1. and 2. Point 2 replaced by 2 P1 - P0, which is collinear with P0 and P1, in the 12-point and in the 14-point witness.
3. and 4. Point 0 moved by (-1000, +37). The 12-point total becomes 477 and the 14-point total 827.
5. The 12-point witness claimed to have 407.

## Not checked

The lower bounds, that no configuration has fewer convex polygons. These rest on the enumeration of candidate order
types (deletion averaging and radial-sign extension, the C++ runs) and on the Grassmann-Plücker product certificates and
Lean files, and we did not check them. Also not checked: the classification (five 13-point order types, one 14-point
order type), whether the five 13-point realizations are distinct order types, and the exterior-twin theorem.

## Authorship and the role of AI

This check is the work of Arhan Canli (Quantavane). The programs, runs and this README were written by AI agents
(Claude, by Anthropic) working for him. The two programs were written separately. Nothing is claimed beyond what they
check, and the mathematics is the claimant's. Questions and corrections: arhan@quantavane.com.
