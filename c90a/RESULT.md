# 90a: the page's elementary covering system, $C_{90} \ge 2$ (14 November 2026)

Page: [90a, the maximum minimum modulus of a distinct covering system](https://teorth.github.io/optimizationproblems/constants/90a.html)
(main at `87b4ebb`). The lower-bound row "2, elementary example above" refers to the system printed in the
description:

$$0 \pmod 2,\quad 0 \pmod 3,\quad 1 \pmod 4,\quad 1 \pmod 6,\quad 11 \pmod{12}.$$

## What was checked

With two separately written programs, both exact:

- `c90a_cover.py` (A): the moduli 2, 3, 4, 6, 12 are greater than 1 and pairwise distinct; the period is 12; every
  residue 0, ..., 11 lies in at least one class. So the system is a distinct covering system with least modulus 2.
- `c90a_density.py` (B): no residue is enumerated. The density of the union is computed by inclusion-exclusion over
  all 31 non-empty sub-families (an intersection has density 0 if two classes are incompatible modulo the gcd of
  their moduli, and 1/lcm otherwise), in exact rationals: it is exactly 1. A periodic set of density 1 is all of
  $\mathbb{Z}$.

Result: the system covers $\mathbb{Z}$, so $C_{90} \ge 2$, as the row says. There is no external input (the object
is five residue classes printed on the page), hence no download and no hash.

## Red controls (all refused by both programs)

- 11 mod 12 replaced by 10 mod 12: the residue 11 is uncovered (density 11/12).
- 1 mod 6 replaced by 5 mod 6: the residue 7 is uncovered (density 11/12).
- The class 1 mod 6 removed: the residue 7 is uncovered.
- 0 mod 2 with 1 mod 2: it covers, but the moduli are not distinct.

Full output: `results/elementary.txt`.

## What was NOT checked

Only this one row. The other lower-bound rows of the page (9, 18, 20, 24, 25, 40, 42 and the starred 44) are
constructions with very many congruences, described in theses and papers or in a Lean development; none of them is
checked here. No upper bound is checked.

## Reproduce

`./run_all.sh` (Python 3.9+, standard library only; a second).
