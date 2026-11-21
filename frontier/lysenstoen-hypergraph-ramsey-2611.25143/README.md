# Six hypergraph Ramsey lower bounds of Christian Lysenstøen (arXiv 2611.25143), checked

**The result is Christian Lysenstøen's.** In *New lower bounds for six small hypergraph Ramsey numbers*,
[arXiv 2611.25143v1](https://arxiv.org/abs/2611.25143) (math.CO, submitted 18 November 2026, 20:23 UTC), he gives
explicit 2-colourings proving

| Bound | Previous bound (as the paper states) | Colouring |
| --- | --- | --- |
| $R(4,5;3) \ge 38$ | 35 (Radziszowski's dynamic survey) | triples of 37 points, invariant under $x \mapsto 10^k x + b$ on $\mathbb{Z}_{37}$ |
| $R(4,6;3) \ge 71$ | 63 | triples of 70 points ($\mathbb{Z}_{69}$ and one fixed point) |
| $R(5,5;3) \ge 99$ | 88 | triples of 98 points ($\mathbb{Z}_{97}$ and one fixed point) |
| $R(5,6;4) \ge 201$ | 67 (survey), 159 (a 2026 data deposit) | 4-subsets of the projective line over $\mathbb{F}_{199}$, $\mathrm{PSL}(2,199)$-invariant |
| $R(6,6;4) \ge 601$ | 133 (survey), 317 (a 2026 data deposit) | 4-subsets of the projective line over $\mathbb{F}_{599}$, given by a table of cross-ratios |
| $R(6,6;5) \ge 84$ | 72 | 5-subsets of $\mathbb{Z}_{83}$, invariant under every affine map |

The colourings are the paper's ancillary files, read here as data. The paper's own programs (`verify.py`,
`verify_crossratio.c`, `verify_projective.cpp`, `rebuild_*.py`, `check_projective.py`) were neither read nor run.

## Verdict: CONFIRMED, all six (checked by two independent programs; not peer reviewed)

Verdict time 21 November 2026, 15:28 UTC: about 38 hours after the paper was announced (arXiv listing of 20 November,
01:00 UTC) and 67 hours after submission.

| Colouring | Program A (C) | Program B (JavaScript) |
| --- | --- | --- |
| 37 points, triples | no colour-0 $K_4^{(3)}$, no colour-1 $K_5^{(3)}$ (exhaustive, no symmetry) | rebuilt from the 32 orbit representatives and 111 maps: equals the stored string; invariant; nothing through point 0 |
| 70 points, triples | no colour-0 $K_4^{(3)}$, no colour-1 $K_6^{(3)}$ (exhaustive) | rebuilt (27 representatives, 759 maps), invariant, nothing through 0 |
| 98 points, triples | no monochromatic $K_5^{(3)}$ (exhaustive) | rebuilt (36 representatives, 2,328 maps), invariant, nothing through 0 |
| 83 points, 5-subsets | full colouring built from the short form; all 20 affine normalisations of each of the 29,034,396 5-subsets agree; all 377,447,148 6-sets examined, none monochromatic | full colouring's SHA-256 equals the stated digest `b196450e...`; affine consistency on the 85,320 sets through 0, 1; all 1,663,740 6-sets through {0, 1}: none monochromatic |
| 200 points, 4-subsets | invariance of all 64,684,950 4-subsets under $x+1$ and $-1/x$; the orbit of $\{0,1,\infty\}$ is all 1,313,400 3-subsets; nothing of size 5 (colour 0) or 6 (colour 1) through $\{0,1,\infty\}$ | colouring from the rule has the stated SHA-256 `621a303b...`; invariance on all 4-subsets; union-find: all 3-subsets in one class; nothing through $\{0,1,2\}$ |
| 600 points, 4-subsets | the table is constant on the six-element cross-ratio classes (so the colouring is $\mathrm{PGL}(2,599)$-invariant); the orbit of $\{0,1,\infty\}$ is all 35,820,200 3-subsets; no monochromatic 6-set through $\{0,1,\infty\}$ | table constant on the orbits of $x \mapsto 1/x$, $x \mapsto 1-x$; no monochromatic 6-set through $\{0,1,2\}$ ($\mathrm{PGL}(2,q)$ is 3-transitive) |

The two programs use different reductions: A searches the three triple colourings and the 5-subset colouring with no
symmetry at all, and for the two projective colourings fixes $\{0,1,\infty\}$ after checking invariance and
transitivity itself; B rebuilds every colouring from its stated rule, checks the stated digests, and searches through
other fixed points ($0$; $\{0,1\}$; $\{0,1,2\}$).

**Red controls**, all refused: A on a colour-0 $K_4$ planted in the 37-point colouring, a colour-1 $K_6$ in the 70-point
one, a colour-0 $K_5$ through the fixed point in the 98-point one, an inconsistent short form for 83 points, one table
entry flipped for 600 points, and the 200-point colouring asked for size 4 instead of 5 (85 colour-0 4-sets found);
B on planted sets in the four symmetric colourings and one table entry flipped for 600 and for 200 points.

**Not checked:** that the previous bounds are as the paper states (we did not re-read the survey for each cell), and
the paper's further consequences ($R(5,5,5;3) \ge 9605$, $R_k(5;3) \ge 98^{2^{k-2}}+1$, the stepped-up 4-uniform bounds),
which follow from known constructions.

## Check it yourself

```sh
./fetch.sh && ./run_all.sh      # or, from the repository root: ./check.sh frontier/lysenstoen-hypergraph-ramsey-2611.25143
```

`fetch.sh` downloads the six data files from the arXiv ancillary files of version 1 and checks them against
`SHA256SUMS`. `run_all.sh` compiles `lys_a.c`, runs both programs on all six colourings and then the red controls, and
exits 0 only if everything agrees (about two minutes; C compiler, Python 3, node 18 or later). Outputs are in `results/`.

## Authorship and the role of AI

This check is the work of Arhan Canli (Quantavane). The programs, the runs and this README were written by AI agents
(Claude, by Anthropic) working for him. The two programs were written separately, in different languages, and share
no code. Nothing is claimed beyond what they check; the mathematics and the colourings are Christian Lysenstøen's.
Questions and corrections: arhan@quantavane.com.

Data and text: CC BY 4.0 for our files; the colourings are the author's.
