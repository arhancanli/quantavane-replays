# 13a: the covering certificate of the upper bound [K2026], $C_{13a} \le 0.2463322372\ldots$\* (14 November 2026)

The upper-bound row that [#314](https://github.com/teorth/optimizationproblems/pull/314) adds to the 13a page:
Keller, arXiv 2609.21968v1, Theorem 1.1, $C_{13a} \le \mathrm{area}(K) = 7170601298323360123534337/29109471743520000000000000
= 0.2463322372\ldots$, where $K$ is the quadrilateral of (4.1) and Theorem 4.1 says that $K$ covers every unit arc.
On 8 November we replayed $K$ and its area (`../c13a`); the covering proof is a finite certificate whose encoding was
not stated anywhere we could read. On 14 November the author stated on #314 the three missing definitions (the 60
contact directions, the matrix of each orientation spec, and the rows with the expansion rule of $M n_i$). This
folder replays the whole certificate from those definitions and the paper (Sections 4.1 to 4.4, read as text).

Certificate: ethan-keller/moser-worm-improved-bounds @ `9a33941` (the commit the paper cites),
`certificates/upper/`: `header.json`, `tree.json.gz` and 50 multiplier files (23 MB), read as data, each file's
SHA-256 equal to the certificate's own `manifest.json`. No program of that repository (Python or Lean) was read or
run.

Two separately written exact programs:

- `c13acov_a.py` (A): Python fractions; directions ordered by half-plane and $x$-coordinate; tree splits checked by
  areas and pairwise disjointness; per leaf the vectors $f_\lambda$ brought to one denominator and every family sum
  accumulated as a running sum.
- `c13acov_b.py` (B): each row stored once as integer vectors over its own denominator, each leaf summed in integers
  over the least common multiple of its rows' denominators; directions ordered by the half-angle tangent
  $y/(1+x)$; splits checked by marking a grid of $(p, q)$ points; family sums as differences of cumulative sums;
  the area of $K$ as $\frac12\sum_i h_i \ell_i$ (A uses the shoelace formula).

| check | result (both programs) |
|---|---|
| side data and $K$ | $\sigma$, $h$ as in (4.1); the area is the stated fraction ($0.24633223720109565\ldots$; 0.24634 and 0.2463322373 are upper bounds) |
| normal relations | $a_{20} = 5541045/8027539$, $a_{21} = 24427652/40137695$, $a_{30} = 47299645/118570699$, $a_{31} = 502233812/592853495$, positive, equal to the values the author states |
| specs | 1,463 (16 `flush`, 9 `rot`, 1,438 `vrot`), every matrix orthogonal |
| contact directions | the 25 `flush`/`rot` specs give $(0, \pm1)$ and exactly 30 + 30 rational unit directions $u^R$, $u^L$, strictly ordered in their open half-planes |
| expansions | every $M_s n_i$ of a `flush`/`rot` spec is a label direction; every other one of a `vrot` spec is $\kappa_1 w_1 + \kappa_2 w_2$ with adjacent label directions and $\kappa \ge 0$, exactly (5,744 two-term expansions) |
| tree | 24,898 nodes, 12,663 leaves, depth 26 (the paper's table); every split partitions its rectangle; no spec branched twice on a path; root $[0,30]^2$ |
| leaves | all 12,663 pass (4.7): every multiplier $> 0$, escape rows only from the leaf's path, $B \ge 1$, $\sum_\lambda f_\lambda = 0$, and $\lvert\sum_{\lambda \in I} f_\lambda\rvert^2 \le 1$ for every $I$ in the four families $\mathcal F(\mathcal R)$; 6,778,020 multipliers (the manifest's count) |
| margins | the largest squared norm over all leaves and families is the same fraction in both programs, $0.99999985031\ldots$; the least $B$ is exactly 1 |

So the covering certificate is accepted exactly as the paper and the author's definitions describe it, and together
with the area this gives the row's bound.

**Red controls**, each refused (a control stops at its first refusal, as the paper's checker rejects). A: $h_2$ in
the escape rows one unit lower in its last place (leaf 0, where $B = 1$, fails); the first two right directions swapped;
the root split's first cell widened so that it overlaps the second; leaf 0 (where $B = 1$ exactly) with its escape
multiplier one unit lower; leaf 0 given an escape row of a spec not on its path. B: $h_3$ in the escape rows one
unit lower (leaf 5 fails); the direction $u^R_6$ reflected in the $x$-axis (a `flush`/`rot` direction is then no label
direction); the root split's first cell narrowed (a gap); every multiplier of leaf 0
doubled ($B = 2$ and the sum still zero, but a family sum's squared norm exceeds 1); a direction $(0, -1)$ from a
spec given to $P_-$ instead of $P_+$, against the stated rule (leaf 9 is the first to fail).

**Not checked (mathematics, read):** Lemmas 4.2 to 4.4 and the reduction of Section 4.3 (finite witness, normal
form, translation criterion, prefix families), which turn an accepted certificate into Theorem 4.1; the author's
Lean formalization, which checks them, was not built. Two readings of the definitions that the data decide: a
direction $(0, -1)$ inside a two-term expansion is given to $P_+$, as for one-term expansions; and the four families
are taken with $k_R = k_L = 30$ exactly as in the paper's table. Under these readings every leaf passes.

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+, standard library only; about three minutes, at most two
processes, under 1 GB). Output: `results/k2026cov.txt`.
