# Euclidean TSP (Beardwood-Halton-Hammersley) constant: beta >= 0.64254, a sharper first term in Dong and Cao's series

Constant 12 of the [optimization constants repository](https://teorth.github.io/optimizationproblems/constants/12a.html):
the length of the shortest tour through n uniform random points of the unit square is asymptotic to beta sqrt(n).

## The result

Dong and Cao ([arXiv:2610.02809](https://arxiv.org/abs/2610.02809), Theorem 3.11) prove, for every K >= 3,

$$\beta \geq \frac58 + \sum_{k=3}^{K} C_k,\qquad C_k = \frac{\Gamma(k-\frac12)}{2\,k!}\int_{(\mathbb R^2)^{k-1}} \frac{g(G_1(y))}{A(y)^{k-1/2}}\,dy_2\cdots dy_k ,$$

where $y_1 = 0$, $G_1(y)$ joins two points at distance at most 1, $g = q - D$ is their lower bound for the excess of
the path-cover number of a connected component over its degree count, and $A(y)$ is the area of the union of the unit
discs centred at the $y_i$. Each $C_k$ is bounded from below on its own; they certify $C_3 \geq 0.01113528859$,
$C_4 \geq 0.005040573276$ and $C_5 \geq 0.001015487669$, giving $\beta \geq 0.6421913$.

We certify a sharper value of the first term:

$$C_3 \geq 0.0114906 \quad\text{(program B: } 0.011490649265622\text{; program A: } 0.011490649266\text{)},$$

so that, with Dong and Cao's own $C_4$ and $C_5$ unchanged,

$$\beta \geq \frac58 + 0.0114906 + 0.005040573276 + 0.001015487669 \geq 0.64254.$$

The best lower bound listed on the repository's page today is 0.6277 (Gaudio and Jaillet 2020); Dong and Cao's
0.6421913 is not listed there. The improvement over Dong and Cao is about $3.5 \cdot 10^{-4}$. The lower-bound formula
(their Theorem 3.11) and the values of $C_4$ and $C_5$ (their Section 3.3.1, Table 1) are Dong and Cao's; only the
evaluation of $C_3$ is new. The combined bound relies on their $C_4$ and $C_5$ computations, whose code is not public
("available upon request"); our $C_3$ alone gives $5/8 + C_3 \geq 0.6364906$.

## Why C_3 is a three-dimensional integral

For three points, a connected $G_1$ is a path ($q = D = 1$, so $g = 0$) or a triangle ($q = 1$, $D = 0$, so $g = 1$).
Hence $C_3 = \frac{\sqrt\pi}{16}\int A^{-5/2}$ over pairs $(y_2, y_3)$ with $|y_2|, |y_3|, |y_2 - y_3| \leq 1$. With
$a = |y_2|$, $b = |y_3|$ and $\theta \in [0, \pi]$ the angle between them ($dy_2\,dy_3 = 4\pi\, a\, b\, da\, db\, d\theta$) and
$G(A) = \int_0^\infty r^4 e^{-A r^2} dr = 3\sqrt\pi / (8 A^{5/2})$:

$$C_3 = 4\pi \int_0^1\!\!\int_0^1\!\!\int_0^\pi \tfrac16\,\mathbf 1[c \leq 1]\; G(A(a,b,c))\; a\, b\; d\theta\, db\, da, \qquad c^2 = a^2 + b^2 - 2ab\cos\theta .$$

## How each cell is bounded

The cube $[0,1]^2 \times [0,\pi]$ is cut into cells. On a cell, the largest $c$ is computed from rational bounds of
$\cos\theta$ at the cell's corners. If it exceeds 1 the cell may contain paths, and it counts as **zero**. Otherwise every
configuration in it is a triangle, and the cell is bounded below by its weight $\frac{a_1^2 - a_0^2}{2}\cdot\frac{b_1^2 - b_0^2}{2}\cdot(\theta_1 - \theta_0)$
times $\frac16 G(\bar A)$. Here $\bar A$ is the exact union area at $(a_1, b_1, \bar c)$, with $\bar c$ the cell's largest $c$
rounded up to the grid $1/n_c$. If that triple is not a triangle, it is expanded until it is one.

By the planar Kneser-Poulsen theorem (K. Bezdek and R. Connelly, J. reine angew. Math. 553 (2002)), the area of a union
of unit discs does not decrease when no pairwise distance decreases. So $A \leq \bar A$ on the whole cell, and $G$
decreases in $A$.

## Two programs that share no code

    python3 -I tsp_c3_a.py --na 120 --nt 300 --nc 160 --claim 0.6425 --controls
    python3 -I tsp_c3_b.py --na 120 --nt 300 --nc 160 --claim 0.6425 --controls

| program | arithmetic | union of three unit discs | C_3 at 120/300/160 | run time |
| --- | --- | --- | --- | --- |
| `tsp_c3_a.py` (with `lower_pairs_a.py`, `lower_triples_a.py`) | integer intervals | Green's theorem over the uncovered arcs | 0.011490649266 | about 31 min (load 60-120) |
| `tsp_c3_b.py` (with `lower_pairs_b.py`, `lower_triples_b.py`) | 40-digit decimals, outward 1e-30 margins | inclusion-exclusion, triple intersection from arcs | 0.011490649265622 | about 13 min |

Both programs use only the Python standard library. Their arithmetic and their union-area routines share no code; the
discretisation (cells, corner rule, rounding of $c$ to the grid) is the same design in both. On the coarse grids
12/24/12 and 24/48/24 they agree to the 12th digit (0.00953280475, 0.01057339274), and on the final grid program A's
exact rational lower bound 0.01149064926562172913... agrees with program B's 0.011490649265622 in every printed digit.
A third program written separately inside the lab (interval arithmetic with mpmath, its own union routine by polar
integration over Voronoi sectors) reproduces the cell sums exactly on the grids 24/48/24 and 60/120/80 and on the first
16 rows of the final grid.

Red controls (both programs):
- the claim moved up past the certified value is refused;
- a coarse grid gives a smaller value;
- counting the cells that mix triangles and paths as triangles raises the value, so the zero weight is in effect.

Outputs: `results/`.

The four `lower_*` files are the programs of our certified lower bound for the minimum spanning tree constant (the
isolated-triples term), reused unchanged. The C_3 programs take the triple term with cut 0 and with the path weight 0.
Program A computes $G$ in closed form ($3\sqrt\pi / (8A^{5/2})$ in interval arithmetic), because the MST routine's
series at cut 0 is very slow.

## Not claimed

- No peer review yet.
- Dong and Cao's preprint has not been refereed either. This bound inherits their Theorem 3.11 and their certified
  $C_4$ and $C_5$.
- Non-rigorous estimates of the exact $C_3$: about 0.011715 (adaptive quadrature in the angle with Gauss-Legendre in
  $a$, $b$; Monte Carlo 0.011725 $\pm$ 0.000015). The certified value is 98.1% of it; the ceiling of this route with
  Dong and Cao's $C_4$, $C_5$ is about 0.64277. Their importance-sampling estimator, reimplemented, gives
  $C_4 \approx 0.00639$ and $C_5 \approx 0.00392$ (their certified 0.00504 and 0.00102): the larger room is there.

Prepared with AI assistance (Claude, Anthropic).
