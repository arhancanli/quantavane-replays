# Euclidean TSP (Beardwood-Halton-Hammersley) constant: beta >= 0.64309, sharper terms in Dong and Cao's series

Constant 12 of the [optimization constants repository](https://teorth.github.io/optimizationproblems/constants/12a.html).
The length of the shortest tour through n uniform random points of the unit square is asymptotic to beta sqrt(n).

## The result

Dong and Cao ([arXiv:2610.02809](https://arxiv.org/abs/2610.02809), Theorem 3.11) prove, for every K >= 3,

$$\beta \geq \frac58 + \sum_{k=3}^{K} C_k,\qquad C_k = \frac{\Gamma(k-\frac12)}{2\,k!}\int_{(\mathbb R^2)^{k-1}} \frac{g(G_1(y))}{A(y)^{k-1/2}}\,dy_2\cdots dy_k .$$

The symbols are:
- $y_1 = 0$;
- $G_1(y)$ is the graph joining two points at distance at most 1;
- $g = q - D$ is their lower bound for the excess of a connected component's path-cover number over its degree count;
- $A(y)$ is the area of the union of the unit discs at the $y_i$.

Each $C_k$ is bounded from below separately. Dong and Cao certify $C_3 \geq 0.01113528859$, $C_4 \geq 0.005040573276$ and
$C_5 \geq 0.001015487669$, giving $\beta \geq 0.6421913$. We certify sharper values of the first two:

| term | Dong-Cao | ours (program A / program B) | method |
| --- | --- | --- | --- |
| $C_3$ | 0.01113528859 | 0.0114906492656 / 0.011490649265622 | exact union areas on (a, b, angle) cells |
| $C_4$ | 0.005040573276 | 0.005584568980 / 0.005584568975 | DC's cubes, adaptive to side 1/1024, exact three-disc areas |
| $C_5$ | 0.001015487669 | (theirs is used; see below) | |

Using the smaller value of each pair:

$$\beta \geq \frac58 + 0.0114906492 + 0.0055845689 + 0.0010154876 = 0.64309070\ldots \geq 0.64309,$$

with Dong and Cao's $C_5$ unchanged.

This is $8.99 \cdot 10^{-4}$ above Dong and Cao's 0.6421913; the repository's page lists 0.6277. The lower-bound
formula, the cube method and the value of $C_5$ are Dong and Cao's. What is new is the evaluation of $C_3$ and $C_4$:
exact areas of unions of three discs, planar Kneser-Poulsen monotonicity, and adaptive refinement. Each of the two
improves theirs on its own. Our own $C_5$ computation (7 dimensions, cubes down to side 1/16) reached only 0.000945
at the resolution we could afford today, below their 0.0010155, so theirs is used.

## C_3: a three-dimensional integral

For three points, a connected $G_1$ is either a path ($q = D = 1$, so $g = 0$) or a triangle ($g = 1$). Hence
$C_3 = \frac{\sqrt\pi}{16}\int A^{-5/2}$ over pairs $(y_2, y_3)$ with all three distances at most 1. In the coordinates
$a = |y_2|$, $b = |y_3|$ and $\theta$ the angle between them, $dy_2\,dy_3 = 4\pi\,a\,b\,da\,db\,d\theta$.

Each cell of $[0,1]^2 \times [0,\pi]$ is handled as follows:
- if its largest $c = |y_2 - y_3|$ exceeds 1, it counts as zero;
- otherwise it is bounded below by its weight times $\frac16 G(\bar A)$, where $G(A) = 3\sqrt\pi/(8A^{5/2})$ and
  $\bar A$ is the exact union area at the cell's largest triangle.

The union area never decreases when no pairwise distance decreases. This is the planar Kneser-Poulsen theorem
(K. Bezdek and R. Connelly, J. reine angew. Math. 553 (2002)).

Programs: `tsp_c3_a.py` (integer intervals, Green's theorem) and `tsp_c3_b.py` (40-digit decimals,
inclusion-exclusion). On the grid 120/300/160 they agree to 12 digits.

## C_4: Dong and Cao's cubes with sharper areas

The integrals use Dong and Cao's coordinates (their (20)) on $|y_j| \leq k-2$ (their (25)). Each cube is handled as
follows:
- exact integer bounds of every squared distance give the edges that are present on the whole cube, absent on the
  whole cube, or undetermined;
- the cube's $g$ is the minimum of $g$ over every graph compatible with them;
- the area is bounded above by the smaller of two forms, using the tables:
  - $A_3(T) + c(d_{lm})$;
  - $A_3(T) + A_3(l, m_1, m_2) - A_2(m_1, m_2)$.

Here $T$ is any triple of the four points and $l$ the remaining point. $A_3$ is the union area of three unit discs at the rounded-up largest distances, by Kneser-Poulsen.
$c(d) = \pi - \text{lens}(d)$ is the part of one disc outside another. $A_2$ is the union area of two discs, at the
rounded-down smallest distance. The second form uses
$|B_l \setminus (B_{m_1}\cup B_{m_2})| = |B_l\cup B_{m_1}\cup B_{m_2}| - |B_{m_1}\cup B_{m_2}|$.

$A_3$, $c$ and $A_2$ come from a table at spacing 1/64 (`c4_tables.py`):
- $A_3$ and $c$ are rounded up at 20 digits;
- $A_2$ is rounded down.

The table is computed twice, once with integer intervals and Green's theorem and once with 40-digit decimals and
inclusion-exclusion. The two files are byte-identical. It is stored as `tables/tables_64.txt.gz` (SHA-256 of the
uncompressed file `2261b3183f538e18f15eafb5017f50ea81a69f218c21401bf51b6f5296b3e320`) and `tables/tables_64_L.txt`.

Cubes are split adaptively. A heuristic estimate of the possible gain decides where to split, and validity does not
depend on those decisions.

The two implementations, `c4_a.c` and `c4_b.c`, are written separately from each other:
- program A uses recursion, interval-of-differences distances and subset enumeration;
- program B uses an explicit stack, clamped box distances, union-find for $g$, a scan of all edge masks and
  bisection for the table indices.

The arithmetic is IEEE double with outward margins:
- table values are moved one ulp outward;
- every area sum is moved up;
- `pow` is lowered by $10^{-13}$ relative;
- the totals are lowered by $10^{-12}$ (A) or $10^{-9}$ (B) relative.

Commands:

    gunzip -k tables/tables_64.txt.gz
    cc -O2 -o c4_a c4_a.c -lm && cc -O2 -o c4_b c4_b.c -lm
    C4_LTAB=tables/tables_64_L.txt C4_DEPTH=7 C4_TAU=3e-12 ./c4_a tables/tables_64.txt 8 0.0055
    C4_LTAB=tables/tables_64_L.txt C4_DEPTH=7 C4_TAU=3e-12 ./c4_b tables/tables_64.txt 8 0.0055

Run times on one core at load about 60: about 1 h 50 min each. Both programs report 7,001,145,723 positive cubes and
299,028,686 splits, and $C_4 \geq 0.005584568980$ (A), $0.005584568975$ (B); the difference is B's larger final
safety margin. An earlier, coarser run (depth 6, threshold $10^{-11}$: $C_4 \geq 0.005450602910$ / $0.005450602904$,
2,827,471,112 positive cubes) is kept in `results/` as well; the referee's sampling checks were made on that version
of the code, which is unchanged (only the depth and threshold differ).

Outputs are in `results/`.

## Checks inside the lab

Two fresh referees worked separately from the author.

**$C_3$:**
- They re-derived the reduction, the constant and the change of variables.
- An independent mpmath interval program with its own union routine reproduces the cell sums exactly on two full
  grids and on a slice of the final grid.

**$C_4$:**
- All 357,760 table entries were checked at 40 digits.
- 20,000 random cubes at depths 0 to 6 were run through both programs; in each one, 262 sampled configurations plus a
  local ascent were checked against the bound. There were no area, $g$ or weight violations, and A and B were equal on
  every cube.
- A float reimplementation of the cube rule matches A to 12 digits.

The referee's notes on $C_4$: the two programs share the (byte-identical) table and the same cube rule, so their
agreement checks the implementations and the referee's review checks the rule; the arithmetic is IEEE double with
outward margins and relies on `pow` being accurate to a few ulps (it is not interval arithmetic throughout).

## Not claimed, and what is left

- Nothing here is peer reviewed. Dong and Cao's preprint is not refereed either. This bound inherits their Theorem 3.11.
- Non-rigorous estimates of the exact terms are:
  - $C_3 \approx 0.011715$ (quadrature);
  - $C_4 \approx 0.00641$ (Monte Carlo);
  - $C_5 \approx 0.0039$ (Dong and Cao's importance sampling, reimplemented).

  There is room left, mostly in $C_5$ (Dong and Cao certify about a quarter of it) and in further terms.

Prepared with AI assistance (Claude, Anthropic).
