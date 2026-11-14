# Five elementary rows on four pages without a replay folder: 20b, 25a, 39a, 11a (14 November 2026)

These rows rest on an elementary object named in the row itself (teorth/optimizationproblems at main `87b4ebb`; the
four pages are pinned by `fetch_page.sh`, but the programs do not read them). They are exact and short; we list them
apart from the certificate replays, so that they can be counted apart.

| row (page text) | object | program A (`celem_a.py`) | program B (`celem_b.py`) | result |
|---|---|---|---|---|
| 20b, "$C_{20b}\ge 1/\sqrt{12}$, Classical, achieved by the cube $[-\frac12,\frac12]^n$ (volume 1)" | the cube | moments from the one-dimensional integrals, Gaussian elimination | moments by iterated integration of monomials over the box, Bareiss elimination | $\det\mathrm{Cov}/\mathrm{Vol}^2 = 12^{-n}$ for $n \le 12$: $L = 1/\sqrt{12}$ in every dimension checked |
| 20b, "$C_{20b}\ge 1/e$, Classical, achieved asymptotically by simplices" | the simplices $\Delta_n$ | Dirichlet moments $\prod a_i!/(n+\lvert a\rvert)!$ | moments by peeling coordinates, Beta integrals by binomial expansion | $L_{\Delta_n}^{2n} = (n!)^2/((n+1)^{n+1}(n+2)^n)$ exactly for $n \le 24$; $L_n$ increasing and below $1/e$ there ($L_{24} = 0.359569\ldots$) |
| 25a, "$4$, Trivial, the cube $B_\infty^n$: $\mathrm{Vol} = 2^n$, $\mathrm{Vol}((B_\infty^n)^\circ) = 2^n/n!$" | the cube and the cross-polytope | $\mathrm{Vol}(K^\circ)$ from its $2^n$ orthant simplices | $\mathrm{Vol}(K^\circ)$ as the leading coefficient of its Ehrhart polynomial (lattice-point counts of $tB_1^n$, $t \le n$, by a recursion; Lagrange interpolation) | $n!\,\mathrm{Vol}(K)\mathrm{Vol}(K^\circ) = 4^n$ for $n \le 12$ |
| 39a, "$8$, Classical (cube)" | the cube $[-1,1]^3$ | any two vertices differ by 2 in a coordinate | for any two vertices the sets of admissible translations are disjoint (a projection does not overlap) | a translate of the interior holds at most one of the 8 vertices: $H_3 \ge 8$ |
| 11a, "$1$, Trivial: for $n = 1$, take $f(x) = x$" | $f(x) = x$ on $\{-1, 1\}$ | the two expectations | the same by enumeration of $\{-1,1\}^n$ | $\mathbb E\lvert f - \mathbb E f\rvert = 1 = \mathbb E\lvert\nabla f\rvert$: ratio 1 |

**Not checked (mathematics, read):** for the simplices, the limit $L_{\Delta_n} \to 1/e$ (Stirling's formula), so
that $C_{20b} = \sup_n L_n^{\mathrm{lc}} \ge 1/e$; the identity of the polar of the cube with the cross-polytope
(program A checks that $\pm e_i$ lie on the boundary of $K^\circ$); for 39a, that $H_3$ is at least the number of
points no two of which fit in one translate of the interior (the definition). The printed decimals
$0.367879$ and $0.288675$ are truncations of $1/e$ and $1/\sqrt{12}$, so valid as lower bounds.

**Red controls**, each refused by both programs: the 20b cube with side 2 and its volume left at 1; one second moment
of the 5-simplex raised by $1/1000$; the cube taken as its own polar (25a); the 16 vertices of a 4-cube projected to
three coordinates (repeated points, 39a); $D_j f$ taken without its factor $1/2$ (11a).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; about ten seconds); `./fetch_page.sh` pins the pages.
Output: `results/celem.txt`.
