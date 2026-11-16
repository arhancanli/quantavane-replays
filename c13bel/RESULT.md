# 13b: the lower bound [Elek1994] 0.8257 (16 November 2026)

A convex set covering a congruent copy of every set of diameter 1 contains a disk $D$ of diameter 1 and an
equilateral triangle $T$ of side 1, so its area is at least $\min_v m(v)$, $m(v) = \operatorname{area}\operatorname{conv}(D \cup (T+v))$.
Two separately written programs show $\min_v m(v) \ge 0.8257$:

- `c13bel_a.py` (the argument, in exact rational intervals): a vertex $p$ at distance $d > 1/2$ adds a cap of area
  $g(d) = \frac12\sqrt{d^2 - \frac14} - \frac14\arccos\frac{1}{2d}$, inside the angular sector of half-angle
  $\arccos\frac{1}{2d}$; $g$ is convex and non-decreasing, so $L(v) = \pi/4 + \sum_k g(\lvert q_k + v\rvert)$ is convex
  and, by the $120°$ symmetry, $L(v) \ge L(0)$. If all vertices are within $d_1 = 7/10 \le 1/\sqrt2$ of the centre,
  the caps' sectors are disjoint and $m(v) \ge L(v) \ge L(0)$; otherwise $m(v) \ge \pi/4 + g(7/10) = 0.83654\ldots$.
  $L(0) = \pi/8 + \sqrt3/4 = 0.825711783590943478\ldots$ (attained with the centres together).
- `c13bel_b.py` (no convexity): branch and bound over $v \in [-\frac12, \frac12]^2$ (99,322 cells, depth up to 16),
  each cell certified $\ge 0.8257$ by the single-cap bound or by the disjoint-caps bound, the sectors' disjointness
  checked on each cell; $\lvert v\rvert \ge 1/2$ by a far-field bound (at least 0.8598).

The printed 0.8257 is a truncation of $\pi/8 + \sqrt3/4$. Red controls, refused by each: the claim 0.8258; $d_1 = 3/4$
and a triangle of side 0.9 (A); a disk of radius 0.49 and the disjoint-caps bound removed (B). (B's radius control
first passed, which exposed that the disk's area had been written as $\pi/4$ for every radius; fixed before any result
was kept.)

Reproduce: `./run_all.sh` (Python 3.9+; about two minutes). Output: `results/c13bel.txt`.
