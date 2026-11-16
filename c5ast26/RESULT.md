# 5a: the bound of Carter, Hunter and O'Bryant maximised exactly, [St26] $C_{5a} \le 0.98173$\* and [CHO25] $C_{5a} \le 0.98183$ (16 November 2026)

Both rows rest on the same published framework: Carter, Hunter and O'Bryant, *On the diameter of finite Sidon sets*
(arXiv 2310.20032, Acta Math. Hungar. 175 (2025), 108-126), Section 3, read as text. For levels
$0 < \alpha_1 < \dots < \alpha_6$, a scale $\tau$ and three window ratios $c_1, c_2, c_3$, Lemma 3.1 and three instances
of Lemma 3.2 give, for every $w$ with $0 = w_0 \le w_1 \le \dots \le w_6 \le w_7 = 1$, four upper bounds on
$b_\infty$; so $b_\infty \le b^* := \max_w \min(f_1, f_{c_1}, f_{c_2}, f_{c_3})(w)$, and $C_{5a} \le b^*/2$
(Proposition 1.3). Here

- $f_1(w) = \tau + \frac1\tau - 2\tau\sum_{j=0}^{6}(w_{j+1}-w_j)\min_{\alpha_j\le z\le\alpha_{j+1}}(z-1)^2$,
- $f_c(w) = c\tau + \frac1{c\tau} - \frac{2\tau}{c^2}\sum_{j=1}^{L}(p_j-p_{j-1})\min_{\alpha_{\eta_j-1}-\alpha_{\zeta_j}\le z\le\alpha_{\eta_j}-\alpha_{\zeta_j-1}}(z-c)^2$,

with $p_0 < \dots < p_L$ the points $w_j$ and $w_j + c$ in $[0, 1]$, $\zeta_j$ and $\eta_j$ the least indices with
$w_{\zeta_j} + c \ge p_j$ and $w_{\eta_j} \ge p_j$, $\alpha_{-1} = \alpha_0 = 0$ and $\alpha_7 = \infty$.

- **[St26]** (Stoltz; rynald0cst0ltziam/Erdos-Sidon-Records-30 @ `aa40ce9`, the certified record
  `records/certified/record_1771350836_46537e4f6f8d/certified_record.json`, read as data): $\tau = 1.128028$,
  $\alpha = (0.708402, 0.785148, 0.869892, 1.124418, 1.182121, 1.243221)$, $c = (0.664209, 0.664134, 0.855601)$.
  Its 19,736 LP certificates were neither read nor used, and none of its programs was read or run: this replay
  recomputes the bound from the parameters alone.
- **[CHO25]**: the parameters printed in the proof of Theorem 3.3.

Two separately written exact programs (Python `fractions`, standard library only):

- `c5a_st26_a.py` (A): cells on which every difference $w_b - w_a$ lies strictly between two consecutive thresholds
  of $\{0, c_1, c_2, c_3\}$, enumerated depth-first with an incremental Floyd-Warshall test of the difference
  constraints in integers; the four bounds on a cell by counting the points $w_j$ and $w_j + c$ in merged order; on
  each cell $\max t$ subject to $t \le f_k(w)$ by a two-phase primal simplex (Bland's rule).
- `c5a_st26_b.py` (B): cells as the paper builds them, from the interlacings $s$ of $q_j = w_j + c$ and $r_j = w_j$,
  chosen jointly for the three windows (warm-started Bellman-Ford); $f_c$ from the paper's $p$, $\zeta$, $\eta$; on
  each cell the dual linear program $\min \sum_k \lambda_k F_k + \sum_i \mu_i h_i$ subject to $\sum_k \lambda_k = 1$,
  $\sum_k \lambda_k a_k = \sum_i \mu_i g_i$, $\lambda, \mu \ge 0$.
- Each skips a cell only on a valid upper bound of its max-min below the best value found so far (A: a box from the
  shortest-path distances; B: the ranges of the increments $w_{j+1} - w_j$ with $\sum = 1$; both also with the dual
  weights of the best cell). Self-tests: A checks its cell formulas against the definition at 200 random points per
  window; B checks its cell functions at an interior point of 300 cells.

| | [St26] | [CHO25] |
|---|---|---|
| full-dimensional cells (A and B, two enumerations) | 19,658 | 24,736 |
| $b^*$ (the same fraction in A and B) | $1.96345059371304719729\ldots$ | $1.96364548408134071674\ldots$ |
| attained at $w_1, \dots, w_6$ | $0.0943871961, 0.335791, 0.3880977121, 0.8328981937, 1, 1$ | $0.1339827349, 0.3001547028, 0.4621967345, 0.9647647028, 0.9779547028, 1$ |
| the four bounds there | all equal to $b^*$ | all equal to $b^*$ |
| dual weights of $(f_1, f_{c_1}, f_{c_2}, f_{c_3})$ | $0.78752, 0.20545, 0.00703, 0.0000071$ | $0.78217, 0.11247, 0.10536, 0.00000016$ |
| $b^*/2$ | $0.98172529685652359864\ldots \le 0.98173$: **holds** | $0.98182274204067035837\ldots \le 0.98183$: **holds** |

The exact fractions are in `results/c5ast26.txt`. St26's record states $b_\infty = 1.9634505937130473$ and
$0.9817252968565237$ in floating point, and its witness $w$ is the maximiser above. For CHO25 the maximiser is the
worst case the paper prints ($w \approx (0.13398, 0.30015, 0.46220, 0.96476, 0.97795, 1)$, all four bounds
$1.963645$); the paper counts 24,822 nonempty triples of cells, which presumably includes lower-dimensional
intersections (both our enumerations require nonempty interiors). With the single window $c_3 = 0.855601$ (red
control 3) program A finds 127 cells, the number the paper gives for each instance of Lemma 3.2 with its own $c$'s.

**Red controls**, each refused by both programs: the St26 claim lowered to 0.98172; Lemma 3.2 with $\tau/c^2$ in place
of $2\tau/c^2$ ($b^* = 1.98406\ldots$); the two near-equal windows $c_1 = 0.664209$ and $c_2 = 0.664134$ removed ($b^* =
1.99584\ldots$). (Removing $c_3$ alone was not used as a control: its dual weight at the maximum is $7 \cdot 10^{-6}$.)

**Not checked (mathematics, read):** Theorem 1.1 (the Erdős-Turán Sidon set equality of [OBr22]), Corollary 1.2,
Lemmas 3.1 and 3.2 and Proposition 1.3, which turn $b^*$ into the bound; they are in the refereed paper.

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; one process at a time, about 35 minutes). Output:
`results/c5ast26.txt`.
