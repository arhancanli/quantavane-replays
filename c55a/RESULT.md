# 55a, the acyclic chromatic index coefficient: the computational inputs of three upper bounds (15 November 2026)

Values recomputed from their definitions, not certificates. Nothing is downloaded and no program of any source was
read or run; the three papers' TeX sources were read as text.

**[KLS2026]** (L. Kirousis, J. Livieratos, A. Singh, arXiv 2602.14859v1, the only version; e-print SHA-256
`c6785ab3…8709`): $\chi'_a(G) \le \lceil 3.142(\Delta-1)\rceil + 1$. With $(2+\gamma)(\Delta-1)+1$ colours the
algorithm's witness forests are counted by $t_n$, the number of rooted unordered trees with $n$ nodes whose internal
nodes all have even outdegree $\ge 4$, weighted by at most $((1 - e^{-1/\gamma})^{3/4})^n$; Lemma 8 needs
$\limsup t_n^{1/n} = 1/\rho$ with $\rho \ge 0.6677$ (their Appendix: $0.6677 \le \rho \le 0.6678$, by Maple) and
then $\gamma = 1.142$ works because $(1 - e^{-1/1.142})^{3/4} < 0.6677$.

The generating function satisfies $T = \Phi(z; T(z), T(z^2), \ldots)$ with
$\Phi(z; w) = z[(e^{\sum_j w_j/j} + e^{\sum_j (-1)^j w_j/j})/2 - (w_1^2 + w_2)/2]$ (the paper's Lemma 9), and $\Phi$
has nonnegative coefficients in every $w_j$. So a **supersolution**, $\Phi(x; Y, W_2, W_3, \ldots) \le Y$ with
$W_j \ge T(x^j)$, gives $T(x) \le Y < \infty$ and $\rho \ge x$ (the monotone iteration from 0 stays below $Y$). The
bounds $W_j$ come from a chain of supersolutions at smaller points (Cauchy: $t_n \le T(v)/v^n$), and
$T(v^j) \le v^{j-1}T(v)$ at the first point ($T(v)/v$ increases).

- `c55a_a.py`: $t_n$, $n \le 200$, from a multiset table (multisets of $k$ trees with $m$ nodes, built one tree size
  at a time); exact rational intervals rounded outward to $2^{-300}$, $e^x$ by its Taylor series with a remainder
  bound; the chain $v = 3/10, 1/2, 3/5$ (bounds $T(v) \le 0.304, 0.547, 0.761$), then $x = 0.6677$.
- `c55a_b.py`: $t_n$, $n \le 100$, from the functional equation itself, coefficient by coefficient (every $t_n$ comes
  out an integer); decimal intervals at 60 digits (exp and ln of the decimal module widened by one unit in the last
  place); the chain $v = 0.25, 0.45, 0.62$, then $x = 0.6677$; $\gamma$ through logarithms.

| statement | program A | program B |
|---|---|---|
| $t_n$ | $1, 0, 0, 0, 1, 0, 1, 0, 2, 0, 3, 0, 6, 0, 10, 0, 20, \ldots$, $t_{99} = 317145363238384$ | the same |
| supersolution at $x = 0.6677$ | $Y_0 = 1.343718520791$, margin $2.97 \cdot 10^{-5}$ | $Y_0 = 1.343718346535$, margin $2.94 \cdot 10^{-5}$ |
| so $\rho \ge 0.6677$ (Lemma 8) | yes | yes |
| $(1 - e^{-1/1.142})^{3/4} < 0.6677$ | $(1-e^{-1/1.142})^3 = 0.1985753410\ldots < 0.6677^4 = 0.1987584059\ldots$ | $\frac34\ln(1-e^{-1/1.142}) \le -0.4041466747 < \ln 0.6677$ |

Information: with the truncated series (coefficientwise below $T$) the fixed-point equation has no solution at
$x = 0.66771$ (program A, a convexity argument), so $0.6677 \le \rho \le 0.66771$: the paper's lower bound holds with
about $5 \cdot 10^{-6}$ to spare. The least $\gamma$ with $(1 - e^{-1/\gamma})^{3/4} < 0.6677$ is $1.1414390\ldots$.

**[FLM2020]** (P. M. S. Fialho, B. N. B. de Lima, A. Procacci, arXiv 1912.04436): $a'(G) \le 3.569(\Delta - 1)$,
from $\rho_\gamma = \min_{x>0}\phi_E(x)/x < 1$ at $\gamma = 1.569$, $\phi_E(x) = \gamma^{-1} a^3 (x+1)^4/(1 - a^2(x+1)^2)$,
$a = 1 - e^{-1/\gamma}$. **[GKPT2017]** (I. Giotis, L. Kirousis, K. I. Psaromiligkos, D. M. Thilikos, arXiv 1407.5374):
$\lceil 3.74(\Delta-1)\rceil + 1$, from the same test with $(x+1)^6$ in place of $(x+1)^4$, at $\gamma = 1.74$. Both
programs exhibit a point $x$ with $\phi(x)/x < 1$ in interval arithmetic: $0.9979145732\ldots$ (FLM2020, $x \approx 0.2466$)
and $0.9806102671\ldots$ (GKPT2017, $x \approx 0.1749$). Information (floating point): the thresholds are
$\gamma = 1.5681601\ldots$ and $1.7309403\ldots$; at GKPT2017's $\gamma = 1.73095$ the minimum is $0.99997890$, as their
Section 6 prints.

**Red controls**, each refused by both programs: the supersolution at $x = 0.66771$ (above $\rho$; margin
$-2.99 \cdot 10^{-5}$); $\gamma = 1.141$; the tree class with outdegree 2 allowed (no supersolution even at
$v = 0.62$ in B, none at $0.6677$ in A).

**Not checked:** the papers' probabilistic analyses (the witness-forest arguments and the weight per internal node),
read only; the [EP2013] 4, [NPS2012] 9.62 and earlier rows.

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; about five minutes). Output: `results/c55a.txt`.
