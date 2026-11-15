# 50a, the Quantum Max Cut approximation ratio: [GP19] 0.498 from its definition (15 November 2026)

A value recomputed from its definition, not a certificate. [GP19] (S. Gharibian, O. Parekh, *Almost optimal classical
approximation algorithms for a quantum generalization of Max-Cut*, arXiv 1909.08846v1, TeX `GP2019_arxiv_v1.tex` read
as text) proves Theorem 4: rounding the SDP vectors to Bloch vectors in $\mathbb R^r$ ($r = 3$ for the
anti-ferromagnetic Heisenberg model, the 50a problem) gives, on an edge whose SDP vectors have inner product $t$, the
ratio $(1 - F_r(t))/(1 - rt)$, where (the lemma of Briët, de Oliveira Filho and Vallentin, in their appendix)
$$F_r(t) = \frac{2}{r}\Big(\frac{\Gamma(\frac{r+1}{2})}{\Gamma(\frac r2)}\Big)^2 t\, {}_2F_1\big(\tfrac12, \tfrac12; \tfrac r2 + 1; t^2\big),$$
and the algorithm's ratio is the minimum over $t \in [-1, 1/r)$. The paper evaluates it in Mathematica on a grid of
step 0.01 (its appendix) and prints $0.878$, $0.649$, $0.498$ for $r = 1, 2, 3$. The page's lower bound for $C_{50}$
is the $r = 3$ value, 0.498.

Two separately written programs bound the minimum over the whole interval from below; nothing is downloaded and no
program of any source was read or run. $F_r$ is odd with positive Taylor coefficients, so increasing; on a cell
$[a, b] \subset [-1, 0]$ the ratio is at least $(1 - F_r(b))/(1 - ra)$, and for $t \in [0, 1/r)$ it is at least
$1 - F_r(1/r)$.

- `c50a_a.py`: $F_3(t) = \frac{8}{3\pi} t\, {}_2F_1(\frac12,\frac12;\frac52;t^2)$ by 600 terms of the series in exact
  rationals (rounded outward to $2^{-200}$), the tail bounded by $c_n \le \binom{2n}{n}4^{-n}\cdot\frac{3/4}{n^2}$ and
  $\binom{2n}{n}4^{-n} \le 1/\sqrt{3n}$; $\pi$ by Machin's formula; cells bisected until each clears the target (2,103
  cells, least lower bound $0.4987004\ldots$).
- `c50a_b.py`: the closed form from Gauss's contiguous relation (with ${}_2F_1(\frac12,\frac12;\frac12;z) = (1-z)^{-1/2}$
  and ${}_2F_1(\frac12,\frac12;\frac32;t^2) = \arcsin(t)/t$):
  $F_3(t) = \frac{2}{\pi t}\big(\sqrt{1-t^2} - (1 - 2t^2)\frac{\arcsin t}{t}\big)$, $F_3(-1) = -1$; decimal intervals
  at 50 digits, $\arcsin$ through $\arctan$ (argument halving and the alternating series), $\pi = 4\arctan 1$; 20,000
  equal cells on $[-1, 0]$ (least lower bound $0.4987482\ldots$ on $[-0.96600, -0.96595]$).

| statement | both programs |
|---|---|
| $\min_{t \in [-1, 1/3)} (1 - F_3(t))/(1 - 3t) \ge 0.4987$ | yes (for $t \ge 0$ the ratio is above $0.7137$) |
| the page's 0.498 | holds |

Information (floating point): the minimum is $0.49876742\ldots$ at $t \approx -0.965953$; the paper's 0.01 grid gives
$0.4987794$; with $r = 2$ the minimum is $0.64991475\ldots$ (the paper's 0.649).

**Anchor** (program B): with $r = 1$, $F_1(t) = \frac2\pi\arcsin t$ and the minimum is Goemans and Williamson's
$0.87856\ldots$: at least 0.8785 on 20,000 cells (for $t \ge 0$ the ratio is at least 1, since $\arcsin t \le \frac\pi2 t$),
and a claim of 0.8786 is refused.

**Red controls**, each refused by both programs: the minimum claimed to be at least 0.4988 (program A stops at a cell
near $t = -0.959$ where the ratio is about 0.4988; program B's least cell bound is 0.49875); $(\Gamma(2)/\Gamma(3/2))$
not squared, that is $8/(3\pi)$ replaced by $4/(3\sqrt\pi)$ (the minimum drops to about 0.471).

**Not checked:** the rounding analysis of Theorem 4 and the BFV lemma (read).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; about three minutes). Output: `results/c50a.txt`.
